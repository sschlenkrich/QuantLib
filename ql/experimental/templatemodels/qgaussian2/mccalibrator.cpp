/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2019 Sebastian Schlenkrich

*/


#include <ql/experimental/templatemodels/qgaussian2/mccalibrator.hpp>

#include <ql/math/optimization/constraint.hpp>
#include <ql/math/optimization/levenbergmarquardt.hpp>
#include <ql/pricingengines/blackformula.hpp>


namespace QuantLib {


    QGMonteCarloCalibrator::Objective::CalibSwaption::CalibSwaption (
        Date                                expiryDate,
        const ext::shared_ptr<SwapIndex>&   swapindex,
        const Handle<YieldTermStructure>&   discountCurve,
        const Handle<SwaptionVolatilityStructure> volTS,
        bool                                contTenorSpread)
        : SwapCashFlows(swapindex->underlyingSwap(expiryDate),discountCurve,contTenorSpread) {
        expiryTime_ = Actual365Fixed().yearFraction(discountCurve->referenceDate(), expiryDate); 
        //S0_ = swapindex->fixing(expiryDate, true);
        // we use cash flows to calculate swap rate and annuity to be fully consistent to MC calculations
        annuity_ = 0.0;
        for (size_t k = 0; k < this->fixedTimes().size(); ++k) annuity_ += this->annuityWeights()[k] * discountCurve->discount(this->fixedTimes()[k]);
        S0_ = 0.0;
        for (size_t k = 0; k < this->floatTimes().size(); ++k) S0_ += this->floatWeights()[k] * discountCurve->discount(this->floatTimes()[k]);
        S0_ /= annuity_;  // S0 = FloatLeg / Annuity
        // now we can calculate vols, prices and vegas
        sigmaATM_ = volTS->volatility(expiryDate, swapindex->tenor(), S0_, true);
        Real sigmap1Std = volTS->volatility(expiryDate, swapindex->tenor(), S0_ + sigmaATM_*sqrt(expiryTime_), true);
        Real sigmam1Std = volTS->volatility(expiryDate, swapindex->tenor(), S0_ - sigmaATM_*sqrt(expiryTime_), true);
        // calculate prices... 
        atmCall_  = bachelierBlackFormula(Option::Call, S0_, S0_, sigmaATM_*sqrt(expiryTime_));
        highCall_ = bachelierBlackFormula(Option::Call, S0_ + sigmaATM_*sqrt(expiryTime_), S0_, sigmap1Std*sqrt(expiryTime_));
        lowPut_   = bachelierBlackFormula(Option::Put,  S0_ - sigmaATM_*sqrt(expiryTime_), S0_, sigmam1Std*sqrt(expiryTime_));
        // ... and vegas
        atmVega_  = bachelierBlackFormulaStdDevDerivative(S0_, S0_, sigmaATM_*sqrt(expiryTime_)) * sqrt(expiryTime_);
        highVega_ = bachelierBlackFormulaStdDevDerivative(S0_ + sigmaATM_*sqrt(expiryTime_), S0_, sigmap1Std*sqrt(expiryTime_)) * sqrt(expiryTime_);
        lowVega_  = bachelierBlackFormulaStdDevDerivative(S0_ - sigmaATM_*sqrt(expiryTime_), S0_, sigmam1Std*sqrt(expiryTime_)) * sqrt(expiryTime_);
    }


    QGMonteCarloCalibrator::Objective::Objective(
                  QGMonteCarloCalibrator                    *calibrator,
                  const std::vector< std::vector< Real > >&  isInput,
                  const std::vector< std::vector< Real > >&  isOutput ) 
        : calibrator_(calibrator), isInput_(isInput), isOutput_(isOutput) {
        // checking dimensions
        QL_REQUIRE(isInput_.size()==calibrator_->model_->times().size(),"QGCalibrator::Objective: wrong input dimension.");
        for (size_t i=0; i<isInput_.size(); ++i) {
            QL_REQUIRE(isInput_[i].size()==(3*(calibrator_->model_->factors()-1)),"QGCalibrator::Objective: wrong input dimension.");
        }
        QL_REQUIRE(isOutput_.size()== calibrator_->model_->times().size(),"QGCalibrator::Objective: wrong output dimension.");
        for (size_t i=0; i<isOutput_.size(); ++i) {
            QL_REQUIRE(isOutput_[i].size()==(3*calibrator_->swapIndices_.size()),"QGCalibrator::Objective: wrong output dimension.");
        }
        // count inputs and outputs
        inputSize_=0;
        for (size_t i=0; i<isInput_.size(); ++i)
            for (size_t j=0; j<isInput_[i].size(); ++j) if (isInput_[i][j]>0.0) ++inputSize_;
        outputSize_=0;
        for (size_t i=0; i<isOutput_.size(); ++i)
            for (size_t j=0; j<isOutput_[i].size(); ++j) if (isOutput_[i][j]>0.0) ++outputSize_;
        QL_REQUIRE(inputSize_>0,"QGCalibrator::Objective: inputSize_>0 required");
        QL_REQUIRE(outputSize_>0, "QGCalibrator::Objective: outputSize_>0 required");

        // set up target swaptions
        // we may have a set of calibration swaptions for each model time
        calibSwaptions_.resize(calibrator_->model_->times().size());
        for (size_t i=0; i<calibSwaptions_.size(); ++i) {
            calibSwaptions_[i].resize(calibrator_->swapIndices_.size());
            // find corresponding expiry date
            BigInteger nMonths = (BigInteger)(calibrator_->model_->times()[i] * 12.0);
            Calendar cal = calibrator_->swapIndices_[0]->fixingCalendar();  // assume all swap indices use thame calendar
            Date expiryDate = cal.advance(calibrator_->model_->termStructure()->referenceDate(), Period(nMonths, Months), Preceding);
            for (size_t j=0; j<calibSwaptions_[i].size(); ++j) {
                // allocate only if at least one of sigmaATM, skew or smile is required for calibration
                if (isOutput_[i][j]>0.0                                       ||
                    isOutput_[i][j+     calibrator_->swapIndices_.size()]>0.0 ||
                    isOutput_[i][j+ 2 * calibrator_->swapIndices_.size()]>0.0 ) {
                    calibSwaptions_[i][j] = ext::shared_ptr<CalibSwaption>(
                        new CalibSwaption(expiryDate, calibrator_->swapIndices_[j], calibrator_->model_->termStructure(), calibrator_->volTS_, true));
                }
            }
        }

        // we need to identify the simulation indices for this optimisation problem
        // first index is determined by independent variables
        size_t firstModelIdx = calibrator_->model_->times().size();
        for (size_t i=isInput.size(); i>0 ; --i) {
            for (size_t j=0; j<isInput[i-1].size(); ++j) {
                if (isInput[i-1][j]>0.0) firstModelIdx = i-1;
                break;  // one input suffices
            }
        }
        // Now we need to find the smallest simulation time s.t. T[firstModelIdx-1]<T
        firstSimulationIdx_ = calibrator_->mcSimulation_->simTimes().size();  // this should be overwritten
        if (firstModelIdx == 0) firstSimulationIdx_ = 1;
        else {
            for (size_t k = 1; k < calibrator_->mcSimulation_->simTimes().size(); ++k) {
                if (calibrator_->mcSimulation_->simTimes()[k] > calibrator_->model_->times()[firstModelIdx - 1]) {
                    firstSimulationIdx_ = k;
                    break;  // we want to have the minimum
                }
            }
        }
        QL_REQUIRE(firstSimulationIdx_ < calibrator_->mcSimulation_->simTimes().size(),
            "QGCalibrator::Objective: firstSimulationIdx_<calibrator_->mcSimulation_->simTimes().size() required");

        // last index is determined by calibration targets
        size_t lastModelIdx = calibrator_->model_->times().size();
        for (size_t i = 0; i < isOutput_.size(); ++i) {
            for (size_t j=0; j<isOutput_[i].size(); ++j)
                if (isOutput_[i][j] > 0.0) {
                    lastModelIdx = i;
                    break; // one output suffices
                }
        }
        // Now we need to find the largest simulation time s.t. T <= T[lastModelIdx]
        lastSimulationIdx_ = calibrator_->mcSimulation_->simTimes().size();  // this should be overwritten
        for (size_t k = calibrator_->mcSimulation_->simTimes().size(); k > 0; --k) {
            if (calibrator_->mcSimulation_->simTimes()[k - 1] <= calibrator_->model_->times()[lastModelIdx]) {
                lastSimulationIdx_ = k - 1;
                break;  // we want to have the maximum
            }
        }
        QL_REQUIRE(lastSimulationIdx_ < calibrator_->mcSimulation_->simTimes().size(),
            "QGCalibrator::Objective: lastSimulationIdx_<calibrator_->mcSimulation_->simTimes().size() required");

        QL_REQUIRE(firstSimulationIdx_<=lastSimulationIdx_, "QGCalibrator::Objective: firstSimulationIdx_<=lastSimulationIdx_ required");

        // ready for optimisation...
    }

    Array QGMonteCarloCalibrator::Objective::initialise() {
        Array X(inputSize_);
        size_t d = calibrator_->model_->factors()-1;
        size_t idx = 0;
        for (size_t i=0; i<isInput_.size(); ++i) {
            for (size_t j=0; j<d; ++j) {
                if (isInput_[i][j]>0.0) {
                    Real sigma = calibrator_->model_->sigma()[j][i];
                    Real x = QGMonteCarloCalibrator::inverse(sigma,calibrator_->sigmaMin_,calibrator_->sigmaMax_) / isInput_[i][j];
                    X[idx] = x;
                    ++idx;
                }
            }
            for (size_t j=0; j<d; ++j) {
                if (isInput_[i][d+j]>0.0) {
                    Real slope = calibrator_->model_->slope()[j][i];
                    Real x = QGMonteCarloCalibrator::inverse(slope,calibrator_->slopeMin_,calibrator_->slopeMax_) / isInput_[i][d+j];
                    X[idx] = x;
                    ++idx;
                }
            }
            for (size_t j=0; j<d; ++j) {
                if (isInput_[i][2*d+j]>0.0) {
                    Real curve = calibrator_->model_->curve()[j][i];
                    Real x = QGMonteCarloCalibrator::inverse(curve,calibrator_->curveMin_,calibrator_->curveMax_) / isInput_[i][2*d+j];
                    X[idx] = x;
                    ++idx;
                }
            }
        }
        return X;
    }

    void QGMonteCarloCalibrator::Objective::update(const Array& X) const {
        std::vector< std::vector< Real > > m_sigma = calibrator_->model_->sigma();
        std::vector< std::vector< Real > > m_slope = calibrator_->model_->slope();
        std::vector< std::vector< Real > > m_curve = calibrator_->model_->curve();
        size_t d = calibrator_->model_->factors()-1;
        size_t idx = 0;
        for (size_t i=0; i<isInput_.size(); ++i) {
            for (size_t j=0; j<d; ++j) {
                if (isInput_[i][j]>0.0) {
                    Real x = X[idx];
                    Real sigma = QGMonteCarloCalibrator::direct(x*isInput_[i][j],calibrator_->sigmaMin_,calibrator_->sigmaMax_);
                    m_sigma[j][i] = sigma;
                    ++idx;
                }
            }
            for (size_t j=0; j<d; ++j) {
                if (isInput_[i][d+j]>0.0) {
                    Real x = X[idx];
                    Real slope = QGMonteCarloCalibrator::direct(x*isInput_[i][d+j],calibrator_->slopeMin_,calibrator_->slopeMax_);
                    m_slope[j][i] = slope;
                    ++idx;
                }
            }
            for (size_t j=0; j<d; ++j) {
                if (isInput_[i][2*d+j]>0.0) {
                    Real x = X[idx];
                    Real curve = QGMonteCarloCalibrator::direct(x*isInput_[i][2*d+j],calibrator_->curveMin_,calibrator_->curveMax_);
                    m_curve[j][i] = curve;
                    ++idx;
                }
            }
        }
        calibrator_->model_->update(m_sigma,m_slope, m_curve);  // Note: we work directly on the calibrator model
    }

    Disposable<Array> QGMonteCarloCalibrator::Objective::values(const Array& x) const {
        update(x);  // first we update the model with the new state
        for (size_t k = firstSimulationIdx_; k <= lastSimulationIdx_; ++k)  // then we simulate the model again... this is expensive!
            calibrator_->mcSimulation_->simulate(k,true);  // we clear additional states since we need to allow interpolation due to good fixing dates
        // MC sampling and objective calculation
        Array objective(outputSize_);
        size_t idx=0;
        for (size_t i=0; i<calibSwaptions_.size(); ++i) {
            for (size_t j=0; j<calibSwaptions_[i].size(); ++j) {
                Real atmCall=0.0, highCall=0.0, lowPut=0.0;
                if ((isOutput_[i][j] > 0.0) ||
                    (isOutput_[i][j + calibSwaptions_[i].size()] > 0.0) ||
                    (isOutput_[i][j + 2 * calibSwaptions_[i].size()] > 0.0)) {
                    // for convenience we pre-calculate some quantities
                    Real expiryTime = calibSwaptions_[i][j]->expiryTime();
                    Real S0 = calibSwaptions_[i][j]->S0();
                    Real annuity = calibSwaptions_[i][j]->annuity();
                    size_t nPaths = calibrator_->mcSimulation_->nPaths();
                    // gather swap rate samples etc. and calculate adjusters
                    std::vector<Real> oneOverBSample(nPaths, 0.0);
                    std::vector<Real> annuitySample(nPaths, 0.0);
                    std::vector<Real> swapRateSample(nPaths, 0.0);
                    ext::shared_ptr<QGMonteCarloCalibrator::MCPayoff> mcFloatLeg(new MCAnnuity(expiryTime, calibSwaptions_[i][j]->floatTimes(), calibSwaptions_[i][j]->floatWeights()));
                    ext::shared_ptr<QGMonteCarloCalibrator::MCPayoff> mcFixedLeg(new MCAnnuity(expiryTime, calibSwaptions_[i][j]->fixedTimes(), calibSwaptions_[i][j]->annuityWeights()));
                    ext::shared_ptr<QGMonteCarloCalibrator::MCPayoff> one(new QGMonteCarloCalibrator::MCBase::FixedAmount(1.0));
                    ext::shared_ptr<QGMonteCarloCalibrator::MCPayoff> oneAtT(new QGMonteCarloCalibrator::MCBase::Pay(one, expiryTime));
                    for (size_t k = 0; k < nPaths; ++k) {
                        ext::shared_ptr<MCSimulation::Path> p = calibrator_->mcSimulation_->path(k);
                        oneOverBSample[k] = oneAtT->discountedAt(p);
                        annuitySample[k] = mcFixedLeg->at(p);
                        swapRateSample[k] = mcFloatLeg->at(p) / annuitySample[k];
                    }
                    // calculate adjusters suksessively
                    Real mcDF = 0.0;
                    for (size_t k = 0; k < nPaths; ++k) mcDF += oneOverBSample[k];
                    mcDF /= nPaths;
                    Real adjOneOverB = calibrator_->model_->termStructure()->discount(expiryTime) / mcDF;
                    for (size_t k = 0; k < nPaths; ++k) oneOverBSample[k] *= adjOneOverB;
                    Real mcAnnuity = 0.0;
                    for (size_t k = 0; k < nPaths; ++k) mcAnnuity += (annuitySample[k] * oneOverBSample[k]);
                    mcAnnuity /= nPaths;
                    Real adjAnnuity = annuity / mcAnnuity;
                    for (size_t k = 0; k < nPaths; ++k) annuitySample[k] *= adjAnnuity;
                    Real mcFloat = 0.0;
                    for (size_t k = 0; k < nPaths; ++k) mcFloat += (annuitySample[k] * swapRateSample[k] * oneOverBSample[k]);
                    mcFloat /= nPaths;
                    Real adjSwapRate = S0 - mcFloat / annuity;
                    for (size_t k = 0; k < nPaths; ++k) swapRateSample[k] += adjSwapRate;
                    // now we can calculate the option prices
                    atmCall = 0.0;
                    highCall = 0.0;
                    lowPut = 0.0;
                    Real highStrike = S0 + calibSwaptions_[i][j]->sigmaATM()*sqrt(expiryTime);
                    Real lowStrike  = S0 - calibSwaptions_[i][j]->sigmaATM()*sqrt(expiryTime);
                    for (size_t k = 0; k < nPaths; ++k) {
                        Real radonNikodymDerivative = annuitySample[k] * oneOverBSample[k] / annuity;
                        if (swapRateSample[k] > S0)
                            atmCall += radonNikodymDerivative * (swapRateSample[k] - S0);
                        if (swapRateSample[k] > highStrike)
                            highCall += radonNikodymDerivative * (swapRateSample[k] - highStrike);
                        if (swapRateSample[k] < lowStrike)
                            lowPut += radonNikodymDerivative * (lowStrike - swapRateSample[k]);
                    }
                    atmCall /= nPaths;
                    highCall /= nPaths;
                    lowPut /= nPaths;
                }
                // TBD: calculate market prices and vegas!!
                if (isOutput_[i][j]>0.0) {
                    objective[idx] = (atmCall - calibSwaptions_[i][j]->atmCall()) / calibSwaptions_[i][j]->atmVega() * isOutput_[i][j];
                    ++idx;
                }
                if (isOutput_[i][j+calibSwaptions_[i].size()]>0.0) {
                    objective[idx] = ( (highCall-calibSwaptions_[i][j]->highCall())/calibSwaptions_[i][j]->highVega() -
                                       (lowPut-calibSwaptions_[i][j]->lowPut())/calibSwaptions_[i][j]->lowVega() )
                                     * isOutput_[i][j + calibSwaptions_[i].size()];
                    ++idx;
                }
                if (isOutput_[i][j+2*calibSwaptions_[i].size()]>0.0) {
                    objective[idx] = ( (highCall-calibSwaptions_[i][j]->highCall())/calibSwaptions_[i][j]->highVega() +
                                       (lowPut-calibSwaptions_[i][j]->lowPut())/calibSwaptions_[i][j]->lowVega()      - 
                                       2*(atmCall-calibSwaptions_[i][j]->atmCall())/calibSwaptions_[i][j]->atmVega()  )
                                     * isOutput_[i][j + 2*calibSwaptions_[i].size()];
                    ++idx;
                }
            }
        }
        // Logging
        // std::string logString("objective = ");
        // for (size_t k = 0; k < objective.size(); ++k) logString = logString + " " + std::to_string(objective[k]) + ",";
        // calibrator_->debugLog_.push_back(logString);
        size_t d = calibrator_->model_->factors() - 1;
        if (d < 2) return objective;  // we don't include regularisation for one-factor models
        if ((calibrator_->penaltySigma_<=0.0)&&(calibrator_->penaltySlope_<=0.0)&&(calibrator_->penaltyCurve_<=0.0))
            return objective;  // we only include regularisation for non-trivial parameters
        // add regularisation here...
        idx = 0;
        // count relevant entries
        for (size_t i = 0; i<isInput_.size(); ++i) {
            for (size_t j = 0; j<d; ++j) {
                if (isInput_[i][j]>0.0) {
                    ++idx;
                    break; // one sigma suffices
                }
            }
            for (size_t j = 0; j<d; ++j) {
                if (isInput_[i][d + j]>0.0) {
                    ++idx;
                    break; // one slope suffices
                }
            }
            for (size_t j = 0; j<d; ++j) {
                if (isInput_[i][2*d + j]>0.0) {
                    ++idx;
                    break; // one curve suffices
                }
            }
        }
        // calculate actual penalty term
        Array reg((d - 1)*idx);
        idx = 0;
        for (size_t i = 0; i<isInput_.size(); ++i) {
            for (size_t j = 0; j<d; ++j) {
                if (isInput_[i][j]>0.0) {
                    for (size_t k=1; k<d; ++k) {
                        reg[idx] = (calibrator_->model_->sigma()[k][i] - calibrator_->model_->sigma()[k - 1][i]) * calibrator_->penaltySigma_;
                        ++idx;
                    }
                    break; // one sigma suffices
                }
            }
            for (size_t j = 0; j<d; ++j) {
                if (isInput_[i][d + j]>0.0) {
                    for (size_t k=1; k<d; ++k) {
                        reg[idx] = (calibrator_->model_->slope()[k][i] - calibrator_->model_->slope()[k - 1][i]) * calibrator_->penaltySlope_;
                        ++idx;
                    }
                    break; // one slope suffices
                }
            }
            for (size_t j = 0; j<d; ++j) {
                if (isInput_[i][2*d + j]>0.0) {
                    for (size_t k=1; k<d; ++k) {
                        reg[idx] = (calibrator_->model_->curve()[k][i] - calibrator_->model_->curve()[k - 1][i]) * calibrator_->penaltyCurve_;
                        ++idx;
                    }
                    break; // one slope suffices
                }
            }
        }
        // Logging
        // logString = "reg = ";
        // for (size_t k = 0; k < reg.size(); ++k) logString = logString + " " + std::to_string(reg[k]) + ",";
        // calibrator_->debugLog_.push_back(logString);
        // concatenate data
        Array objectiveAndRegularisation(objective.size() + reg.size());
        for (size_t k = 0; k < objective.size(); ++k) objectiveAndRegularisation[k] = objective[k];
        for (size_t k = 0; k < reg.size(); ++k) objectiveAndRegularisation[objective.size() + k] = reg[k];
        return objectiveAndRegularisation;
    }

    Real QGMonteCarloCalibrator::Objective::value(const Array& x) const {
        Real sum = 0.0;
        for (size_t i=0; i<x.size(); ++i) sum += x[i]*x[i];
        return 0.5*sum;
    }

    const ext::shared_ptr<QGMonteCarloCalibrator::QuasiGaussianModel> QGMonteCarloCalibrator::Objective::model(const Array& x) {
        update(x);
        return calibrator_->model_;
    }

    // constructor
    QGMonteCarloCalibrator::QGMonteCarloCalibrator(
        const ext::shared_ptr<QGMonteCarloCalibrator::QuasiGaussianModel>&  model,
        const Handle<SwaptionVolatilityStructure>&            volTS,
        const std::vector< ext::shared_ptr<SwapIndex> >&    swapIndices,
        const Real                                            monteCarloStepSize,
        const Size                                            monteCarloPaths,
        const Real                                            sigmaMax,
        const Real                                            slopeMax,
        const Real                                            curveMax,
        const Real                                            sigmaWeight,
        const Real                                            slopeWeight,
        const Real                                            curveWeight,
        const Real                                            penaltySigma,
        const Real                                            penaltySlope,
        const Real                                            penaltyCurve,
        const EndCriteria&                                    endCriteria )
        : model_(model->clone()), volTS_(volTS), swapIndices_(swapIndices),  // we clone the model to be able to work on it directly
          sigmaMin_(0.0), sigmaMax_(sigmaMax), slopeMin_(0.0), slopeMax_(slopeMax), curveMin_(0.0), curveMax_(curveMax),
          sigmaWeight_(sigmaWeight), slopeWeight_(slopeWeight), curveWeight_(curveWeight),
          penaltySigma_(penaltySigma), penaltySlope_(penaltySlope), penaltyCurve_(penaltyCurve),
          endCriteria_(endCriteria) {
        // maybe check Feller condition: eta^2 < 2 theta
        // we need to set up the Monte-Carlo simulation...
        // the critical step here is to determine the simulation times appropriately
        QL_REQUIRE(monteCarloStepSize > 0.0, "monteCarloStepSize > 0.0 requires");
        QL_REQUIRE(monteCarloPaths > 0, "monteCarloPaths > 0 requires");
        std::vector<Real> simulationTimes(1, 0.0);
        Size idx = 0;
        while (simulationTimes.back() < model->times().back()) {
            Real Tnext = simulationTimes.back() + monteCarloStepSize;  // that is the default
            if (Tnext > model->times()[idx] - 1.0 / 365) {  // we want to avoid MC steps smaller than a single day
                Tnext = model->times()[idx];
                ++idx;
            }
            simulationTimes.push_back(Tnext);
        }
        mcSimulation_ = ext::shared_ptr<MCSimulation>(new MCSimulation(model_,
            simulationTimes, simulationTimes, monteCarloPaths, 1234, false, true, true));  // we need to allow time interpolation
        mcSimulation_->prepareForSlicedSimulation();  // we need this for sliced simulation later on
        mcSimulation_->simulate(0);  // we want to initialise the initial states

        // now we can do the bootstrapping
        for (size_t k = 0; k < model_->times().size(); ++k) {  // we bootstrap parameters
            std::vector< std::vector< Real > > isInput(model_->times().size(), std::vector< Real >(3*(model_->factors()-1),0.0));
            std::vector< std::vector< Real > > isOutput(model_->times().size(), std::vector< Real >(3 * swapIndices_.size(), 0.0));
            if (sigmaMax_ > 0.0) {
                for (size_t j = 0; j < model_->factors() - 1; ++j) {
                    isInput[k][j] = (sigmaWeight_>0.0) ? (1.0) : (0.0);  // we do not scale inputs
                }
                for (size_t j = 0; j < swapIndices_.size(); ++j) {
                    isOutput[k][j] = (sigmaWeight_>0.0) ? (sigmaWeight_) : (0.0);
                }
            }
            if (slopeMax_ > 0.0) {
                for (size_t j = 0; j < model_->factors() - 1; ++j) {
                    isInput[k][j + (model_->factors() - 1)] = (slopeWeight_>0.0) ? (1.0) : (0.0);  // we do not scale inputs
                }
                for (size_t j = 0; j < swapIndices_.size(); ++j) {
                    isOutput[k][j + swapIndices_.size()] = (slopeWeight_>0.0) ? (slopeWeight_) : (0.0); 
                }
            }
            if (curveMax_ > 0.0) {
                for (size_t j = 0; j < model_->factors() - 1; ++j) {
                    isInput[k][j + 2*(model_->factors() - 1)] = (curveWeight_>0.0) ? (1.0) : (0.0);  // we do not scale inputs
                }
                for (size_t j = 0; j < swapIndices_.size(); ++j) {
                    isOutput[k][j + 2*swapIndices_.size()] = (curveWeight_>0.0) ? (curveWeight_) : (0.0);
                }
            }
            Real epsfcn = 1.0e-4;  // this is not very sophisticated
            Integer info = calibrate(isInput, isOutput, epsfcn);		
            debugLog_.push_back(std::string("k = " + std::to_string(k) + ", info = " + std::to_string(info)));
            acceptCalibration();  // maybe better check info
        }
    }

    Integer QGMonteCarloCalibrator::calibrate(
                                   const std::vector< std::vector< Real > >&  isInput,
                                   const std::vector< std::vector< Real > >&  isOutput,
                                    // optimization parameters
                                   Real                                       epsfcn ) {  // delta for finite differences
        // set up
        NoConstraint constraint;
        Objective obj(this, isInput, isOutput);
        Array x = obj.initialise();
        Problem problem(obj, constraint, x);
        LevenbergMarquardt optimizationMethod(epsfcn, endCriteria_.rootEpsilon(), endCriteria_.gradientNormEpsilon());  // (epsfcn, xtol, gtol)
        // EndCriteria endCriteria(maxfev, 100 /* unused */, 0 /* unused */, ftol, 0 /* unused */);
        // calibrate
        optimizationMethod.minimize(problem,endCriteria_);  // here we use maxfev and ftol
        calibratedModel_ = obj.model(problem.currentValue());
        return optimizationMethod.getInfo();
    }


}

