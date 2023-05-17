/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2018 Sebastian Schlenkrich

*/


#include <ql/experimental/templatemodels/qgaussian2/qgcalibrator.hpp>

#include <ql/math/optimization/constraint.hpp>
#include <ql/math/optimization/levenbergmarquardt.hpp>
#include <ql/pricingengines/blackformula.hpp>


namespace QuantLib {


    QGCalibrator::Objective::CalibSwaption::CalibSwaption (
        Date                                expiryDate,
        const ext::shared_ptr<SwapIndex>&   swapindex,
        const Handle<YieldTermStructure>&   discountCurve,
        const Handle<SwaptionVolatilityStructure> volTS,
        bool                                contTenorSpread,
        Real                                modelTimesStepSize )
        : SwapCashFlows(swapindex->underlyingSwap(expiryDate),discountCurve,contTenorSpread) {
        // first we need to set up the model times for the swap rate model
        Time Tend = Actual365Fixed().yearFraction(discountCurve->referenceDate(), expiryDate);
        Size N = (Size) (Tend/modelTimesStepSize);  // this is rounded towards the lower
        modelTimes_.resize(N+1); // we need to store zero as well
        modelTimes_[0] = 0.0;
        for (Size i=1; i<=N; ++i) {
            modelTimes_[i] = modelTimes_[i-1] + modelTimesStepSize;
        }
        if (modelTimes_[N] < Tend - 1.0 / 365.0) modelTimes_.push_back(Tend);
        else modelTimes_[N] = Tend;
        // then we may calculate volatilities
        Rate S0 = swapindex->fixing(expiryDate, true);
        sigmaATM_ = volTS->volatility(expiryDate, swapindex->tenor(), S0, true);
        Real sigmap1Std = volTS->volatility(expiryDate, swapindex->tenor(), S0 + sigmaATM_*sqrt(Tend), true);
        Real sigmam1Std = volTS->volatility(expiryDate, swapindex->tenor(), S0 - sigmaATM_*sqrt(Tend), true);
        skew_ = sigmap1Std - sigmam1Std;
        smile_ = sigmap1Std + sigmam1Std - 2.0 * sigmaATM_;
    }


    QGCalibrator::Objective::Objective(
                  QGCalibrator                               *calibrator,
                  const std::vector< std::vector< Real > >&  isInput,
                  const std::vector< std::vector< Real > >&  isOutput ) 
        : calibrator_(calibrator), isInput_(isInput), isOutput_(isOutput) {
        // copy model initial values
        model_ = calibrator->model_->clone();
        // checking dimensions
        QL_REQUIRE(isInput_.size()==model_->times().size(),"QGCalibrator::Objective: wrong input dimension.");
        for (size_t i=0; i<isInput_.size(); ++i) {
            QL_REQUIRE(isInput_[i].size()==(2*model_->factors()-1),"QGCalibrator::Objective: wrong input dimension.");
        }
        QL_REQUIRE(isOutput_.size()== model_->times().size(),"QGCalibrator::Objective: wrong output dimension.");
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
                        new CalibSwaption(expiryDate, calibrator_->swapIndices_[j], calibrator_->model_->termStructure(), calibrator_->volTS_, true, calibrator_->modelTimesStepSize_));
                }
            }
        }
        // ready for optimisation...
    }

    Array QGCalibrator::Objective::initialise() {
        Array X(inputSize_);
        size_t d = model_->factors()-1;
        size_t idx = 0;
        for (size_t i=0; i<isInput_.size(); ++i) {
            for (size_t j=0; j<d; ++j) {
                if (isInput_[i][j]>0.0) {
                    Real sigma = model_->sigma()[j][i];
                    Real x = QGCalibrator::inverse(sigma,calibrator_->sigmaMin_,calibrator_->sigmaMax_) / isInput_[i][j];
                    X[idx] = x;
                    ++idx;
                }
            }
            for (size_t j=0; j<d; ++j) {
                if (isInput_[i][d+j]>0.0) {
                    Real b = model_->slope()[j][i];
                    Real x = QGCalibrator::inverse(b,calibrator_->slopeMin_,calibrator_->slopeMax_) / isInput_[i][d+j];
                    X[idx] = x;
                    ++idx;
                }
            }
            if (isInput_[i][d+d]>0.0) {
                Real eta = model_->eta()[i];
                Real x = QGCalibrator::inverse(eta,calibrator_->etaMin_,calibrator_->etaMax_) / isInput_[i][d+d];
                X[idx] = x;
                ++idx;
            }
        }
        return X;
    }

    void QGCalibrator::Objective::update(const Array& X) const {
        std::vector< std::vector< Real > > m_sigma = model_->sigma();
        std::vector< std::vector< Real > > m_slope = model_->slope();
        std::vector< Real >                m_eta   = model_->eta();
        size_t d = model_->factors()-1;
        size_t idx = 0;
        for (size_t i=0; i<isInput_.size(); ++i) {
            for (size_t j=0; j<d; ++j) {
                if (isInput_[i][j]>0.0) {
                    Real x = X[idx];
                    Real sigma = QGCalibrator::direct(x*isInput_[i][j],calibrator_->sigmaMin_,calibrator_->sigmaMax_);
                    m_sigma[j][i] = sigma;
                    ++idx;
                }
            }
            for (size_t j=0; j<d; ++j) {
                if (isInput_[i][d+j]>0.0) {
                    Real x = X[idx];
                    Real slope = QGCalibrator::direct(x*isInput_[i][d+j],calibrator_->slopeMin_,calibrator_->slopeMax_);
                    m_slope[j][i] = slope;
                    ++idx;
                }
            }
            if (isInput_[i][d+d]>0.0) {
                Real x = X[idx];
                Real eta = QGCalibrator::direct(x*isInput_[i][d+d],calibrator_->etaMin_,calibrator_->etaMax_);
                m_eta[i] = eta;
                ++idx;
            }
        }
        model_->update(m_sigma,m_slope,m_eta);
    }

    Array QGCalibrator::Objective::values(const Array& x) const {
        update(x);
        // we may have a swaption model for each swaption
        std::vector< std::vector< ext::shared_ptr<QGAverageSwaprateModel> > > swaprateModels;
        // however we build the model only if neccessary
        swaprateModels.resize(calibSwaptions_.size());
        for (size_t i=0; i<swaprateModels.size(); ++i) {
            swaprateModels[i].resize(calibSwaptions_[i].size());
            for (size_t j=0; j<swaprateModels[i].size(); ++j) {
                if (calibSwaptions_[i][j]) {
                    ext::shared_ptr<QGSwaprateModel> tmp(new QGSwaprateModel(model_, calibSwaptions_[i][j]->floatTimes(), calibSwaptions_[i][j]->floatWeights(), calibSwaptions_[i][j]->fixedTimes(), calibSwaptions_[i][j]->annuityWeights(), calibSwaptions_[i][j]->modelTimes(), calibrator_->useExpectedXY_));
                    swaprateModels[i][j] = ext::shared_ptr<QGAverageSwaprateModel>( new QGAverageSwaprateModel( tmp ) );
                }
            }
        }
        // averaging and objective calculation
        Array objective(outputSize_);
        size_t idx=0;
        for (size_t i=0; i<swaprateModels.size(); ++i) {
            for (size_t j=0; j<swaprateModels[i].size(); ++j) {
                Real sigmaATM=0.0, skew=0.0, smile=0.0;
                if ((isOutput_[i][j] > 0.0) ||
                    (isOutput_[i][j + calibSwaptions_[i].size()] > 0.0) ||
                    (isOutput_[i][j + 2 * calibSwaptions_[i].size()] > 0.0)) {
                    Real S0 = swaprateModels[i][j]->S0();
                    Real callATM  = swaprateModels[i][j]->vanillaOption(S0, 1, 1.0e-6, 1000);
                    Real expTime = swaprateModels[i][j]->modelTimes()[swaprateModels[i][j]->modelTimes().size() - 1];
                    sigmaATM = bachelierBlackFormulaImpliedVol(Option::Call, S0, S0, expTime, callATM);
                    Real putMinus = swaprateModels[i][j]->vanillaOption(S0 - calibSwaptions_[i][j]->sigmaATM()*sqrt(expTime), -1, 1.0e-6, 1000);
                    Real callPlus = swaprateModels[i][j]->vanillaOption(S0 + calibSwaptions_[i][j]->sigmaATM()*sqrt(expTime),  1, 1.0e-6, 1000);
                    Real sigmaMns = bachelierBlackFormulaImpliedVol(Option::Put, S0 - calibSwaptions_[i][j]->sigmaATM()*sqrt(expTime), S0, expTime, putMinus);
                    Real sigmaPls = bachelierBlackFormulaImpliedVol(Option::Call, S0 + calibSwaptions_[i][j]->sigmaATM()*sqrt(expTime), S0, expTime, callPlus);
                    skew = sigmaPls - sigmaMns;
                    smile = sigmaPls + sigmaMns - 2.0*sigmaATM;
                }
                if (isOutput_[i][j]>0.0) {
                    objective[idx] = (sigmaATM - calibSwaptions_[i][j]->sigmaATM()) * isOutput_[i][j];
                    ++idx;
                }
                if (isOutput_[i][j+calibSwaptions_[i].size()]>0.0) {
                    objective[idx] = (skew - calibSwaptions_[i][j]->skew()) * isOutput_[i][j + calibSwaptions_[i].size()];
                    ++idx;
                }
                if (isOutput_[i][j+2*calibSwaptions_[i].size()]>0.0) {
                    objective[idx] = (smile - calibSwaptions_[i][j]->smile()) * isOutput_[i][j + 2 * calibSwaptions_[i].size()];
                    ++idx;
                }
            }
        }
        // Logging
        // std::string logString("objective = ");
        // for (size_t k = 0; k < objective.size(); ++k) logString = logString + " " + std::to_string(objective[k]) + ",";
        // calibrator_->debugLog_.push_back(logString);
        size_t d = model_->factors() - 1;
        if (d < 2) return objective;  // we don't include regularisation for one-factor models
        if ((calibrator_->penaltySigma_<=0.0)&&(calibrator_->penaltySlope_<=0.0))
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
        }
        // calculate actual penalty term
        Array reg((d - 1)*idx);
        idx = 0;
        for (size_t i = 0; i<isInput_.size(); ++i) {
            for (size_t j = 0; j<d; ++j) {
                if (isInput_[i][j]>0.0) {
                    for (size_t k=1; k<d; ++k) {
                        reg[idx] = (model_->sigma()[k][i] - model_->sigma()[k - 1][i]) * calibrator_->penaltySigma_;
                        ++idx;
                    }
                    break; // one sigma suffices
                }
            }
            for (size_t j = 0; j<d; ++j) {
                if (isInput_[i][d + j]>0.0) {
                    for (size_t k = 1; k<d; ++k) {
                        reg[idx] = (model_->slope()[k][i] - model_->slope()[k - 1][i]) * calibrator_->penaltySlope_;
                        ++idx;
                    }
                    break; // one slo�e suffices
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

    Real QGCalibrator::Objective::value(const Array& x) const {
        Real sum = 0.0;
        for (size_t i=0; i<x.size(); ++i) sum += x[i]*x[i];
        return 0.5*sum;
    }

    const ext::shared_ptr<QuasiGaussianModel> QGCalibrator::Objective::model(const Array& x) {
        update(x);
        return model_;
    }

    // constructor
    QGCalibrator::QGCalibrator(
        const ext::shared_ptr<QuasiGaussianModel>&          model,
        const Handle<SwaptionVolatilityStructure>&            volTS,
        const std::vector< ext::shared_ptr<SwapIndex> >&    swapIndices,
        const Real                                            modelTimesStepSize,
        const bool                                            useExpectedXY,
        const Real                                            sigmaMax,
        const Real                                            slopeMax,
        const Real                                            etaMax,
        const Real                                            sigmaWeight,
        const Real                                            slopeWeight,
        const Real                                            etaWeight,
        const Real                                            penaltySigma,
        const Real                                            penaltySlope,
        const EndCriteria&                                    endCriteria )
        : model_(model), volTS_(volTS), swapIndices_(swapIndices), modelTimesStepSize_(modelTimesStepSize), 
          useExpectedXY_(useExpectedXY),
          sigmaMin_(0.0), sigmaMax_(sigmaMax), slopeMin_(0.0), slopeMax_(slopeMax), etaMin_(0.0), etaMax_(etaMax),
          sigmaWeight_(sigmaWeight), slopeWeight_(slopeWeight), etaWeight_(etaWeight),
          penaltySigma_(penaltySigma), penaltySlope_(penaltySlope), endCriteria_(endCriteria) {
        // maybe check Feller condition: eta^2 < 2 theta
        for (size_t k = 0; k < model_->times().size(); ++k) {  // we bootstrap parameters
            std::vector< std::vector< Real > > isInput(model_->times().size(), std::vector< Real >(2*model_->factors()-1,0.0));
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
            if (etaMax_ > 0.0) {
                isInput[k][2 * (model_->factors() - 1)] = (etaWeight_>0.0) ? (1.0) : (0.0); // we do not scale inputs
                for (size_t j = 0; j < swapIndices_.size(); ++j) {  // maybe we should better only calibrate against the fisrt index...?
                    isOutput[k][j + 2* swapIndices_.size()] = (etaWeight_>0.0) ? (etaWeight_) : (0.0);
                }
            }
            Real epsfcn = 1.0e-4;  // this is not very sophisticated
            Integer info = calibrate(isInput, isOutput, epsfcn);		
            debugLog_.push_back(std::string("k = " + std::to_string(k) + ", info = " + std::to_string(info)));
            acceptCalibration();  // maybe better check info
        }
    }

    Integer QGCalibrator::calibrate(
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

