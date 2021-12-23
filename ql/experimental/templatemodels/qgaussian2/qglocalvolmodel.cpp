/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2018 Sebastian Schlenkrich

*/

#include <ql/pricingengines/blackformula.hpp>
#include <ql/math/interpolations/linearinterpolation.hpp>
#include <ql/termstructures/volatility/smilesection.hpp>

#include <ql/experimental/templatemodels/qgaussian2/qglocalvolmodel.hpp>

namespace QuantLib {

    // model implementaion
    QGLocalvolModel::QGLocalvolModel(
        const Handle<YieldTermStructure>&                      termStructure,
        const Handle<SwaptionVolatilityStructure>&             volTS,
        const Real                                             chi,
        const ext::shared_ptr<SwapIndex>&                      swapIndex,
        const std::vector< Real >&                             times,
        const std::vector<Real>&                               stdDevGrid,
        const size_t                                           nPaths,
        const BigNatural                                       seed,
        const size_t                                           debugLevel)
        : QuasiGaussianModel(termStructure, 1, times,
            std::vector< std::vector<Real> >(1, std::vector<Real>(times.size(), 0.0)),  // sigma
            std::vector< std::vector<Real> >(1, std::vector<Real>(times.size(), 0.0)),  // slope
            std::vector< std::vector<Real> >(1, std::vector<Real>(times.size(), 0.0)),  // curve
            std::vector<Real>(times.size(), 0.0),                                       // eta
            std::vector<Real>(1, 1.0),                                                  // delta
            std::vector<Real>(1, chi),                                                  // mean reversion
            std::vector< std::vector<Real> >(1, std::vector<Real>(1, 1.0)),             // Gamma
            0.1                                                                         // theta
        ),
        volTS_(volTS), swapIndex_(swapIndex), stdDevGrid_(stdDevGrid), nPaths_(nPaths), seed_(seed),
        calcStochVolAdjustment_(false), kernelWidth_(0.0), debugLevel_(debugLevel) {
        // we can't calibrate here because *this needs to be assigned to shared pointer first
    }

    QGLocalvolModel::QGLocalvolModel(
        const Handle<YieldTermStructure>&                      termStructure,
        const Handle<SwaptionVolatilityStructure>&             volTS,
        const Real                                             chi,
        const Real                                             theta,
        const Real                                             eta,
        const ext::shared_ptr<SwapIndex>&                      swapIndex,
        const std::vector< Real >&                             times,
        const std::vector<Real>&                               stdDevGrid,
        const bool                                             calcStochVolAdjustment,
        const Real                                             kernelWidth,
        const size_t                                           nPaths,
        const BigNatural                                       seed,
        const size_t                                           debugLevel)
        : QuasiGaussianModel(termStructure, 1, times,
            std::vector< std::vector<Real> >(1, std::vector<Real>(times.size(), 0.0)),  // sigma
            std::vector< std::vector<Real> >(1, std::vector<Real>(times.size(), 0.0)),  // slope
            std::vector< std::vector<Real> >(1, std::vector<Real>(times.size(), 0.0)),  // curve
            std::vector<Real>(times.size(), eta),                                       // eta
            std::vector<Real>(1, 1.0),                                                  // delta
            std::vector<Real>(1, chi),                                                  // mean reversion
            std::vector< std::vector<Real> >(1, std::vector<Real>(1, 1.0)),             // Gamma
            theta                                                                       // theta
        ),
        volTS_(volTS), swapIndex_(swapIndex), stdDevGrid_(stdDevGrid), nPaths_(nPaths), seed_(seed),
        calcStochVolAdjustment_(calcStochVolAdjustment), kernelWidth_(kernelWidth), debugLevel_(debugLevel) {
        // we can't calibrate here because *this needs to be assigned to shared pointer first
    }

    inline std::vector< std::vector<Real> >
    QGLocalvolModel::sigma_xT(const Real t, const State& s) {
        if (sigmaMode_ == Parent) return QuasiGaussianModel::sigma_xT(t, s);  // this is more like a fall back if volatility is irrelevant
        if (sigmaMode_ == Calibration) {
            // the current implementation is intended only for calibration phase
            // we rely on the swap rate model beeing initialised and updated properly before this function is called
            // this will most likely not work if model is used with external MC simulation
            QL_REQUIRE(swapRateModel_ != 0, "swapRateModel_!=0 required");
            Real observationTime = swapRateModel_->modelTimes().back();  // this should be times()[idx] from calibrateAndSimulate()
            std::vector<Real> x(1, s.x(0));
            std::vector< std::vector<Real> > y(1, std::vector<Real>(1, s.y(0, 0)));
            Real swapRate = swapRateModel_->swapRate(observationTime, x, y);
            Real swapGradient = (swapRateModel_->swapGradient(observationTime, x, y))[0];
            Real lvol;
            if (sigmaS_.size() == 0) {
                lvol = volTS_->volatility(observationTime, swapIndex_->tenor(), swapRate); // initial approximation
            }
            else {  // ensure flat extrapolation..
                bool warning = false;
                if (swapRate < sigmaS_.back().xMin()) {
                    if (debugLevel_ > 2) debugLog_.push_back("Warning: xMin = " + std::to_string(sigmaS_.back().xMin()) + ", swapRate = " + std::to_string(swapRate));
                    swapRate = sigmaS_.back().xMin();
                    warning = true;
                }
                if (swapRate > sigmaS_.back().xMax()) {
                    if (debugLevel_ > 2) debugLog_.push_back("Warning: xMax = " + std::to_string(sigmaS_.back().xMax()) + ", swapRate = " + std::to_string(swapRate));
                    swapRate = sigmaS_.back().xMax();
                    warning = true;
                }
                lvol = sigmaS_.back()(swapRate);  // this relies on the usage within calibrateAndSimulate()													   
                if ((warning)&&(debugLevel_ > 2)) {
                    debugLog_.push_back("Warning: lvol = " + std::to_string(lvol) + ", swapGradient = " + std::to_string(swapGradient));
                }
            }
            if (debugLevel_ > 3) debugLog_.push_back("t = " + std::to_string(t) +
                                                     ", x = " + std::to_string(x[0]) +
                                                     ", y = " + std::to_string(y[0][0]) +
                                                     ", swapRate = " + std::to_string(swapRate) +
                                                     ", swapGradient = " + std::to_string(swapGradient) +
                                                     ", lvol = " + std::to_string(lvol) );
            return  std::vector< std::vector<Real> >(1, std::vector<Real>(1, lvol / swapGradient));
        }
        QL_REQUIRE(false, "post-calibration mode for sigma_x not implemented yet");
    }

    // we cache the float leg and fixed leg valuation to speed-up MC simulation
    QGLocalvolModel::SwaptionFactory::SwaptionFactory(const Time obsTime, const SwapCashFlows& scf)
        : floatLeg_(ext::shared_ptr<QGLocalvolModel::MCPayoff>(new QGLocalvolModel::MCBase::Cache(
            ext::shared_ptr<QGLocalvolModel::MCPayoff>(new MCAnnuity(obsTime, scf.floatTimes(), scf.floatWeights()))))),
        annuityLeg_(ext::shared_ptr<QGLocalvolModel::MCPayoff>(new QGLocalvolModel::MCBase::Cache(
            ext::shared_ptr<QGLocalvolModel::MCPayoff>(new MCAnnuity(obsTime, scf.fixedTimes(), scf.annuityWeights())))))
    {}

    ext::shared_ptr<QGLocalvolModel::MCPayoff> QGLocalvolModel::SwaptionFactory::swaption(const Real strike, const Real callOrPut) {
        ext::shared_ptr<QGLocalvolModel::MCPayoff> fixedRate(new QGLocalvolModel::MCBase::FixedAmount(strike));
        ext::shared_ptr<QGLocalvolModel::MCPayoff> fixedLeg(new QGLocalvolModel::MCBase::Mult(fixedRate, annuityLeg_));
        ext::shared_ptr<QGLocalvolModel::MCPayoff> swap;
        if (callOrPut == 1.0) {
            swap = ext::shared_ptr<QGLocalvolModel::MCPayoff>(new QGLocalvolModel::MCBase::Axpy(-1.0, fixedLeg, floatLeg_));
        }
        else {
            swap = ext::shared_ptr<QGLocalvolModel::MCPayoff>(new QGLocalvolModel::MCBase::Axpy(-1.0, floatLeg_, fixedLeg));
        }
        ext::shared_ptr<QGLocalvolModel::MCPayoff> zero(new QGLocalvolModel::MCBase::FixedAmount(0.0));
        ext::shared_ptr<QGLocalvolModel::MCPayoff> swpt(new QGLocalvolModel::MCBase::Max(swap, zero));
        ext::shared_ptr<QGLocalvolModel::MCPayoff> pay(new QGLocalvolModel::MCBase::Pay(swpt, floatLeg_->observationTime()));
        return pay;
    }

    std::vector< std::vector<Real> > QGLocalvolModel::calibrationTest(const std::vector<Date>& exerciseDates,
                                                                      const std::vector<Real>&  stdDevStrikes) {
        // add adjusters...
        std::vector< std::vector<Real> > resultTable;
        // we need to do some sanity checks 
        QL_REQUIRE(simulation_, "non-empty simulation_ required");
        // derive adjuster dates
        std::vector<Real> adjObsTimes;
        for (size_t i = 0; i < exerciseDates.size(); ++i) {
            if (exerciseDates[i] <= termStructure()->referenceDate()) continue; // skip dates in the past
            Date fixingDate = exerciseDates[i];
            if (!swapIndex_->isValidFixingDate(fixingDate)) fixingDate = swapIndex_->fixingCalendar().adjust(fixingDate, Following);
            Time exerciseTime = Actual365Fixed().yearFraction(termStructure()->referenceDate(), fixingDate);
            adjObsTimes.push_back(exerciseTime);
        }
        // derive offset dates
        std::vector<Real> adjOffsetTimes;
        if (swapIndex_->tenor().units() == Months) {
            adjOffsetTimes.push_back(swapIndex_->tenor().length() / 12.0);
        }
        if (swapIndex_->tenor().units() == Years) {
            for (Integer i=1; i<=swapIndex_->tenor().length(); ++i) adjOffsetTimes.push_back((Real)i);
        }
        QL_REQUIRE(adjOffsetTimes.size() > 0, "adjOffsetTimes.size()>0 required");
        // calculate adjusters
        simulation_->calculateNumeraireAdjuster(adjObsTimes);
        simulation_->calculateZCBAdjuster(adjObsTimes, adjOffsetTimes);
        // now we may calculate instruments...
        for (size_t i = 0; i < exerciseDates.size(); ++i) {
            if (exerciseDates[i] <= termStructure()->referenceDate()) continue; // skip dates in the past
            Date fixingDate = exerciseDates[i];
            if (!swapIndex_->isValidFixingDate(fixingDate)) fixingDate = swapIndex_->fixingCalendar().adjust(fixingDate, Following);
            Time exerciseTime = Actual365Fixed().yearFraction(termStructure()->referenceDate(), fixingDate);
            SwapCashFlows scf(swapIndex_->underlyingSwap(fixingDate), termStructure_, true);
            SwaptionFactory testFactory(exerciseTime, scf);

            // swap rate etc.
            Real annuity = 0.0;
            for (size_t k = 0; k < scf.fixedTimes().size(); ++k) annuity += scf.annuityWeights()[k] * termStructure()->discount(scf.fixedTimes()[k]);
            Real floatLeg = 0.0;
            for (size_t k = 0; k < scf.floatTimes().size(); ++k) floatLeg += scf.floatWeights()[k] * termStructure()->discount(scf.floatTimes()[k]);
            Real swapRate = floatLeg / annuity;

            // set up smile section and strike grid
            ext::shared_ptr<SmileSection> smileSection = volTS_->smileSection(exerciseTime, swapIndex_->tenor(), true);
            Real stdDev = smileSection->optionPrice(swapRate, Option::Call) / M_1_SQRTPI / M_SQRT_2;
            std::vector<Real> strikeGrid(stdDevStrikes);
            for (size_t k = 0; k < strikeGrid.size(); ++k) strikeGrid[k] = strikeGrid[k] * stdDev + swapRate;

            std::vector< ext::shared_ptr<MCPayoff> > options;
            for (size_t k = 0; k < strikeGrid.size(); ++k) {
                options.push_back(testFactory.swaption(strikeGrid[k], 1.0));   // call
                options.push_back(testFactory.swaption(strikeGrid[k], -1.0));  // put
            }
            // cached version of MC simulation
            std::vector<Real> optionNPVs = MCPayoff::Pricer::NPVs(options, simulation_);

            for (size_t k = 0; k < strikeGrid.size(); ++k) {
                Real call = smileSection->optionPrice(strikeGrid[k], Option::Call);
                Real testCall = optionNPVs[2*k];
                Real testPut = optionNPVs[2*k+1];
                testCall /= annuity;
                testPut /= annuity;
                Real callVol = 0.0;
                Real putVol = 0.0;
                Real vanVol = 0.0;
                try {
                    callVol = bachelierBlackFormulaImpliedVol(Option::Call, strikeGrid[k], swapRate, exerciseTime, testCall);
                } catch (std::exception const& e) {}
                try {
                    putVol = bachelierBlackFormulaImpliedVol(Option::Put, strikeGrid[k], swapRate, exerciseTime, testPut);
                } catch (std::exception const& e) {}
                try {
                    vanVol = bachelierBlackFormulaImpliedVol(Option::Call, strikeGrid[k], swapRate, exerciseTime, call);
                } catch (std::exception const& e) {}
                std::vector<Real> resultRow = { (double)exerciseDates[i].serialNumber(),
                                                (double)fixingDate.serialNumber(),
                                                exerciseTime,
                                                swapRate,
                                                (double)k,
                                                stdDevStrikes[k],
                                                strikeGrid[k],
                                                call,
                                                testCall,
                                                testPut,
                                                vanVol,
                                                callVol,
                                                putVol                                       };
                resultTable.push_back(resultRow);
            }
        }
        return resultTable;
    }


    inline size_t QGLocalvolModel::minIdx(const std::vector<Real>& X, const Real x) {
        QL_REQUIRE(X.size() > 0, "X.size>0 required");
        if (x <= X[0]) return 0;
        if (x > X.back()) return X.size();
        // bisection search
        size_t a = 0, b = X.size() - 1;
        while (b - a > 1) {
            size_t s = (a + b) / 2;
            if (x <= X[s]) b = s;
            else           a = s;
        }
        return b;
    }

    ext::shared_ptr<QGLocalvolModel::QGSwaprateModel> QGLocalvolModel::qGSwapRateModel(const SwapCashFlows& scf, const Real obsTime) {
        std::vector<Real> swapRateModelTimes{ 0.0, obsTime };
        sigmaMode_ = Parent; // switch-off sigma_x calculation; we only need yield curve information
        ext::shared_ptr<QGSwaprateModel> swapRateModel(new QGSwaprateModel(ext::static_pointer_cast<QuasiGaussianModel>(shared_from_this()),
            scf.floatTimes(), scf.floatWeights(), scf.fixedTimes(), scf.annuityWeights(), swapRateModelTimes, false));
        sigmaMode_ = Calibration; // switch-on sigma_x calculation for simulation and evolve calls
        return swapRateModel;
    }


    QGLocalvolModel::Initialiser::Initialiser(ext::shared_ptr<QGLocalvolModel> model) {
        // reset local volatility attributes and debugging
        model->debugLog_.clear();
        model->sigmaS_.clear();
        model->strikeGrid_.clear();
        model->locvolGrid_.clear();

        // initialise MC simulation
        model->simulation_ = ext::shared_ptr<MCSimulation>(new MCSimulation(model, model->times(), model->times(), model->nPaths_, model->seed_, false, true, true));
        model->simulation_->prepareForSlicedSimulation();
        model->simulation_->simulate(0);
        QL_REQUIRE(model->simulation_->simTimes().size() == model->times().size() + 1, "simulation_->simTimes().size()==times().size()+1 required.");

        // prepare for first simulation step using approximate local vol
        today_ = model->termStructure_->referenceDate(); // check if this is the correct date...
                                                  // we need to set up an initial swap rate model for gradient calculation
        Date fixingDate = model->swapIndex_->fixingCalendar().advance(today_, 1, Days, Following);
        SwapCashFlows scf(model->swapIndex_->underlyingSwap(fixingDate), model->termStructure_, true);
        model->swapRateModel_ = model->qGSwapRateModel(scf, Actual365Fixed().yearFraction(today_, fixingDate));
    }

    // setup swap rate etc. for given time point
    QGLocalvolModel::SwapRate::SwapRate(const QGLocalvolModel *model, const Date today, const Real fixingTime) {
        fixingDate_ = today + ((BigInteger)ClosestRounding(0)(fixingTime * 365.0)); // assuming act/365 day counting
        if (!model->swapIndex_->isValidFixingDate(fixingDate_)) fixingDate_ = model->swapIndex_->fixingCalendar().adjust(fixingDate_, Following);
        scf_ = SwapCashFlows(model->swapIndex_->underlyingSwap(fixingDate_), model->termStructure_, true);        // assume continuous tenor spreads
        annuity_ = 0.0;
        for (size_t k = 0; k < scf_.fixedTimes().size(); ++k) annuity_ += scf_.annuityWeights()[k] * model->termStructure_->discount(scf_.fixedTimes()[k]);
        Real floatLeg = 0.0;
        for (size_t k = 0; k < scf_.floatTimes().size(); ++k) floatLeg += scf_.floatWeights()[k] * model->termStructure_->discount(scf_.floatTimes()[k]);
        swapRate_ = floatLeg / annuity_;
    }

    QGLocalvolModel::McCalculator::McCalculator(
        QGLocalvolModel          *model,  // model cannot be const coz we push into debugLog_
        const Real               obsTime,
        const SwapCashFlows&     scf,
        const Real               annuity,
        const Real               swapRate,
        const std::vector<Real>& smileStrikeGrid) 
        : oneOverBSample_(model->simulation_->nPaths(),0.0),
        annuitySample_(model->simulation_->nPaths(), 0.0),
        swapRateSample_(model->simulation_->nPaths(), 0.0),
        vanillaOptions_(smileStrikeGrid.size(), 0.0),
        avgCalcStrikes_(0.0) {
        ext::shared_ptr<QGLocalvolModel::MCPayoff> mcFloatLeg(new MCAnnuity(obsTime, scf.floatTimes(), scf.floatWeights()));
        ext::shared_ptr<QGLocalvolModel::MCPayoff> mcFixedLeg(new MCAnnuity(obsTime, scf.fixedTimes(), scf.annuityWeights()));
        ext::shared_ptr<QGLocalvolModel::MCPayoff> one(new QGLocalvolModel::MCBase::FixedAmount(1.0));
        ext::shared_ptr<QGLocalvolModel::MCPayoff> oneAtT(new QGLocalvolModel::MCBase::Pay(one, obsTime));
        for (size_t k = 0; k < model->simulation_->nPaths(); ++k) {
            ext::shared_ptr<MCSimulation::Path> p = model->simulation_->path(k);
            oneOverBSample_[k] = oneAtT->discountedAt(p);
            annuitySample_[k] = mcFixedLeg->at(p);
            swapRateSample_[k] = mcFloatLeg->at(p) / annuitySample_[k];
        }
        // calculate adjusters suksessively
        Real mcDF = 0.0;
        for (size_t k = 0; k < model->simulation_->nPaths(); ++k) mcDF += oneOverBSample_[k];
        mcDF /= model->simulation_->nPaths();
        Real adjOneOverB = model->termStructure_->discount(obsTime) / mcDF;
        for (size_t k = 0; k < model->simulation_->nPaths(); ++k) oneOverBSample_[k] *= adjOneOverB;
        Real mcAnnuity = 0.0;
        for (size_t k = 0; k < model->simulation_->nPaths(); ++k) mcAnnuity += (annuitySample_[k] * oneOverBSample_[k]);
        mcAnnuity /= model->simulation_->nPaths();
        Real adjAnnuity = annuity / mcAnnuity;
        for (size_t k = 0; k < model->simulation_->nPaths(); ++k) annuitySample_[k] *= adjAnnuity;
        Real mcFloat = 0.0;
        for (size_t k = 0; k < model->simulation_->nPaths(); ++k) mcFloat += (annuitySample_[k] * swapRateSample_[k] * oneOverBSample_[k]);
        mcFloat /= model->simulation_->nPaths();
        Real adjSwapRate = swapRate - mcFloat / annuity;
        for (size_t k = 0; k < model->simulation_->nPaths(); ++k) swapRateSample_[k] += adjSwapRate;
        // find index s.t. strike[idx]>=swapRate
        size_t callIdx = smileStrikeGrid.size();
        for (size_t j = smileStrikeGrid.size(); j > 0; --j)
            if (smileStrikeGrid[j - 1] >= swapRate) callIdx = j - 1; // not very efficient
        
        // calculate out-of-the-money option prices
        for (size_t k = 0; k < model->simulation_->nPaths(); ++k) {
            // find index in ascending vector, evaluate n s.t. strike[n-1] < s <= strike[n]
            size_t strikeIdx = minIdx(smileStrikeGrid, swapRateSample_[k]);
            avgCalcStrikes_ += fabs(1.0*strikeIdx - 1.0*callIdx);
            // calculate put option prices
            // for (size_t j = callIdx; j > strikeIdx; --j) {
            for (size_t j = strikeIdx; j < callIdx; ++j) {
                Real value = smileStrikeGrid[j] - swapRateSample_[k];
                if (value < 0.0) {
                    value = 0.0; // better safe
                    if (model->debugLevel_ > 2) model->debugLog_.push_back("Warning (put): callIdx = " + std::to_string(callIdx) + ", k = " + std::to_string(k) + ", strikeIdx = " + std::to_string(strikeIdx) + ", j = " + std::to_string(j));
                }
                vanillaOptions_[j] += (annuitySample_[k] * value*oneOverBSample_[k]);
            }
            // calculate call option prices
            for (size_t j = callIdx; j < strikeIdx; ++j) {
                Real value = swapRateSample_[k] - smileStrikeGrid[j];
                if (value < 0.0) {
                    value = 0.0; // better safe
                    if (model->debugLevel_ > 2) model->debugLog_.push_back("Warning (call): callIdx = " + std::to_string(callIdx) + ", k = " + std::to_string(k) + ", strikeIdx = " + std::to_string(strikeIdx) + ", j = " + std::to_string(j));
                }
                vanillaOptions_[j] += (annuitySample_[k] * value*oneOverBSample_[k]);
            }
        }
        avgCalcStrikes_ = avgCalcStrikes_ / model->simulation_->nPaths();
        for (size_t j = 0; j < vanillaOptions_.size(); ++j) vanillaOptions_[j] = vanillaOptions_[j] / model->simulation_->nPaths() / annuity;
        // translate put into call prices
        for (size_t j = 0; j < vanillaOptions_.size(); ++j) {
            if (smileStrikeGrid[j] < swapRate) {
                vanillaOptions_[j] = vanillaOptions_[j] + swapRate - smileStrikeGrid[j];
            }
        }

        if (model->debugLevel_ > 2) {
            // brute-force double-check call and put calculation
            std::vector<Real> vanillaCalls(smileStrikeGrid.size(), 0.0);
            std::vector<Real> vanillaPuts(smileStrikeGrid.size(), 0.0);
            for (size_t k = 0; k < model->simulation_->nPaths(); ++k) {
                for (size_t j = 0; j < smileStrikeGrid.size(); ++j) {
                    Real call = swapRateSample_[k] - smileStrikeGrid[j];
                    Real put = -call;
                    if (call < 0.0) call = 0.0;
                    if (put < 0.0) put = 0.0;
                    vanillaCalls[j] += (annuitySample_[k] * call * oneOverBSample_[k]);
                    vanillaPuts[j] += (annuitySample_[k] * put * oneOverBSample_[k]);
                }
            }
            for (size_t j = 0; j < smileStrikeGrid.size(); ++j) {
                vanillaCalls[j] = vanillaCalls[j] / model->simulation_->nPaths() / annuity;
                vanillaPuts[j] = vanillaPuts[j] / model->simulation_->nPaths() / annuity;
            }
            for (size_t j = 0; j < smileStrikeGrid.size(); ++j) {
                model->debugLog_.push_back("obsTime = " + std::to_string(obsTime) + ", swapRate = " + std::to_string(swapRate) + ", j = " + std::to_string(j) + ", strike = " + std::to_string(smileStrikeGrid[j]) + ", vanilla = " + std::to_string(vanillaOptions_[j]) + ", call = " + std::to_string(vanillaCalls[j]) + ", put = " + std::to_string(vanillaPuts[j]) + ", P/C(bp) = " + std::to_string((vanillaCalls[j] - vanillaPuts[j] - (swapRate - smileStrikeGrid[j]))*1.0e4) + ", (V-C)bp = " + std::to_string((vanillaOptions_[j] - vanillaCalls[j])*1.0e4));
            }
        }
    }

    QGLocalvolModel::StochvolExpectation::StochvolExpectation(
            const QGLocalvolModel  *model,
            const size_t           simIdx,
            const Real             lambda,   // lambda = 1.0 / kernelWidth / stdDev
            const Real             annuity,
            const McCalculator&    mcCalc,
            const std::vector<Real>&  strikeGrid,
            Real                  (*kernel)(const Real)
        )
            : expectationZCondS_(strikeGrid.size(), 0.0) {
        // first we need to extract z(T) from the simulation
        std::vector<Real> stochVarianceSample(model->simulation_->nPaths());
        for (size_t k = 0; k < model->simulation_->nPaths(); ++k) stochVarianceSample[k] = model->simulation_->observedPath(k)[simIdx][2];
        // now we can iterate the paths and calculate conditional expectations
        // E^A[.] = E^Q[ z(T)*q(T) | S(T)=K ] / E^Q[ q(T) | S(T)=K ]
        // E^Q[.] = sum{ z_i*q_i*Kernel(S_i) } / n
        // R.-N.-Derivative q(T) = N(0)/An(0) * An(T)/N(T)
        std::vector<Real> zTimesQGrid(strikeGrid.size(), 0.0);  // numerator
        std::vector<Real> qGrid(strikeGrid.size(), 0.0);  // denumerator
        for (size_t k = 0; k < model->simulation_->nPaths(); ++k) {
            size_t startIdx = minIdx(strikeGrid, mcCalc.swapRateSample()[k] - 1.0/lambda);
            size_t endIdx = minIdx(strikeGrid, mcCalc.swapRateSample()[k] + 1.0 / lambda);
            Real q = mcCalc.annuitySample()[k] * mcCalc.oneOverBSample()[k] / annuity;  // N(0)=1
            for (size_t j = startIdx; j < endIdx; ++j) {
                Real kernelAtStrike = lambda * (*kernel)(lambda*(mcCalc.swapRateSample()[k] - strikeGrid[j]));
                zTimesQGrid[j] += stochVarianceSample[k] * q * kernelAtStrike;
                qGrid[j] += q * kernelAtStrike;
            }
        }
        // finally adjust sigma_SV = sigma_LV / E^A[ z(T) | S(T) = K ]
        for (size_t j = 0; j < strikeGrid.size(); ++j) expectationZCondS_[j] = zTimesQGrid[j] / qGrid[j];
    }



    inline void QGLocalvolModel::checkMCPrices( const Real               obsTime,
                                                const SwapCashFlows&     scf,
                                                const Real               annuity,
                                                const Real               swapRate,
                                                const std::vector<Real>& smileStrikeGrid ) {
        SwaptionFactory testFactory(obsTime, scf);
        for (size_t k = 0; k < smileStrikeGrid.size(); ++k) {
            ext::shared_ptr<MCPayoff> mcCall = testFactory.swaption(smileStrikeGrid[k], 1.0);
            ext::shared_ptr<MCPayoff> mcPut = testFactory.swaption(smileStrikeGrid[k], -1.0);
            Real testCall = MCPayoff::Pricer::NPV(std::vector< ext::shared_ptr<MCPayoff> >(1, mcCall), simulation_);
            Real testPut = MCPayoff::Pricer::NPV(std::vector< ext::shared_ptr<MCPayoff> >(1, mcPut), simulation_);
            testCall /= annuity;
            testPut /= annuity;
            debugLog_.push_back("T = " + std::to_string(obsTime) + ", swapRate = " + std::to_string(swapRate) + ", k = " + std::to_string(k) + ", strike = " + std::to_string(smileStrikeGrid[k]) + ", testCall = " + std::to_string(testCall) + ", testPut = " + std::to_string(testPut) );
            try {
                Real callVol = bachelierBlackFormulaImpliedVol(Option::Call, smileStrikeGrid[k], swapRate, obsTime, testCall);
                Real putVol = bachelierBlackFormulaImpliedVol(Option::Put, smileStrikeGrid[k], swapRate, obsTime, testPut);
                Real vanVol = volTS_->volatility(obsTime, swapIndex_->tenor(), smileStrikeGrid[k], true);
                debugLog_.push_back("obsTime = " + std::to_string(obsTime) + ", swapRate = " + std::to_string(swapRate) + ", k = " + std::to_string(k) + ", strike = " + std::to_string(smileStrikeGrid[k]) + ", callVol = " + std::to_string(callVol) + ", putVol = " + std::to_string(putVol) + ", vanVol = " + std::to_string(vanVol));
            }
            catch (std::exception const& e) {
                std::string what = e.what();
                debugLog_.push_back("Error: " + what);
            }
        }
    }


    void QGLocalvolModelBackwardFlavor::simulateAndCalibrate() {
        Initialiser init(ext::static_pointer_cast<QGLocalvolModelBackwardFlavor>(shared_from_this()));
        // run simulation and calibration
        for (size_t idx = 0; idx < times().size(); ++idx) {
            // simulate next step T[idx-1] to T[idx], note X[] starts at 0
            simulation_->simulate(idx + 1);

            // specify swap rate at T[idx]
            SwapRate swapRate(this, init.today(), times()[idx]);

            // set up smile section and strike grid
            ext::shared_ptr<SmileSection> smileSection = volTS_->smileSection(times()[idx], swapIndex_->tenor(), true);
            Real stdDev = smileSection->optionPrice(swapRate.swapRate(), Option::Call) / M_1_SQRTPI / M_SQRT_2;
            std::vector<Real> initialStrikeGrid(stdDevGrid_);
            for (size_t k = 0; k < initialStrikeGrid.size(); ++k) initialStrikeGrid[k] = initialStrikeGrid[k] * stdDev + swapRate.swapRate();

            // these are the grids used for smile section calculation
            std::vector<Real> smileStrikeGrid;
            std::vector<Real> callGrid;
            std::vector<Real> d2CalldK2Grid;
            for (size_t k = 0; k < initialStrikeGrid.size(); ++k) {

                Real h = 1.0e-5; // 0.1bp shift size
                Real call = smileSection->optionPrice(initialStrikeGrid[k], Option::Call);
                Real callph = smileSection->optionPrice(initialStrikeGrid[k] + h, Option::Call);
                Real callmh = smileSection->optionPrice(initialStrikeGrid[k] - h, Option::Call);

                // we only calculate local vol if within [1e-6, 1 - 1e-6] quantile 
                Real quantile = 1.0 + (callph - callmh) / 2.0 / h;
                Real d2CalldK2 = (callph + callmh - 2.0*call) / h / h;
                if (debugLevel_>1) debugLog_.push_back("T = " + std::to_string(times()[idx]) + ", swapRate = " + std::to_string(swapRate.swapRate()) + ", k = " + std::to_string(k) + ", strike = " + std::to_string(initialStrikeGrid[k]) + ", call = " + std::to_string(call) + ", callph = " + std::to_string(callph) + ", callmh = " + std::to_string(callmh) + ", quantile = " + std::to_string(quantile) + ", d2CalldK2 = " + std::to_string(d2CalldK2));
                if ((quantile < 1.0e-6) || (quantile > 1.0 - 1.e-6)) {
                    if (debugLevel_>1) debugLog_.push_back("Warning: Skip local vol.");
                    continue;  // skip calculation
                }
                if (d2CalldK2 < 1.0e-8) {
                    if (debugLevel_>1) debugLog_.push_back("Warning: d2CalldK2 = " + std::to_string(d2CalldK2) + "Skip local vol.");
                    continue;  // skip calculation
                }
                smileStrikeGrid.push_back(initialStrikeGrid[k]);
                callGrid.push_back(call);
                d2CalldK2Grid.push_back(d2CalldK2);
            }
            // test calibration
            if (debugLevel_ > 2) checkMCPrices(times()[idx], swapRate.scf(), swapRate.annuity(), swapRate.swapRate(), initialStrikeGrid); // we might want to evaluate MC swaptions for debugging

            // calculate MC option prices
            Time obsTime = (idx>0) ? (times()[idx - 1]) : (0.0);
            McCalculator mcCalc(this, obsTime, swapRate.scf(), swapRate.annuity(), swapRate.swapRate(), smileStrikeGrid);

            std::vector<Real> strikeGrid;
            std::vector<Real> dCalldTGrid;
            std::vector<Real> localVol;
            // calculate dC/dT and setup final strikes and vols
            for (size_t j = 0; j < mcCalc.vanillaOptions().size(); ++j) {
                Real dCalldT = (callGrid[j] - mcCalc.vanillaOptions()[j]) / (times()[idx] - obsTime);
                if (dCalldT < 1.0e-8) {
                    if (debugLevel_>1) debugLog_.push_back("Warning: dCalldT = " + std::to_string(dCalldT) + ", Skip local vol.");
                    continue;  // skip calculation
                }
                Real lvol = sqrt(2 * dCalldT / d2CalldK2Grid[j]);
                strikeGrid.push_back(smileStrikeGrid[j]);
                dCalldTGrid.push_back(dCalldT);
                localVol.push_back(lvol);
            }

            // logging
            if (debugLevel_ > 0) debugLog_.push_back("T = " + std::to_string(times()[idx]) + ", fixingDate = " + std::to_string(swapRate.fixingDate().serialNumber()) + ", swapRate = " + std::to_string(swapRate.swapRate()) + ", annuity = " + std::to_string(swapRate.annuity()) + ", stdDev = " + std::to_string(stdDev) + ", avgCalcStrikes = " + std::to_string(mcCalc.avgCalcStrikes()));
            if (debugLevel_ > 1) {
                for (size_t j = 0; j < strikeGrid.size(); ++j) {
                    debugLog_.push_back("strike = " + std::to_string(strikeGrid[j]) + ", dCalldT = " + std::to_string(dCalldTGrid[j]) + ", localVol = " + std::to_string(localVol[j]));
                }
            }

            // calculate E^A[ z(T) | S(T) = K ] and adjust local vol
            if (calcStochVolAdjustment_) {
                StochvolExpectation expZ(this, idx, 1.0 / kernelWidth_ / stdDev, swapRate.annuity(), mcCalc, strikeGrid, &kernel);				
                for (size_t j = 0; j < strikeGrid.size(); ++j) localVol[j] = localVol[j] / sqrt(expZ.expectationZCondS()[j]);  // finally adjust sigma_SV = sigma_LV / sqrt(E^A[ z(T) | S(T) = K ])
            }

            // set up interpolation
            // sigmaS_.push_back(LinearInterpolation(strikeGrid.begin(), strikeGrid.end(), localVol.begin()));
            strikeGrid_.push_back(strikeGrid);
            locvolGrid_.push_back(localVol);
            sigmaS_.push_back(Linear().interpolate(strikeGrid_.back().begin(), strikeGrid_.back().end(), locvolGrid_.back().begin()));

            // set up swap rate model here such that we don't need to care about it during sigma_x calculation
            swapRateModel_ = qGSwapRateModel(swapRate.scf(), times()[idx]);
        }
    }

    void QGLocalvolModelForwardFlavor::simulateAndCalibrate() {
        QL_REQUIRE(volTS_->volatilityType() == VolatilityType::Normal, "Normal volatilities required.");
        Initialiser init(ext::static_pointer_cast<QGLocalvolModelBackwardFlavor>(shared_from_this()));
        simulation_->simulate(1);

        // run actual simulation and calibration
        for (size_t idx = 1; idx < times().size(); ++idx) {
            // specify swap rate fixing at T[idx], but observed at T[idx-1]
            SwapRate swapRate(this, init.today(), times()[idx]);

            // set up smile section and strike grid
            ext::shared_ptr<SmileSection> smileSection = volTS_->smileSection(times()[idx], swapIndex_->tenor(), true);  // the vanilla model only provides terminal distribution of swap rate at T[idx]; at T[idx-1] it is a different swap rate
            Real stdDev = smileSection->volatility(swapRate.swapRate()) * sqrt(times()[idx - 1]);  // approximated ATM forward vol
            std::vector<Real> initialStrikeGrid(stdDevGrid_);
            for (size_t k = 0; k < initialStrikeGrid.size(); ++k) initialStrikeGrid[k] = initialStrikeGrid[k] * stdDev + swapRate.swapRate();

            // these are the grids used for smile section calculation
            std::vector<Real> smileStrikeGrid;
            std::vector<Real> callGrid;
            std::vector<Real> d2CalldK2Grid;
            for (size_t k = 0; k < initialStrikeGrid.size(); ++k) {

                Real h = 1.0e-5; // 0.1bp shift size
                // basic idea is to approximate forward smile based on terminal smile from vanilla model
                Real call = bachelierBlackFormula(Option::Call, initialStrikeGrid[k], swapRate.swapRate(), smileSection->volatility(initialStrikeGrid[k])*sqrt(times()[idx - 1]));
                Real callph = bachelierBlackFormula(Option::Call, initialStrikeGrid[k]+h, swapRate.swapRate(), smileSection->volatility(initialStrikeGrid[k]+h)*sqrt(times()[idx - 1]));
                Real callmh = bachelierBlackFormula(Option::Call, initialStrikeGrid[k]-h, swapRate.swapRate(), smileSection->volatility(initialStrikeGrid[k]-h)*sqrt(times()[idx - 1]));

                // we only calculate local vol if within [1e-6, 1 - 1e-6] quantile 
                Real quantile = 1.0 + (callph - callmh) / 2.0 / h;
                Real d2CalldK2 = (callph + callmh - 2.0*call) / h / h;
                if (debugLevel_>1) debugLog_.push_back("T = " + std::to_string(times()[idx]) + ", swapRate = " + std::to_string(swapRate.swapRate()) + ", k = " + std::to_string(k) + ", strike = " + std::to_string(initialStrikeGrid[k]) + ", call = " + std::to_string(call) + ", callph = " + std::to_string(callph) + ", callmh = " + std::to_string(callmh) + ", quantile = " + std::to_string(quantile) + ", d2CalldK2 = " + std::to_string(d2CalldK2));
                if ((quantile < 1.0e-6) || (quantile > 1.0 - 1.e-6)) {
                    if (debugLevel_>1) debugLog_.push_back("Warning: Skip local vol.");
                    continue;  // skip calculation
                }
                if (d2CalldK2 < 1.0e-8) { // this check is not trivial since our forward vol/price approximation is not guarantied arb-free even if terminal smile is arb-free
                    if (debugLevel_>1) debugLog_.push_back("Warning: d2CalldK2 = " + std::to_string(d2CalldK2) + "Skip local vol.");
                    continue;  // skip calculation
                }
                // we need to re-calculate call price at the next grid point because this is needed for dCall /dT
                call = smileSection->optionPrice(initialStrikeGrid[k], Option::Call);

                smileStrikeGrid.push_back(initialStrikeGrid[k]);
                callGrid.push_back(call);
                d2CalldK2Grid.push_back(d2CalldK2);
            }

            // calculate MC option prices
            Time obsTime = (idx>0) ? (times()[idx - 1]) : (0.0);  // this should always be T[idx-1] here
            McCalculator mcCalc(this, obsTime, swapRate.scf(), swapRate.annuity(), swapRate.swapRate(), smileStrikeGrid);

            std::vector<Real> strikeGrid;
            std::vector<Real> dCalldTGrid;
            std::vector<Real> localVol;
            // calculate dC/dT and setup final strikes and vols
            for (size_t j = 0; j < mcCalc.vanillaOptions().size(); ++j) {
                Real dCalldT = (callGrid[j] - mcCalc.vanillaOptions()[j]) / (times()[idx] - obsTime);
                if (dCalldT < 1.0e-8) {
                    if (debugLevel_>1) debugLog_.push_back("Warning: dCalldT = " + std::to_string(dCalldT) + ", Skip local vol.");
                    continue;  // skip calculation
                }
                Real lvol = sqrt(2 * dCalldT / d2CalldK2Grid[j]);
                strikeGrid.push_back(smileStrikeGrid[j]);
                dCalldTGrid.push_back(dCalldT);
                localVol.push_back(lvol);
            }

            // logging
            if (debugLevel_ > 0) debugLog_.push_back("T = " + std::to_string(times()[idx]) + ", fixingDate = " + std::to_string(swapRate.fixingDate().serialNumber()) + ", swapRate = " + std::to_string(swapRate.swapRate()) + ", annuity = " + std::to_string(swapRate.annuity()) + ", stdDev = " + std::to_string(stdDev) + ", avgCalcStrikes = " + std::to_string(mcCalc.avgCalcStrikes()));
            if (debugLevel_ > 1) {
                for (size_t j = 0; j < strikeGrid.size(); ++j) {
                    debugLog_.push_back("strike = " + std::to_string(strikeGrid[j]) + ", dCalldT = " + std::to_string(dCalldTGrid[j]) + ", localVol = " + std::to_string(localVol[j]));
                }
            }

            // calculate E^A[ z(T) | S(T) = K ] and adjust local vol
            if (calcStochVolAdjustment_) {
                StochvolExpectation expZ(this, idx, 1.0 / kernelWidth_ / stdDev, swapRate.annuity(), mcCalc, strikeGrid, &kernel);
                for (size_t j = 0; j < strikeGrid.size(); ++j) localVol[j] = localVol[j] / sqrt(expZ.expectationZCondS()[j]);  // finally adjust sigma_SV = sigma_LV / sqrt(E^A[ z(T) | S(T) = K ])
            }

            // set up interpolation
            // sigmaS_.push_back(LinearInterpolation(strikeGrid.begin(), strikeGrid.end(), localVol.begin()));
            strikeGrid_.push_back(strikeGrid);
            locvolGrid_.push_back(localVol);
            sigmaS_.push_back(Linear().interpolate(strikeGrid_.back().begin(), strikeGrid_.back().end(), locvolGrid_.back().begin()));

            // set up swap rate model here such that we don't need to care about it during sigma_x calculation
            swapRateModel_ = qGSwapRateModel(swapRate.scf(), times()[idx-1]);

            // simulate the next step T[idx-1] to T[idx]
            simulation_->simulate(idx + 1);

            // test calibration
            if (debugLevel_ > 2) checkMCPrices(times()[idx], swapRate.scf(), swapRate.annuity(), swapRate.swapRate(), smileStrikeGrid); // we might want to evaluate MC swaptions for debugging
        }
    }

    void QGLocalvolModelAnalyticFlavor::simulateAndCalibrate() {
        QL_REQUIRE(volTS_->volatilityType() == VolatilityType::Normal, "Normal volatilities required.");
        Initialiser init(ext::static_pointer_cast<QGLocalvolModelBackwardFlavor>(shared_from_this()));
        simulation_->simulate(1);

        // run actual simulation and calibration
        for (size_t idx = 1; idx < times().size(); ++idx) {
            // specify swap rate fixing at T[idx], but observed at T[idx-1]
            SwapRate swapRate(this, init.today(), times()[idx]);

            // set up smile section and strike grid
            ext::shared_ptr<SmileSection> smileSection = volTS_->smileSection(times()[idx], swapIndex_->tenor(), true);  // the vanilla model only provides terminal distribution of swap rate at T[idx]; at T[idx-1] it is a different swap rate
            Real stdDev = smileSection->volatility(swapRate.swapRate()) * sqrt(times()[idx - 1]);  // approximated ATM forward vol
            std::vector<Real> initialStrikeGrid(stdDevGrid_);
            for (size_t k = 0; k < initialStrikeGrid.size(); ++k) initialStrikeGrid[k] = initialStrikeGrid[k] * stdDev + swapRate.swapRate();

            // these are the grids used for smile section calculation
            std::vector<Real> smileStrikeGrid;
            std::vector<Real> callGrid0, callGrid1;  // T[idx-1] and T[idx]
            std::vector<Real> d2CalldK2Grid;
            for (size_t k = 0; k < initialStrikeGrid.size(); ++k) {

                Real h = 1.0e-5; // 0.1bp shift size
                                 // basic idea is to approximate forward smile based on terminal smile from vanilla model
                Real call0  = bachelierBlackFormula(Option::Call, initialStrikeGrid[k], swapRate.swapRate(), smileSection->volatility(initialStrikeGrid[k])*sqrt(times()[idx - 1]));
                Real callph = bachelierBlackFormula(Option::Call, initialStrikeGrid[k] + h, swapRate.swapRate(), smileSection->volatility(initialStrikeGrid[k] + h)*sqrt(times()[idx - 1]));
                Real callmh = bachelierBlackFormula(Option::Call, initialStrikeGrid[k] - h, swapRate.swapRate(), smileSection->volatility(initialStrikeGrid[k] - h)*sqrt(times()[idx - 1]));

                // we only calculate local vol if within [1e-6, 1 - 1e-6] quantile 
                Real quantile = 1.0 + (callph - callmh) / 2.0 / h;
                Real d2CalldK2 = (callph + callmh - 2.0*call0) / h / h;
                if (debugLevel_>1) debugLog_.push_back("T = " + std::to_string(times()[idx]) + ", swapRate = " + std::to_string(swapRate.swapRate()) + ", k = " + std::to_string(k) + ", strike = " + std::to_string(initialStrikeGrid[k]) + ", call0 = " + std::to_string(call0) + ", callph = " + std::to_string(callph) + ", callmh = " + std::to_string(callmh) + ", quantile = " + std::to_string(quantile) + ", d2CalldK2 = " + std::to_string(d2CalldK2));
                if ((quantile < 1.0e-6) || (quantile > 1.0 - 1.e-6)) {
                    if (debugLevel_>1) debugLog_.push_back("Warning: Skip local vol.");
                    continue;  // skip calculation
                }
                if (d2CalldK2 < 1.0e-8) { // this check is not trivial since our forward vol/price approximation is not guarantied arb-free even if terminal smile is arb-free
                    if (debugLevel_>1) debugLog_.push_back("Warning: d2CalldK2 = " + std::to_string(d2CalldK2) + "Skip local vol.");
                    continue;  // skip calculation
                }
                // we need to re-calculate call price at the next grid point because this is needed for dCall /dT
                Real call1 = smileSection->optionPrice(initialStrikeGrid[k], Option::Call);

                smileStrikeGrid.push_back(initialStrikeGrid[k]);
                callGrid0.push_back(call0);
                callGrid1.push_back(call1);
                d2CalldK2Grid.push_back(d2CalldK2);
            }

            std::vector<Real> strikeGrid;
            std::vector<Real> dCalldTGrid;
            std::vector<Real> localVol;
            // calculate dC/dT and setup final strikes and vols
            for (size_t j = 0; j < smileStrikeGrid.size(); ++j) {
                Real dCalldT = (callGrid1[j] - callGrid0[j]) / (times()[idx] - times()[idx-1]);
                if (dCalldT < 1.0e-8) {
                    if (debugLevel_>1) debugLog_.push_back("Warning: dCalldT = " + std::to_string(dCalldT) + ", Skip local vol.");
                    continue;  // skip calculation
                }
                Real lvol = sqrt(2 * dCalldT / d2CalldK2Grid[j]);
                strikeGrid.push_back(smileStrikeGrid[j]);
                dCalldTGrid.push_back(dCalldT);
                localVol.push_back(lvol);
            }

            // logging
            if (debugLevel_ > 0) debugLog_.push_back("T = " + std::to_string(times()[idx]) + ", fixingDate = " + std::to_string(swapRate.fixingDate().serialNumber()) + ", swapRate = " + std::to_string(swapRate.swapRate()) + ", annuity = " + std::to_string(swapRate.annuity()) + ", stdDev = " + std::to_string(stdDev) );
            if (debugLevel_ > 1) {
                for (size_t j = 0; j < strikeGrid.size(); ++j) {
                    debugLog_.push_back("strike = " + std::to_string(strikeGrid[j]) + ", dCalldT = " + std::to_string(dCalldTGrid[j]) + ", localVol = " + std::to_string(localVol[j]));
                }
            }

            // set up interpolation
            // sigmaS_.push_back(LinearInterpolation(strikeGrid.begin(), strikeGrid.end(), localVol.begin()));
            strikeGrid_.push_back(strikeGrid);
            locvolGrid_.push_back(localVol);
            sigmaS_.push_back(Linear().interpolate(strikeGrid_.back().begin(), strikeGrid_.back().end(), locvolGrid_.back().begin()));

            // set up swap rate model here such that we don't need to care about it during sigma_x calculation
            swapRateModel_ = qGSwapRateModel(swapRate.scf(), times()[idx - 1]);

            // simulate the next step T[idx-1] to T[idx]
            simulation_->simulate(idx + 1);

            // test calibration
            if (debugLevel_ > 2) checkMCPrices(times()[idx], swapRate.scf(), swapRate.annuity(), swapRate.swapRate(), smileStrikeGrid); // we might want to evaluate MC swaptions for debugging
        }
    }


    void QGLocalvolModelMonteCarloFlavor::simulateAndCalibrate() {
        QL_REQUIRE(volTS_->volatilityType() == VolatilityType::Normal, "Normal volatilities required.");
        Initialiser init(ext::static_pointer_cast<QGLocalvolModelBackwardFlavor>(shared_from_this()));
        simulation_->simulate(1);

        // run actual simulation and calibration
        for (size_t idx = 1; idx < times().size(); ++idx) {
            // specify swap rate fixing at T[idx], but observed at T[idx-1]
            SwapRate swapRate(this, init.today(), times()[idx]);

            // set up smile section and strike grid
            ext::shared_ptr<SmileSection> smileSection = volTS_->smileSection(times()[idx], swapIndex_->tenor(), true);  // the vanilla model only provides terminal distribution of swap rate at T[idx]; at T[idx-1] it is a different swap rate
            Real stdDev = smileSection->volatility(swapRate.swapRate()) * sqrt(times()[idx - 1]);  // approximated ATM forward vol
            std::vector<Real> initialStrikeGrid(stdDevGrid_);
            for (size_t k = 0; k < initialStrikeGrid.size(); ++k) initialStrikeGrid[k] = initialStrikeGrid[k] * stdDev + swapRate.swapRate();

            // these are the grids used for smile section calculation
            std::vector<Real> smileStrikeGrid;
            std::vector<Real> callGrid;
            std::vector<Real> d2CalldK2Grid;
            for (size_t k = 0; k < initialStrikeGrid.size(); ++k) {

                Real h = 1.0e-5; // 0.1bp shift size
                                 // basic idea is to approximate forward smile based on terminal smile from vanilla model
                Real call = bachelierBlackFormula(Option::Call, initialStrikeGrid[k], swapRate.swapRate(), smileSection->volatility(initialStrikeGrid[k])*sqrt(times()[idx - 1]));
                Real callph = bachelierBlackFormula(Option::Call, initialStrikeGrid[k] + h, swapRate.swapRate(), smileSection->volatility(initialStrikeGrid[k] + h)*sqrt(times()[idx - 1]));
                Real callmh = bachelierBlackFormula(Option::Call, initialStrikeGrid[k] - h, swapRate.swapRate(), smileSection->volatility(initialStrikeGrid[k] - h)*sqrt(times()[idx - 1]));

                // we only calculate local vol if within [1e-6, 1 - 1e-6] quantile 
                Real quantile = 1.0 + (callph - callmh) / 2.0 / h;
                Real d2CalldK2 = (callph + callmh - 2.0*call) / h / h;
                if (debugLevel_>1) debugLog_.push_back("T = " + std::to_string(times()[idx]) + ", swapRate = " + std::to_string(swapRate.swapRate()) + ", k = " + std::to_string(k) + ", strike = " + std::to_string(initialStrikeGrid[k]) + ", call = " + std::to_string(call) + ", callph = " + std::to_string(callph) + ", callmh = " + std::to_string(callmh) + ", quantile = " + std::to_string(quantile) + ", d2CalldK2 = " + std::to_string(d2CalldK2));
                if ((quantile < 1.0e-6) || (quantile > 1.0 - 1.e-6)) {
                    if (debugLevel_>1) debugLog_.push_back("Warning: Skip local vol.");
                    continue;  // skip calculation
                }
                if (d2CalldK2 < 1.0e-8) { // this check is not trivial since our forward vol/price approximation is not guarantied arb-free even if terminal smile is arb-free
                    if (debugLevel_>1) debugLog_.push_back("Warning: d2CalldK2 = " + std::to_string(d2CalldK2) + "Skip local vol.");
                    continue;  // skip calculation
                }
                // we need to re-calculate call price at the next grid point because this is needed for dCall /dT
                call = smileSection->optionPrice(initialStrikeGrid[k], Option::Call);

                smileStrikeGrid.push_back(initialStrikeGrid[k]);
                callGrid.push_back(call);
                d2CalldK2Grid.push_back(d2CalldK2);
            }

            // calculate MC option prices
            Time obsTime = (idx>0) ? (times()[idx - 1]) : (0.0);  // this should always be T[idx-1] here
            McCalculator mcCalc(this, obsTime, swapRate.scf(), swapRate.annuity(), swapRate.swapRate(), smileStrikeGrid);

            std::vector<Real> strikeGrid;
            std::vector<Real> dCalldTGrid;
            std::vector<Real> localVol;
            // calculate dC/dT and d^2C/dK^2 based on MC prices  and setup final strikes and vols
            for (size_t j = 1; j < mcCalc.vanillaOptions().size()-1; ++j) { // we can only calculate d^2C/dK^2 for inner grid points
                Real dCalldT = (callGrid[j] - mcCalc.vanillaOptions()[j]) / (times()[idx] - obsTime);
                if (dCalldT < 1.0e-8) {
                    if (debugLevel_>1) debugLog_.push_back("Warning: dCalldT = " + std::to_string(dCalldT) + ", Skip local vol.");
                    continue;  // skip calculation
                }
                // MC-based calculation
                Real dCalldKph = (mcCalc.vanillaOptions()[j+1] - mcCalc.vanillaOptions()[j]) / (smileStrikeGrid[j+1] - smileStrikeGrid[j]);
                Real dCalldKmh = (mcCalc.vanillaOptions()[j] - mcCalc.vanillaOptions()[j-1]) / (smileStrikeGrid[j] - smileStrikeGrid[j-1]);
                Real d2CalldK2mc = (dCalldKph - dCalldKmh) / (smileStrikeGrid[j+1] - smileStrikeGrid[j-1]) * 2.0;
                if (d2CalldK2mc < 1.0e-8) { // catch numerical issues with mc prices 
                    if (debugLevel_>1) debugLog_.push_back("Warning: d2CalldK2mc = " + std::to_string(d2CalldK2mc) + "Skip local vol.");
                    continue;  // skip calculation
                }
                if (debugLevel_ > 1) {
                    debugLog_.push_back("strike = " + std::to_string(smileStrikeGrid[j]) + ", d2CalldK2 = " + std::to_string(d2CalldK2Grid[j]) + ", d2CalldK2mc = " + std::to_string(d2CalldK2mc));
                }

                Real lvol = sqrt(2 * dCalldT / d2CalldK2mc);
                strikeGrid.push_back(smileStrikeGrid[j]);
                dCalldTGrid.push_back(dCalldT);
                localVol.push_back(lvol);
            }

            // logging
            if (debugLevel_ > 0) debugLog_.push_back("T = " + std::to_string(times()[idx]) + ", fixingDate = " + std::to_string(swapRate.fixingDate().serialNumber()) + ", swapRate = " + std::to_string(swapRate.swapRate()) + ", annuity = " + std::to_string(swapRate.annuity()) + ", stdDev = " + std::to_string(stdDev) + ", avgCalcStrikes = " + std::to_string(mcCalc.avgCalcStrikes()));
            if (debugLevel_ > 1) {
                for (size_t j = 0; j < strikeGrid.size(); ++j) {
                    debugLog_.push_back("strike = " + std::to_string(strikeGrid[j]) + ", dCalldT = " + std::to_string(dCalldTGrid[j]) + ", localVol = " + std::to_string(localVol[j]));
                }
            }

            // calculate E^A[ z(T) | S(T) = K ] and adjust local vol
            if (calcStochVolAdjustment_) {
                StochvolExpectation expZ(this, idx, 1.0 / kernelWidth_ / stdDev, swapRate.annuity(), mcCalc, strikeGrid, &kernel);
                for (size_t j = 0; j < strikeGrid.size(); ++j) localVol[j] = localVol[j] / sqrt(expZ.expectationZCondS()[j]);  // finally adjust sigma_SV = sigma_LV / sqrt(E^A[ z(T) | S(T) = K ])
            }

            // set up interpolation
            // sigmaS_.push_back(LinearInterpolation(strikeGrid.begin(), strikeGrid.end(), localVol.begin()));
            strikeGrid_.push_back(strikeGrid);
            locvolGrid_.push_back(localVol);
            sigmaS_.push_back(Linear().interpolate(strikeGrid_.back().begin(), strikeGrid_.back().end(), locvolGrid_.back().begin()));

            // set up swap rate model here such that we don't need to care about it during sigma_x calculation
            swapRateModel_ = qGSwapRateModel(swapRate.scf(), times()[idx - 1]);

            // simulate the next step T[idx-1] to T[idx]
            simulation_->simulate(idx + 1);
        }
    }

    void QGLocalvolModelForwardStochVolFlavor::simulateAndCalibrate() {
        QL_REQUIRE(volTS_->volatilityType() == VolatilityType::Normal, "Normal volatilities required.");
        Initialiser init(ext::static_pointer_cast<QGLocalvolModelBackwardFlavor>(shared_from_this()));
        simulation_->simulate(1);

        // run actual simulation and calibration
        for (size_t idx = 1; idx < times().size(); ++idx) {
            // specify swap rate fixing at T[idx], but observed at T[idx-1]
            SwapRate swapRate(this, init.today(), times()[idx]);

            // set up smile section and strike grid
            ext::shared_ptr<SmileSection> smileSection = volTS_->smileSection(times()[idx], swapIndex_->tenor(), true);  // the vanilla model only provides terminal distribution of swap rate at T[idx]; at T[idx-1] it is a different swap rate
            Real stdDev = smileSection->volatility(swapRate.swapRate()) * sqrt(times()[idx - 1]);  // approximated ATM forward vol
            std::vector<Real> initialStrikeGrid(stdDevGrid_);
            for (size_t k = 0; k < initialStrikeGrid.size(); ++k) initialStrikeGrid[k] = initialStrikeGrid[k] * stdDev + swapRate.swapRate();

            // these are the grids used for smile section calculation
            std::vector<Real> smileStrikeGrid;
            std::vector<Real> callGrid;
            std::vector<Real> d2CalldK2Grid;
            for (size_t k = 0; k < initialStrikeGrid.size(); ++k) {

                Real h = 1.0e-5; // 0.1bp shift size
                                 // basic idea is to approximate forward smile based on terminal smile from vanilla model
                Real call = bachelierBlackFormula(Option::Call, initialStrikeGrid[k], swapRate.swapRate(), smileSection->volatility(initialStrikeGrid[k])*sqrt(times()[idx - 1]));
                Real callph = bachelierBlackFormula(Option::Call, initialStrikeGrid[k] + h, swapRate.swapRate(), smileSection->volatility(initialStrikeGrid[k] + h)*sqrt(times()[idx - 1]));
                Real callmh = bachelierBlackFormula(Option::Call, initialStrikeGrid[k] - h, swapRate.swapRate(), smileSection->volatility(initialStrikeGrid[k] - h)*sqrt(times()[idx - 1]));

                // we only calculate local vol if within [1e-6, 1 - 1e-6] quantile 
                Real quantile = 1.0 + (callph - callmh) / 2.0 / h;
                Real d2CalldK2 = (callph + callmh - 2.0*call) / h / h;
                if (debugLevel_>1) debugLog_.push_back("T = " + std::to_string(times()[idx]) + ", swapRate = " + std::to_string(swapRate.swapRate()) + ", k = " + std::to_string(k) + ", strike = " + std::to_string(initialStrikeGrid[k]) + ", call = " + std::to_string(call) + ", callph = " + std::to_string(callph) + ", callmh = " + std::to_string(callmh) + ", quantile = " + std::to_string(quantile) + ", d2CalldK2 = " + std::to_string(d2CalldK2));
                if ((quantile < 1.0e-6) || (quantile > 1.0 - 1.e-6)) {
                    if (debugLevel_>1) debugLog_.push_back("Warning: Skip local vol.");
                    continue;  // skip calculation
                }
                if (d2CalldK2 < 1.0e-8) { // this check is not trivial since our forward vol/price approximation is not guarantied arb-free even if terminal smile is arb-free
                    if (debugLevel_>1) debugLog_.push_back("Warning: d2CalldK2 = " + std::to_string(d2CalldK2) + "Skip local vol.");
                    continue;  // skip calculation
                }
                // we need to re-calculate call price at the next grid point because this is needed for dCall /dT
                call = smileSection->optionPrice(initialStrikeGrid[k], Option::Call);

                smileStrikeGrid.push_back(initialStrikeGrid[k]);
                callGrid.push_back(call);
                d2CalldK2Grid.push_back(d2CalldK2);

                // done!
            }

            // calculate MC option prices
            Time obsTime = (idx>0) ? (times()[idx - 1]) : (0.0);  // this should always be T[idx-1] here
            McCalculator mcCalc(this, obsTime, swapRate.scf(), swapRate.annuity(), swapRate.swapRate(), smileStrikeGrid);

            std::vector<Real> strikeGrid;
            std::vector<Real> dCalldTGrid;
            std::vector<Real> localVol;
            // calculate dC/dT and setup final strikes and vols
            for (size_t j = 0; j < mcCalc.vanillaOptions().size(); ++j) {
                Real dCalldT = (callGrid[j] - mcCalc.vanillaOptions()[j]) / (times()[idx] - obsTime);
                if (dCalldT < 1.0e-8) {
                    if (debugLevel_>1) debugLog_.push_back("Warning: dCalldT = " + std::to_string(dCalldT) + ", Skip local vol.");
                    continue;  // skip calculation
                }
                Real lvol = sqrt(2 * dCalldT / d2CalldK2Grid[j]);
                strikeGrid.push_back(smileStrikeGrid[j]);
                dCalldTGrid.push_back(dCalldT);
                localVol.push_back(lvol);
            }

            // logging
            if (debugLevel_ > 0) debugLog_.push_back("T = " + std::to_string(times()[idx]) + ", fixingDate = " + std::to_string(swapRate.fixingDate().serialNumber()) + ", swapRate = " + std::to_string(swapRate.swapRate()) + ", annuity = " + std::to_string(swapRate.annuity()) + ", stdDev = " + std::to_string(stdDev) + ", avgCalcStrikes = " + std::to_string(mcCalc.avgCalcStrikes()));
            if (debugLevel_ > 1) {
                for (size_t j = 0; j < strikeGrid.size(); ++j) {
                    debugLog_.push_back("strike = " + std::to_string(strikeGrid[j]) + ", dCalldT = " + std::to_string(dCalldTGrid[j]) + ", localVol = " + std::to_string(localVol[j]));
                }
            }


            // calculate E^A[ z(T) | S(T) = K ] and adjust local vol
            StochvolExpectation expZ(this, idx, 1.0 / kernelWidth_ / stdDev, swapRate.annuity(), mcCalc, strikeGrid, &kernel);
            // finally adjust sigma_SV = sigma_LV / E^A[ z(T) | S(T) = K ]
            for (size_t j = 0; j < strikeGrid.size(); ++j) localVol[j] = localVol[j] / sqrt(expZ.expectationZCondS()[j]);

            // set up interpolation
            strikeGrid_.push_back(strikeGrid);
            locvolGrid_.push_back(localVol);
            sigmaS_.push_back(Linear().interpolate(strikeGrid_.back().begin(), strikeGrid_.back().end(), locvolGrid_.back().begin()));

            // set up swap rate model here such that we don't need to care about it during sigma_x calculation
            swapRateModel_ = qGSwapRateModel(swapRate.scf(), times()[idx - 1]);

            // simulate the next step T[idx-1] to T[idx]
            simulation_->simulate(idx + 1);

            // test calibration
            if (debugLevel_ > 2) checkMCPrices(times()[idx], swapRate.scf(), swapRate.annuity(), swapRate.swapRate(), smileStrikeGrid); // we might want to evaluate MC swaptions for debugging
        }
    }



}

