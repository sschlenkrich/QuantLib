/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2018 Sebastian Schlenkrich

*/



#ifndef quantlib_qglocalvolmodel_hpp
#define quantlib_qglocalvolmodel_hpp

#include <ql/math/interpolation.hpp>

#include <ql/indexes/swapindex.hpp>
#include <ql/instruments/swaption.hpp>
#include <ql/termstructures/volatility/swaption/swaptionvolstructure.hpp>

//#include <ql/experimental/templatemodels/auxilliaries/auxilliariesT.hpp>

#include <ql/experimental/basismodels/swaptioncfs.hpp>

#include <ql/experimental/templatemodels/qgaussian2/quasigaussianmodel2T.hpp>
#include <ql/experimental/templatemodels/qgaussian2/qgswapratemodelT.hpp>
#include <ql/experimental/templatemodels/montecarlo/mcsimulationT.hpp>
#include <ql/experimental/templatemodels/montecarlo/ratespayoffT.hpp>


namespace QuantLib {

    // calibrate Quasi Gaussian model to implied normal volatilities
    class QGLocalvolModel : public QuasiGaussianModel2T<QuantLib::Time, QuantLib::Real, QuantLib::Real> { 
    public:
        typedef QuasiGaussianModel2T<QuantLib::Time, QuantLib::Real, QuantLib::Real> QuasiGaussianModel;
        typedef QGSwaprateModelT<QuantLib::Time, QuantLib::Real, QuantLib::Real> QGSwaprateModel;
        typedef MCSimulationT<QuantLib::Time,QuantLib::Real,QuantLib::Real> MCSimulation;    
        typedef MCPayoffT<QuantLib::Time, QuantLib::Real, QuantLib::Real> MCPayoff;
        typedef BasePayoffT<QuantLib::Time, QuantLib::Real, QuantLib::Real> MCBase;
        typedef RatesPayoffT<QuantLib::Time, QuantLib::Real, QuantLib::Real>::Annuity MCAnnuity;
        typedef RatesPayoffT<QuantLib::Time, QuantLib::Real, QuantLib::Real>::GeneralSwaption MCSwaption;

    protected:
        // calibration target volatilities
        Handle<SwaptionVolatilityStructure> volTS_;
        
        // we calibrate to a strip of swaption volatilities; maybe also co-terminals can be relevant
        ext::shared_ptr<SwapIndex> swapIndex_;

        // we describe local volatility sigmaS as a set of 1-D interpolations per time step 1 to N (excluding 0)
        std::vector< Interpolation > sigmaS_;
        std::vector < std::vector<Real> > strikeGrid_;  // we need to separately store data for interpolation
        std::vector < std::vector<Real> > locvolGrid_;

        // we have three modes for sigma_x calculation (during and after calibration phase)
        enum { 
            Parent,         // delegate sigma_x calculation to parent class; this typically doesn't do something meaningfull
            Calibration,    // calculate sigma_x as specified for calibration procedure
            Pricing         // use calibrated local vol after it is fully calibrated
        } sigmaMode_;

        // we specify the local vol grid in terms of standard deviations
        std::vector<Real> stdDevGrid_;

        // simulation derived during calibration; we need to store nPaths_ and seed_ separately
        // because we can't initialise simulation within constructor
        size_t     nPaths_;
        BigNatural seed_;
        ext::shared_ptr<MCSimulation> simulation_;
    
        // stochastic volatility uses kernel-based conditional expectation evaluation
        bool calcStochVolAdjustment_;  // switch adjustment on/off independent of vol-of-vol marking
        Real kernelWidth_; // kernel width for conditional expectation calculation in terms of stdDevs
        inline static Real kernel(const Real u) { return ((fabs(u) < 1.0) ? (1.0 - fabs(u)) : (0.0)); }

        // debugging, warning and errors
        std::vector<std::string> debugLog_;
        size_t debugLevel_; // 0 ... no debugging
                            // 1 ... debugging for time steps
                            // 2 ... debugging for strikes
                            // 3 ... warnings for simulation
                            // 4 ... debugging for each path (not recommended)

        // we need to cache the swap rate model and observation time for the current time step
        ext::shared_ptr<QGSwaprateModel> swapRateModel_;

        // we need a factory for MC swaptions for efficient calculation
        class SwaptionFactory {
            ext::shared_ptr<MCPayoff> floatLeg_;
            ext::shared_ptr<MCPayoff> annuityLeg_;
        public:
            SwaptionFactory(const Time obsTime, const SwapCashFlows& scf);
            ext::shared_ptr<MCPayoff> swaption(const Real strike, const Real CallOrPut);
        };

        // some helpers and subroutines to structure the methodologies

        // find idx = min_j{ x <= X[j] }
        inline static size_t minIdx(const std::vector<Real>& X, const Real x);

        ext::shared_ptr<QGSwaprateModel> qGSwapRateModel(const SwapCashFlows& scf, const Real obsTime);

        // encapsulate initial set up
        class Initialiser {
        private:
            Date today_;
        public:
            Initialiser(ext::shared_ptr<QGLocalvolModel> model);
            const Date& today() const { return today_; }
        };

        // set up swap rate and provide corresponding quantities
        class SwapRate {
            Real swapRate_, annuity_;
            Date fixingDate_;
            SwapCashFlows scf_;
        public:
            SwapRate(const QGLocalvolModel *model, const Date today, const Real fixingTime);
            const Real  swapRate()   const { return swapRate_;   }
            const Real  annuity()    const { return annuity_;    }
            const Date& fixingDate() const { return fixingDate_; }
            const SwapCashFlows& scf() const { return scf_; }
        };

        // encapsulate the MC calculation
        class McCalculator {
        private:
            std::vector<Real>  oneOverBSample_;  // assume size simulation_->nPaths()
            std::vector<Real>  annuitySample_;   // assume size simulation_->nPaths()
            std::vector<Real>  swapRateSample_;  // assume size simulation_->nPaths()
            std::vector<Real>  vanillaOptions_;  // assume vanillaOptions(smileStrikeGrid.size(), 0.0);
            Real               avgCalcStrikes_;
        public:
            McCalculator(QGLocalvolModel          *model,
                         const Real               obsTime,
                         const SwapCashFlows&     scf,
                         const Real               annuity,
                         const Real               swapRate,
                         const std::vector<Real>& smileStrikeGrid);
            const std::vector<Real>&  oneOverBSample() const { return oneOverBSample_; }
            const std::vector<Real>&  annuitySample()  const { return annuitySample_;  }
            const std::vector<Real>&  swapRateSample() const { return swapRateSample_; }
            const std::vector<Real>&  vanillaOptions() const { return vanillaOptions_; }
            const Real                avgCalcStrikes() const { return avgCalcStrikes_; }
        };

        // encapsulate stoch vol adjustment; calculate E^A[ z(T) | S(T) = K ] 
        class StochvolExpectation {
        private:
            std::vector<Real>  expectationZCondS_;
        public:
            StochvolExpectation(const QGLocalvolModel     *model,
                                const size_t              simIdx,
                                const Real                lambda,   // lambda = 1.0 / kernelWidth / stdDev
                                const Real                annuity,
                                const McCalculator&       mcCalc,
                                const std::vector<Real>&  strikeGrid,
                                Real                      (*kernel)(const Real) );
            const std::vector<Real>&  expectationZCondS() const { return expectationZCondS_; }
        };

        inline void checkMCPrices( const Real obsTime,
                                   const SwapCashFlows&     scf,
                                   const Real               annuity,
                                   const Real               swapRate,
                                   const std::vector<Real>& smileStrikeGrid );

    public:
        // pure local volatility
        QGLocalvolModel( const Handle<YieldTermStructure>&                      termStructure,
                         const Handle<SwaptionVolatilityStructure>&             volTS,
                         const Real                                             chi,
                         const ext::shared_ptr<SwapIndex>&                    swapIndex,
                         const std::vector<Real>&                               times,
                         const std::vector<Real>&                               stdDevGrid,
                         const size_t                                           nPaths,
                         const BigNatural                                       seed = 1234,
                         const size_t                                           debugLevel = 1);

        // allow for stochastic volatility
        QGLocalvolModel(const Handle<YieldTermStructure>&                      termStructure,
                        const Handle<SwaptionVolatilityStructure>&             volTS,
                        const Real                                             chi,
                        const Real                                             theta,
                        const Real                                             eta,
                        const ext::shared_ptr<SwapIndex>&                    swapIndex,
                        const std::vector<Real>&                               times,
                        const std::vector<Real>&                               stdDevGrid,
                        const bool                                             calcStochVolAdjustment,
                        const Real                                             kernelWidth,
                        const size_t                                           nPaths,
                        const BigNatural                                       seed = 1234,
                        const size_t                                           debugLevel = 1);

        virtual ~QGLocalvolModel() = default;

        // do the actual calculation
        virtual void simulateAndCalibrate() = 0;

        virtual std::vector< std::vector<Real> > sigma_xT(const Real t, const State& s);

        // inspectors
        inline const ext::shared_ptr<MCSimulation> simulation() { return simulation_; }
        inline const Real sigmaS(const size_t idx, const Real s) { return sigmaS_[idx](s); }

        std::vector<std::string> debugLog() { return debugLog_; }

        // test the calibration of the model
        std::vector< std::vector<Real> > calibrationTest(const std::vector<Date>&  exerciseDates,
                                                         const std::vector<Real>&  stdDevStrikes );
    };

    // alternative backward-looking calibration methodology for QGLocalvolModel
    class QGLocalvolModelBackwardFlavor : public QGLocalvolModel {
    public:
        QGLocalvolModelBackwardFlavor(
            const Handle<YieldTermStructure>&                      termStructure,
            const Handle<SwaptionVolatilityStructure>&             volTS,
            const Real                                             chi,
            const Real                                             theta,
            const Real                                             eta,
            const ext::shared_ptr<SwapIndex>&                    swapIndex,
            const std::vector<Real>&                               times,
            const std::vector<Real>&                               stdDevGrid,
            const bool                                             calcStochVolAdjustment,
            const Real                                             kernelWidth,
            const size_t                                           nPaths,
            const BigNatural                                       seed = 1234,
            const size_t                                           debugLevel = 1)
            : QGLocalvolModel(termStructure, volTS, chi, theta, eta, swapIndex, times, stdDevGrid, calcStochVolAdjustment, kernelWidth, nPaths, seed, debugLevel) {}

        virtual ~QGLocalvolModelBackwardFlavor() = default;

        // do the actual calculation
        virtual void simulateAndCalibrate();
    };

    // alternative forward-looking calibration methodology for QGLocalvolModel
    class QGLocalvolModelForwardFlavor : public QGLocalvolModel {
    public:
        QGLocalvolModelForwardFlavor(
            const Handle<YieldTermStructure>&                      termStructure,
            const Handle<SwaptionVolatilityStructure>&             volTS,
            const Real                                             chi,
            const Real                                             theta,
            const Real                                             eta,
            const ext::shared_ptr<SwapIndex>&                    swapIndex,
            const std::vector<Real>&                               times,
            const std::vector<Real>&                               stdDevGrid,
            const bool                                             calcStochVolAdjustment,
            const Real                                             kernelWidth,
            const size_t                                           nPaths,
            const BigNatural                                       seed = 1234,
            const size_t                                           debugLevel = 1)
            : QGLocalvolModel(termStructure, volTS, chi, theta, eta, swapIndex, times, stdDevGrid, calcStochVolAdjustment, kernelWidth, nPaths, seed, debugLevel) {}

        virtual ~QGLocalvolModelForwardFlavor() = default;

        // do the actual calculation
        virtual void simulateAndCalibrate();
    };

    // alternative analytical calibration methodology for QGLocalvolModel
    class QGLocalvolModelAnalyticFlavor : public QGLocalvolModel {
    public:
        QGLocalvolModelAnalyticFlavor(
            const Handle<YieldTermStructure>&                      termStructure,
            const Handle<SwaptionVolatilityStructure>&             volTS,
            const Real                                             chi,
            const ext::shared_ptr<SwapIndex>&                    swapIndex,
            const std::vector<Real>&                               times,
            const std::vector<Real>&                               stdDevGrid,
            const size_t                                           nPaths,
            const BigNatural                                       seed = 1234,
            const size_t                                           debugLevel = 1)
            : QGLocalvolModel(termStructure, volTS, chi, swapIndex, times, stdDevGrid, nPaths, seed, debugLevel) {}

        virtual ~QGLocalvolModelAnalyticFlavor() = default;

        // do the actual calculation
        virtual void simulateAndCalibrate();
    };

    // alternative forward-looking calibration methodology for QGLocalvolModel
    class QGLocalvolModelMonteCarloFlavor : public QGLocalvolModel {
    public:
        QGLocalvolModelMonteCarloFlavor(
            const Handle<YieldTermStructure>&                      termStructure,
            const Handle<SwaptionVolatilityStructure>&             volTS,
            const Real                                             chi,
            const Real                                             theta,
            const Real                                             eta,
            const ext::shared_ptr<SwapIndex>&                    swapIndex,
            const std::vector<Real>&                               times,
            const std::vector<Real>&                               stdDevGrid,
            const bool                                             calcStochVolAdjustment,
            const Real                                             kernelWidth,
            const size_t                                           nPaths,
            const BigNatural                                       seed = 1234,
            const size_t                                           debugLevel = 1)
            : QGLocalvolModel(termStructure, volTS, chi, theta, eta, swapIndex, times, stdDevGrid, calcStochVolAdjustment, kernelWidth, nPaths, seed, debugLevel) {}

        virtual ~QGLocalvolModelMonteCarloFlavor() = default;

        // do the actual calculation
        virtual void simulateAndCalibrate();
    };

    // alternative forward-looking calibration methodology for QGLocalvolModel
    class QGLocalvolModelForwardStochVolFlavor : public QGLocalvolModel {
    private:
    public:
        QGLocalvolModelForwardStochVolFlavor(
            const Handle<YieldTermStructure>&                      termStructure,
            const Handle<SwaptionVolatilityStructure>&             volTS,
            const Real                                             chi,
            const Real                                             theta,
            const Real                                             eta,
            const ext::shared_ptr<SwapIndex>&                    swapIndex,
            const std::vector<Real>&                               times,
            const std::vector<Real>&                               stdDevGrid,
            const Real                                             kernelWidth,
            const size_t                                           nPaths,
            const BigNatural                                       seed = 1234,
            const size_t                                           debugLevel = 1)
            : QGLocalvolModel(termStructure, volTS, chi, theta, eta, swapIndex, times, stdDevGrid, true, kernelWidth, nPaths, seed, debugLevel) {}

        virtual ~QGLocalvolModelForwardStochVolFlavor() = default;

        // do the actual calculation
        virtual void simulateAndCalibrate();
    };



}

#endif  /* ifndef quantlib_qglocalvolmodel_hpp */
