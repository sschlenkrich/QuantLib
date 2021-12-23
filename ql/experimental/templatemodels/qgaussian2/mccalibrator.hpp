/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2019 Sebastian Schlenkrich

*/



#ifndef quantlib_mccalibrator_hpp
#define quantlib_mccalibrator_hpp

#include <ql/indexes/swapindex.hpp>
#include <ql/instruments/swaption.hpp>
#include <ql/termstructures/volatility/swaption/swaptionvolstructure.hpp>

#include <ql/math/optimization/constraint.hpp>
#include <ql/math/optimization/levenbergmarquardt.hpp>

#include <ql/experimental/basismodels/swaptioncfs.hpp>

#include <ql/experimental/templatemodels/qgaussian2/quasigaussianmodel2T.hpp>
#include <ql/experimental/templatemodels/montecarlo/mcsimulationT.hpp>
#include <ql/experimental/templatemodels/montecarlo/ratespayoffT.hpp>


namespace QuantLib {


    // calibrate Quasi Gaussian model to implied normal volatilities
    class QGMonteCarloCalibrator { 

        typedef QuasiGaussianModel2T<QuantLib::Time, QuantLib::Real, QuantLib::Real> QuasiGaussianModel;
        typedef MCSimulationT<QuantLib::Time, QuantLib::Real, QuantLib::Real> MCSimulation;
        typedef MCPayoffT<QuantLib::Time, QuantLib::Real, QuantLib::Real> MCPayoff;
        typedef BasePayoffT<QuantLib::Time, QuantLib::Real, QuantLib::Real> MCBase;
        typedef RatesPayoffT<QuantLib::Time, QuantLib::Real, QuantLib::Real>::Annuity MCAnnuity;
        typedef RatesPayoffT<QuantLib::Time, QuantLib::Real, QuantLib::Real>::GeneralSwaption MCSwaption;

        // the Quasi-Gaussian model aimed to be calibrated
        ext::shared_ptr<QuasiGaussianModel> model_;

        // the Monte-Carlo simulation used for calibration
        ext::shared_ptr<MCSimulation> mcSimulation_;

        // the resulting calibrated Quasi-Gaussian model
        ext::shared_ptr<QuasiGaussianModel> calibratedModel_;

        // calibration target volatilities
        Handle<SwaptionVolatilityStructure> volTS_;
        
        // we calibrate to strips of swaption volatilities; maybe also co-terminals can be relevant
        std::vector< ext::shared_ptr<SwapIndex> > swapIndices_;

        // constraints for input model parameters
        Real  sigmaMin_, sigmaMax_, slopeMin_, slopeMax_, curveMin_, curveMax_;

        // we want to control calibration by individual weights
        Real sigmaWeight_, slopeWeight_, curveWeight_;

        // apply a penalty on differing slope and skew values per factor
        Real penaltySigma_, penaltySlope_, penaltyCurve_;

        // transformation (-inf, +inf) -> (a, b)
        static const Real direct(const Real x, const Real a, const Real b) {
            //return (b-a)*(atan(x)/M_PI + 0.5) + a;
            return TemplateAuxilliaries::direct(x, a, b);
        }

        // transformation (a, b) -> (-inf, +inf)
        static const Real inverse(const Real y, const Real a, const Real b) {
            //return tan( ((y-a)/(b-a)-0.5) * M_PI );
            return TemplateAuxilliaries::inverse(y, a, b);
        }

        // we need to know when to stop iterating
        EndCriteria endCriteria_;

        // we do some logging for degugging purposes
        std::vector< std::string > debugLog_;

        // the objective function controls the actual calibration step
        class Objective : public CostFunction {
        private:
            // reference to access model, simulation, swaptions and targets
            QGMonteCarloCalibrator  *calibrator_;
            // specify parameters used for optimisation; dimension [modeltimes] x [ d (sigma) + d (slope) + d (curve) ]
            std::vector< std::vector< Real > > isInput_;
            // specify targets used for optimisation; dimension [exercises] x [ swapterm (SigmaATM) + swapterm (Skew) + swapterm (Smile) ]
            std::vector< std::vector< Real > > isOutput_;
            // count inputs and outputs
            Size inputSize_;
            Size outputSize_;

            Size firstSimulationIdx_;  // the first simulationTime we simulate to
            Size lastSimulationIdx_;   // the last simulationTime we simulate to

            // container for single curve swaption cash flow representation and calibration targets; initialised at inception
            class CalibSwaption : public SwapCashFlows {
                // cash flows are inherited from base class
                // swap rate details
                Real expiryTime_, S0_, annuity_;
                // targets
                Real atmCall_, highCall_, lowPut_;
                // auxilliaries
                Real sigmaATM_;
                Real atmVega_, highVega_, lowVega_;
            public:
                CalibSwaption ( Date                                expiryDate,
                                const ext::shared_ptr<SwapIndex>& swapindex,
                                const Handle<YieldTermStructure>&   discountCurve,
                                const Handle<SwaptionVolatilityStructure> volTS,
                                bool                                contTenorSpread = true);
                // inspectors
                inline Real expiryTime() const { return expiryTime_; }
                inline Real S0()       const { return S0_; }
                inline Real annuity()  const { return annuity_; }

                inline Real atmCall()  const { return atmCall_;  }
                inline Real highCall() const { return highCall_; }
                inline Real lowPut()   const { return lowPut_;   }

                inline Real atmVega()  const { return atmVega_;  }
                inline Real highVega() const { return highVega_; }
                inline Real lowVega()  const { return lowVega_;  }

                inline Real sigmaATM() const { return sigmaATM_; }
            };
            std::vector< std::vector< ext::shared_ptr<CalibSwaption> > > calibSwaptions_;
            
        public:
            const Size inputSize()  const { return inputSize_; }
            const Size outputSize() const { return outputSize_; }

            Objective( QGMonteCarloCalibrator                    *calibrator,
                       const std::vector< std::vector< Real > >&  isInput,
                       const std::vector< std::vector< Real > >&  isOutput );

            // initialize state X with model parameters and apply inverse transformation
            Array initialise();
                // allocate X with inputSize
                // for all isInput[i]
                //     for j=0.. d-1, isInput[i][j] ? sigma -> inverseSigma -> X[idx], ++idx
                //     for j=d..2d-1, isInput[i][j] ? slope -> inverseSlope -> X[idx], ++idx
                //     for j=2d..3d-1,isInput[i][j] ? curve -> inverseCurve -> X[idx], ++idx

            // apply direct transformation and update model parameters
            void update(const Array& X) const;
                // for all isInput[i]
                //     for j=0.. d-1, isInput[i][j] ? X[idx] -> inverseSigma -> sigma, ++idx
                //     for j=d..2d-1, isInput[i][j] ? X[idx] -> inverseSlope -> slope, ++idx
                //     for j=2d..3d-1,isInput[i][j] ? X[idx] -> inverseCurve -> curve, ++idx
                //  model->update(sigma,slope,curve)

            // CostFunction interface
            virtual Disposable<Array> values(const Array& x) const;
                // update(x)
                // allocate Y with outputSize()
                // for k=0..swaptions_.size() construct QGSwaprateModel swpModel[k]
                // for all isOutput[i], j=0..swapterm-1
                //     isOutput[i][0*swapterm+j] ? swpModel[k] -> sigmaATM -> Y[idx], Y[idx] -= sigmaATM[k], ++idx
                //     isOutput[i][1*swapterm+j] ? swpModel[k] -> skew     -> Y[idx], Y[idx] -= skew[k],     ++idx 
                //     isOutput[i][2*swapterm+j] ? swpModel[k] -> smile    -> Y[idx], Y[idx] -= smile[k],    ++idx
                //     ++k

            // 0.5 ||Y||^2
            virtual Real value(const Array& x) const;

            // return the model for a given input state
            const ext::shared_ptr<QuasiGaussianModel> model(const Array& x);
        };

    public:

        // constructor
        QGMonteCarloCalibrator(
                      const ext::shared_ptr<QGMonteCarloCalibrator::QuasiGaussianModel>&   model,
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
                      const EndCriteria&                                    endCriteria );

        virtual ~QGMonteCarloCalibrator() = default;

        // a single optimisation run
        Integer calibrate( const std::vector< std::vector< Real > >&  isInput,
                           const std::vector< std::vector< Real > >&  isOutput,
                           Real                                       epsfcn = 1.0e-4 );  // delta for finite differences

        // inspectors
        inline const ext::shared_ptr<QuasiGaussianModel> calibratedModel() const { return calibratedModel_; }
        inline const ext::shared_ptr<MCSimulation> mcSimulation() const { return mcSimulation_; }
        inline const std::vector<std::string>& debugLog() const { return debugLog_; }

        inline void acceptCalibration() { model_ = calibratedModel_; }

    };

}


#endif  /* ifndef quantlib_mccalibrator_hpp */
