/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2018 Sebastian Schlenkrich

*/



#ifndef quantlib_qgcalibrator_hpp
#define quantlib_qgcalibrator_hpp

#include <ql/indexes/swapindex.hpp>
#include <ql/instruments/swaption.hpp>
#include <ql/termstructures/volatility/swaption/swaptionvolstructure.hpp>

#include <ql/math/optimization/constraint.hpp>
#include <ql/math/optimization/levenbergmarquardt.hpp>

#include <ql/experimental/basismodels/swaptioncfs.hpp>

#include <ql/experimental/templatemodels/qgaussian2/qgaverageswapratemodelT.hpp>


#define _MIN_( a, b ) ( (a) < (b) ? (a) : (b) )
#define _MAX_( a, b ) ( (a) > (b) ? (a) : (b) )

namespace QuantLib {

    typedef QuasiGaussianModel2T<QuantLib::Time, QuantLib::Real, QuantLib::Real> QuasiGaussianModel;
    typedef QGSwaprateModelT<QuantLib::Time, QuantLib::Real, QuantLib::Real> QGSwaprateModel;
    typedef QGAverageSwaprateModelT<QuantLib::Time, QuantLib::Real, QuantLib::Real> QGAverageSwaprateModel;


    // calibrate Quasi Gaussian model to implied normal volatilities
    class QGCalibrator { 

        // the Quasi-Gaussian model aimed to be calibrated
        ext::shared_ptr<QuasiGaussianModel> model_;

        // the resulting calibrated Quasi-Gaussian model
        ext::shared_ptr<QuasiGaussianModel> calibratedModel_;

        // we may use a MC simulation to adjust for closed form vs MC bias
        // ext::shared_ptr<RealMCSimulation> mcSimulation_;

        // calibration target volatilities
        Handle<SwaptionVolatilityStructure> volTS_;
        
        // we calibrate to strips of swaption volatilities; maybe also co-terminals can be relevant
        std::vector< ext::shared_ptr<SwapIndex> > swapIndices_;

        // parameters for quasi-Gaussian swaption model
        Real   modelTimesStepSize_;  // in year fractions
        bool   useExpectedXY_;

        // constraints for input model parameters
        Real  sigmaMin_, sigmaMax_, slopeMin_, slopeMax_, etaMin_, etaMax_;

        // we want to control calibration by individual weights
        Real sigmaWeight_, slopeWeight_, etaWeight_;

        // apply a penalty on differing slope and skew values per factor
        Real penaltySigma_, penaltySlope_;

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
            // reference to access model, swaptions and targets
            QGCalibrator  *calibrator_;
            // we work on a copy of the initial model (better safe then sorry)
            ext::shared_ptr<QuasiGaussianModel> model_;
            // specify parameters used for optimisation; dimension [modeltimes] x [ d (sigma) + d (slope) + 1 (eta) ]
            std::vector< std::vector< Real > > isInput_;
            // specify targets used for optimisation; dimension [exercises] x [ swapterm (SigmaATM) + swapterm (Skew) + swapterm (Smile) ]
            std::vector< std::vector< Real > > isOutput_;
            // count inputs and outputs
            Size inputSize_;
            Size outputSize_;

            // container for single curve swaption cash flow representation and calibration targets; initialised at inception
            class CalibSwaption : public SwapCashFlows {
                // cash flows are inherited from base class
                // targets
                Real sigmaATM_, skew_, smile_;
                // auxilliaries
                std::vector< Real > modelTimes_;
            public:
                CalibSwaption ( Date                                expiryDate,
                                const ext::shared_ptr<SwapIndex>&   swapindex,
                                const Handle<YieldTermStructure>&   discountCurve,
                                const Handle<SwaptionVolatilityStructure> volTS,
                                bool                                contTenorSpread = true,
                                Real                                modelTimesStepSize = 1.0/12.0);
                // inspectors
                inline Real sigmaATM() const { return sigmaATM_; }
                inline Real skew()     const { return skew_;     }
                inline Real smile()    const { return smile_;    }
                inline const std::vector<Real>& modelTimes() const { return modelTimes_; }
            };
            std::vector< std::vector< ext::shared_ptr<CalibSwaption> > > calibSwaptions_;
            
        public:
            const Size inputSize()  const { return inputSize_; }
            const Size outputSize() const { return outputSize_; }

            Objective( QGCalibrator                               *calibrator,
                       const std::vector< std::vector< Real > >&  isInput,
                       const std::vector< std::vector< Real > >&  isOutput );

            // initialize state X with model parameters and apply inverse transformation
            Array initialise();
                // allocate X with inputSize
                // for all isInput[i]
                //     for j=0.. d-1, isInput[i][j] ? sigma -> inverseSigma -> X[idx], ++idx
                //     for j=d..2d-1, isInput[i][j] ? slope -> inverseSlope -> X[idx], ++idx
                //     for j=2d,      isInput[i][j] ? eta   -> inverseEta   -> X[idx], ++idx

            // apply direct transformation and update model parameters
            void update(const Array& X) const;
                // for all isInput[i]
                //     for j=0.. d-1, isInput[i][j] ? X[idx] -> inverseSigma -> lambda, ++idx
                //     for j=d..2d-1, isInput[i][j] ? X[idx] -> inverseSlope -> b     , ++idx
                //     for j=2d,      isInput[i][j] ? X[idx] -> inverseEta   -> eta   , ++idx
                //  model->update(lambda,b,eta)

            // CostFunction interface
            virtual Array values(const Array& x) const;
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
        QGCalibrator( const ext::shared_ptr<QuasiGaussianModel>&          model,
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
                      const EndCriteria&                                    endCriteria );

        virtual ~QGCalibrator() = default;

        // a single optimisation run
        Integer calibrate( const std::vector< std::vector< Real > >&  isInput,
                           const std::vector< std::vector< Real > >&  isOutput,
                           Real                                       epsfcn = 1.0e-4 );  // delta for finite differences

        // inspectors
        inline const ext::shared_ptr<QuasiGaussianModel> calibratedModel() const { return calibratedModel_; }
        inline const std::vector<std::string>& debugLog() const { return debugLog_; }

        inline void acceptCalibration() { model_ = calibratedModel_; }

    };

}

#undef _MIN_
#undef _MAX_

#endif  /* ifndef quantlib_qgcalibrator_hpp */
