/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2010 Sebastian Schlenkrich

*/



#ifndef quantlib_qgaussiancalibrator_hpp
#define quantlib_qgaussiancalibrator_hpp

#include <ql/instruments/swaption.hpp>

#include <ql/math/optimization/constraint.hpp>
#include <ql/math/optimization/levenbergmarquardt.hpp>

#include <ql/experimental/templatemodels/auxilliaries/auxilliariesT.hpp>

#include <ql/experimental/basismodels/swaptioncfs.hpp>

#include <ql/experimental/templatemodels/qgaussian/quasigaussianmodels.hpp>
#include <ql/experimental/templatemodels/montecarlo/montecarlomodells.hpp>


#define _MIN_( a, b ) ( (a) < (b) ? (a) : (b) )
#define _MAX_( a, b ) ( (a) > (b) ? (a) : (b) )

namespace QuantLib {


    // calibrate stochastic volatility model to implied normal volatilities
    class QuasiGaussianModelCalibrator { 

        // the Quasi-Gaussian model aimed to be calibrated
        ext::shared_ptr<RealQuasiGaussianModel> model_;

        // the resulting calibrated Quasi-Gaussian model
        ext::shared_ptr<RealQuasiGaussianModel> calibratedModel_;

        // we may use a MC simulation to adjust for closed form vs MC bias
        ext::shared_ptr<RealMCSimulation> mcSimulation_;

        // calibration targets

        // grid of reference (ATM) swaptions as calibration instruments, swaptions_[exercise][swapterm]
        std::vector< std::vector< ext::shared_ptr<Swaption> > > swaptions_;

        // grid of reference swaption stoch vol model parameters
        // user supplied or calculated
        // dimensions according to dimensions of swaptions_
        std::vector< std::vector< Real > > lambda_;
        std::vector< std::vector< Real > > b_;
        std::vector< std::vector< Real > > eta_;

        // constraints for input model parameters
        Real  lambdaMin_, lambdaMax_, bMin_, bMax_, etaMin_, etaMax_;
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

        // parameters for quasi-Gaussian swaption model
        Real   modelTimesStepSize_;  // in year fractions
        bool   useExpectedXY_;

        class Objective : public CostFunction {
        private:
            // reference to access model, swaptions and targets
            QuasiGaussianModelCalibrator  *calibrator_;
            // we work on a copy of the initial model (better safe then sorry)
            ext::shared_ptr<RealQuasiGaussianModel> model_;
            // specify parameters used for optimisation; dimension [modeltimes] x [ d (lambda) + d (b) + 1 (eta) ]
            std::vector< std::vector< Real > > isInput_;
            // specify targets used for optimisation; dimension [exercises] x [ swapterm (lambda) + swapterm (b) + swapterm (eta) ]
            std::vector< std::vector< Real > > isOutput_;
            // count inputs and outputs
            Size inputSize_;
            Size outputSize_;

            // container for single curve swaption cash flow representation and calibration targets; initialised at inception
            class CalibSwaption : public SwaptionCashFlows {
                // cash flows are inherited from base class
                // targets
                Real lambda_, b_, eta_;
                // auxilliaries
                std::vector< Real > modelTimes_;
            public:
                CalibSwaption ( Real                               lambda,
                                Real                               b,
                                Real                               eta,
                                const ext::shared_ptr<Swaption>& swaption,
                                const Handle<YieldTermStructure>&  discountCurve,
                                bool                               contTenorSpread = true,
                                Real                               modelTimesStepSize = 1.0/12.0);
                // inspectors
                inline Real lambda() const { return lambda_; }
                inline Real b()      const { return b_;      }
                inline Real eta()    const { return eta_;    }
                inline const std::vector<Real>& modelTimes() const { return modelTimes_; }
            };
            std::vector< std::vector< ext::shared_ptr<CalibSwaption> > > calibSwaptions_;
            
        public:
            const Size inputSize()  const { return inputSize_; }
            const Size outputSize() const { return outputSize_; }

            Objective( QuasiGaussianModelCalibrator              *calibrator,
                       const std::vector< std::vector< Real > >&  isInput,
                       const std::vector< std::vector< Real > >&  isOutput );

            // initialize state X with model parameters and apply inverse transformation
            Array initialise();
                // allocate X with inputSize
                // for all isInput[i]
                //     for j=0.. d-1, isInput[i][j] ? lambda -> inverseLambda -> X[idx], ++idx
                //     for j=d..2d-1, isInput[i][j] ? b      -> inverseB      -> X[idx], ++idx
                //     for j=2d,      isInput[i][j] ? eta    -> inverseEta    -> X[idx], ++idx

            // apply direct transformation and update model parameters
            void update(const Array& X) const;
                // for all isInput[i]
                //     for j=0.. d-1, isInput[i][j] ? X[idx] -> inverseLambda -> lambda, ++idx
                //     for j=d..2d-1, isInput[i][j] ? X[idx] -> inverseB      -> b     , ++idx
                //     for j=2d,      isInput[i][j] ? X[idx] -> inverseEta    -> eta   , ++idx
                //  model->update(lambda,b,eta)

            // CostFunction interface
            virtual Array values(const Array& x) const;
                // update(x)
                // allocate Y with outputSize()
                // for k=0..swaptions_.size() construct RealQGSwaptionModel swpModel[k]
                // for all isOutput[i], j=0..swapterm-1
                //     isOutput[i][0*swapterm+j] ? swpModel[k] -> avLambda -> Y[idx], Y[idx] -= lambda[k], ++idx
                //     isOutput[i][1*swapterm+j] ? swpModel[k] -> avB      -> Y[idx], Y[idx] -= b[k],      ++idx 
                //     isOutput[i][2*swapterm+j] ? swpModel[k] -> avEta    -> Y[idx], Y[idx] -= eta[k],    ++idx
                //     ++k

            // 0.5 ||Y||^2
            virtual Real value(const Array& x) const;

            // return the model for a given input state
            const ext::shared_ptr<RealQuasiGaussianModel> model(const Array& x);
        };

    public:

        // constructor
        QuasiGaussianModelCalibrator( ext::shared_ptr<RealQuasiGaussianModel> model,
                                      ext::shared_ptr<RealMCSimulation>       mcSimulation,
                                      std::vector< std::vector< ext::shared_ptr<Swaption> > > swaptions,
                                      std::vector< std::vector< Real > > lambda,
                                      std::vector< std::vector< Real > > b,
                                      std::vector< std::vector< Real > > eta,
                                      Real                               lambdaMin,
                                      Real                               lambdaMax,
                                      Real                               bMin,
                                      Real                               bMax,
                                      Real                               etaMin,
                                      Real                               etaMax,
                                      Real                               modelTimesStepSize,
                                      bool                               useExpectedXY);

        Integer calibrate( const std::vector< std::vector< Real > >&  isInput,
                           const std::vector< std::vector< Real > >&  isOutput,
                           // optimization parameters
                           Real                                       epsfcn = 1.0e-10,
                           Real                                       ftol   = 1.0e-8,
                           Real                                       xtol   = 1.0e-8,
                           Real                                       gtol   = 1.0e-8,
                           Size                                       maxfev = 10000    );

        // inspectors
        inline const ext::shared_ptr<RealQuasiGaussianModel> calibratedModel() const { return calibratedModel_; }

        inline void acceptCalibration() { model_ = calibratedModel_; }

        // inline const ext::shared_ptr<RealQuasiGaussianModel> model()        { return model_;        }
        // inline const ext::shared_ptr<RealMCSimulation>       mcSimulation() { return mcSimulation_; }
        // inline const std::vector< std::vector< ext::shared_ptr<Swaption> > >& swaptions() { return swaptions_; }
        // inline const std::vector< std::vector< Real > >& lambda() { return lambda_; }
        // inline const std::vector< std::vector< Real > >& b()      { return b_;      }
        // inline const std::vector< std::vector< Real > >& eta()    { return eta_;    }

    };

}

#undef _MIN_
#undef _MAX_

#endif  /* ifndef quantlib_qgaussiancalibrator_hpp */
