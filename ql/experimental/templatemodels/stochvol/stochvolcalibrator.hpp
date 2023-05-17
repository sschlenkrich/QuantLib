/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2010 Sebastian Schlenkrich

*/



#ifndef quantlib_stochvolcalibrator_hpp
#define quantlib_stochvolcalibrator_hpp

#include <ql/math/optimization/constraint.hpp>
#include <ql/math/optimization/levenbergmarquardt.hpp>

#include <ql/experimental/templatemodels/stochvol/stochvolmodels.hpp>


#define _MIN_( a, b ) ( (a) < (b) ? (a) : (b) )
#define _MAX_( a, b ) ( (a) > (b) ? (a) : (b) )

namespace QuantLib {

    // calibrate stochastic volatility model to implied normal volatilities
    class StochVolModelCalibrator : public CostFunction { 
    private:
        // model parameters
        Real  lambda_;
        Real  b_;
        Real  L_;
        Real  theta_;
        Real  m_;
        Real  eta_;
        Real  z0_;
        Real  rho_;
        // calibration flags
        bool  lambdaIsFixed_;
        bool  bIsFixed_;
        bool  LIsFixed_;
        bool  thetaIsFixed_;
        bool  mIsFixed_;
        bool  etaIsFixed_;
        bool  z0IsFixed_;
        bool  rhoIsFixed_;
        // lower boundaries
        Real  lambdaMin_;
        Real  bMin_;
        Real  LMin_;
        Real  thetaMin_;
        Real  mMin_;
        Real  etaMin_;
        Real  z0Min_;
        Real  rhoMin_;
        // upper boundaries
        Real  lambdaMax_;
        Real  bMax_;
        Real  LMax_;
        Real  thetaMax_;
        Real  mMax_;
        Real  etaMax_;
        Real  z0Max_;
        Real  rhoMax_;
        // optimization parameters
        Real epsfcn_, ftol_, xtol_, gtol_, glAbsAcc_;
        Size maxfev_, glMaxEval_;
        // calibration targets
        Real              exercTime_;
        Real              forward_;
        std::vector<Real> strikes_;
        std::vector<Real> vols_;
    public:
        // parameter transformation
        Real  direct(const Real x, const Real a, const Real b) const;
        Real  inverse(const Real y, const Real a, const Real b) const;
        Array direct(const Array& X) const;
        Array inverse(const Array& Y) const;
        // initialize state X with model parameters and apply inverse transformation
        Array initialise();
        // CostFunction interface
        virtual Array values(const Array& x) const;
        virtual Real value(const Array& x) const;

        // inspectors
        const Real  lambda() const { return lambda_; }
        const Real  b()		 const { return b_;		 }
        const Real  L()		 const { return L_;		 }
        const Real  theta()	 const { return theta_;	 }
        const Real  m()		 const { return m_;		 }
        const Real  eta()	 const { return eta_;	 }
        const Real  z0()	 const { return z0_;	 }
        const Real  rho()	 const { return rho_;	 }

        // factory
        ext::shared_ptr<RealStochVolModel> model() {
            return ext::shared_ptr<RealStochVolModel>(new RealStochVolModel(lambda_,b_,L_,theta_,m_,eta_,z0_,rho_));
        }
        ext::shared_ptr<RealStochVolModel> model(const Array& X) const;
        // count number of effective states
        size_t xDim() const;
        // update model parameters with state X
        void update(const Array& X);

        // construct object and calibrate
        StochVolModelCalibrator( // initial model parameters
                                 const Real  lambda,
                                 const Real  b,
                                 const Real  L,
                                 const Real  theta,
                                 const Real  m,
                                 const Real  eta,
                                 const Real  z0,
                                 const Real  rho,
                                 // calibration flags
                                 const bool  lambdaIsFixed,
                                 const bool  bIsFixed,
                                 const bool  LIsFixed,
                                 const bool  thetaIsFixed,
                                 const bool  mIsFixed,
                                 const bool  etaIsFixed,
                                 const bool  z0IsFixed,
                                 const bool  rhoIsFixed,
                                 // calibration targets
                                 const Real  exercTime,
                                 const Real  forward,
                                 const std::vector<Real>& strikes,
                                 const std::vector<Real>& vols,
                                 const std::vector<Real>& optimizationParams  // { [min], [max], epsfcn, ftol, xtol, gtol, maxfev, glAbsAcc, glMaxEval }							 
                                 );

    };

}

#undef _MIN_
#undef _MAX_

#endif  /* ifndef quantlib_stochvolcalibrator_hpp */
