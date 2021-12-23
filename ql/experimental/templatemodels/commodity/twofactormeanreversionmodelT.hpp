/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2015 Sebastian Schlenkrich

*/



#ifndef quantlib_template2fmeanreversionmodel_hpp
#define quantlib_template2fmeanreversionmodel_hpp

#include <ql/shared_ptr.hpp>
#include <ql/errors.hpp>
#include <ql/experimental/templatemodels/auxilliaries/auxilliariesT.hpp>
#include <ql/experimental/templatemodels/auxilliaries/integratorsT.hpp>
#include <ql/experimental/templatemodels/stochasticprocessT.hpp>
#include <ql/experimental/templatemodels/commodity/indextermstructure.hpp>



#define _MIN_( a, b ) ( (a) < (b) ? (a) : (b) )
#define _MAX_( a, b ) ( (a) > (b) ? (a) : (b) )

namespace QuantLib {

    // 2-factor mean reverting model
    //
    //    dY(t) = -a Y(t) dt  +  sigma(t) dW_Y(t), Y(0)=0
    //    dZ(t) = -b Z(t) dt  +    eta(t) dW_Z(t), Z(0)=0
    //    dW_Y(t) dW_Z(t) = rho dt
    //
    template <class DateType, class PassiveType, class ActiveType>
    class TwoFactorMeanReversionModelT : public StochasticProcessT<DateType,PassiveType,ActiveType> {
        // from base class
        using typename StochasticProcessT<DateType,PassiveType,ActiveType>::MatA;


    protected:
    
        // container class definitions
        typedef std::vector<DateType>                      VecD;
        typedef std::vector<PassiveType>                   VecP; 
        typedef std::vector<ActiveType>                    VecA;
    
        // term structure for deterministic part
        Handle<IndexTermStructure> futureTS_;            // deterministic part
        
        // unique grid for time-dependent parameters
        VecD                       times_;   // time-grid of left-constant model parameter values
        // time-dependent parameters, left-piecewise constant on times_-grid
        VecA                       sigma_;  // volatility for Y
        VecA                       eta_;    // volatility for Z
        // time-homogeneous parameters
        PassiveType                a_;      // mean reversion for Y 		
        PassiveType                b_;      // mean reversion for Z
        PassiveType                rho_;    // correlation Y vs Z

        // functor for (co-)variance integration
        struct CovarianceHelper {
            ActiveType a_, b_;
            DateType T_;
            CovarianceHelper(const ActiveType a, const ActiveType b, const DateType T) : a_(a), b_(b), T_(T) {}
            // f(t) = exp{-(a+b)(T-t)}
            // F(t) = exp{-(a+b)(T-t)} / (a+b)
            inline PassiveType operator() (const DateType t) { return exp(-(a_+b_)*(T_-t)) / (a_+b_); }
        };

    public:
        // constructor
        
        TwoFactorMeanReversionModelT( const Handle<IndexTermStructure>&    futureTS,
                                      const VecD&                          times,
                                      const VecA&                          sigma,
                                      const VecA&                          eta,
                                      const PassiveType                     a,
                                      const PassiveType                     b,
                                      const PassiveType                     rho 
                                      )
        : futureTS_(futureTS), times_(times), sigma_(sigma), eta_(eta), a_(a), b_(b),	rho_(rho) {
            // check for valid parameter inputs
        }

        virtual ~TwoFactorMeanReversionModelT() = default;

        // inspectors
        inline const ActiveType sigma(const DateType t) const { return sigma_[TemplateAuxilliaries::idx(times_,t)]; }
        inline const ActiveType eta  (const DateType t) const { return eta_[TemplateAuxilliaries::idx(times_,t)];   }

        inline const PassiveType a()   const { return a_;   }
        inline const PassiveType b()   const { return b_;   }
        inline const PassiveType rho() const { return rho_; }
        
        // analytic formulas

        // deterministic part dependending on future index; overloaded in lognormal model
        inline virtual const ActiveType phi(const DateType t) const { return futureTS_->value(t); }

        // (future) variance of Y process
        inline ActiveType varianceY(const DateType t, const DateType T) const {
            CovarianceHelper F(a_,a_,T);
            VecA sigma2(sigma_.size());
            for (size_t k=0; k<sigma_.size(); ++k) sigma2[k] = sigma_[k]*sigma_[k];
            TemplateAuxilliaries::PieceWiseConstantIntegral<PassiveType,ActiveType,CovarianceHelper> integral(times_, sigma2, F);
            ActiveType variance = integral(t,T);
            return variance;
        }

        // (future) variance of Z process
        inline ActiveType varianceZ(const DateType t, const DateType T) const {
            CovarianceHelper F(b_,b_,T);
            VecA eta2(eta_.size());
            for (size_t k=0; k<sigma_.size(); ++k) eta2[k] = eta_[k]*eta_[k];
            TemplateAuxilliaries::PieceWiseConstantIntegral<PassiveType,ActiveType,CovarianceHelper> integral(times_, eta2, F);
            ActiveType variance = integral(t,T);
            return variance;
        }

        // (future) covariance X-Z process
        inline ActiveType covarianceYZ(const DateType t, const DateType T) const {
            CovarianceHelper F(a_,b_,T);
            VecA sigmaTimesEta(sigma_.size());
            for (size_t k=0; k<sigmaTimesEta.size(); ++k) sigmaTimesEta[k] = sigma_[k]*eta_[k];
            TemplateAuxilliaries::PieceWiseConstantIntegral<PassiveType,ActiveType,CovarianceHelper> integral(times_, sigmaTimesEta, F);
            ActiveType covariance = integral(t,T);
            return covariance;
        }

        
        // future expectation, implementation in derived classes
        inline virtual ActiveType futureAsset(const DateType t, const DateType T, const ActiveType Y, const ActiveType Z) {
            QL_REQUIRE(false,"futureAsset() should be overloaded in derived classes");
            return 0.0;
        }
        
        // variance of asset future, implementation in derived classes
        inline virtual ActiveType varianceFuture( const DateType startTime, const DateType expiryTime, const DateType settlementTime ) {
            QL_REQUIRE(false,"varianceFuture() should be overloaded in derived classes");
            return 0.0;
        }

        // basic instruments

        inline virtual ActiveType averageFuture ( const VecD& settlementTimes, const VecP& settlementWeights) {
            QL_REQUIRE(false,"averageFuture() should be overloaded in derived classes");
            return 0.0;
        }

        inline virtual ActiveType varianceAverageFuture ( const DateType expiryTime, const VecD& settlementTimes, const VecP& settlementWeights) {
            QL_REQUIRE(false,"varianceAverageFuture() should be overloaded in derived classes");
            return 0.0;
        }

        inline virtual ActiveType vanillaOption ( const DateType expiryTime, const VecD& settlementTimes, const VecP& settlementWeights, PassiveType strike, int callOrPut) {
            QL_REQUIRE(false,"vanillaOption() should be overloaded in derived classes");
            return 0.0;
        }

        // stochastic process interface

        // dimension of X
        inline virtual size_t size() { return 2; }
        // stochastic factors of x and z (maybe distinguish if trivially eta=0)
        inline virtual size_t factors() { return 2; }
        // initial values for simulation
        inline virtual VecP initialValues() {
            VecP X(2,0.0);
            return X;
        }

        // a[t,X(t)]
        inline virtual VecA drift( const DateType t, const VecA& X) {
            VecA a(2,0.0);
            // Y-variable
            a[0] = -a_ * X[0];
            // Z-variable
            a[1] = -b_ * X[1];
            return a;
        }
        // b[t,X(t)]
        inline virtual MatA diffusion( const DateType t, const VecA& X) {
            MatA B(2);
            B[0].resize(2,0.0);
            B[1].resize(2,0.0);
            // Y-variable sigma(t) dW_Y(t)
            B[0][0] = sigma(t); 
            B[0][1] = 0.0;
            // Z-variable eta(t) dW_Z(t)
            B[1][0] = eta(t) * rho_             ;
            B[1][1] = eta(t) * sqrt(1-rho_*rho_);
            // finished
            return B;
        }


        // integrate X1 = mu + nu dW
        inline virtual void evolve( const DateType t0, const VecA& X0, const DateType dt, const VecD& dW, VecA& X1 ) {
            // ensure X1 has size of X0
            // Y1 = exp(-a dt) Y0 + sqrt{VarY} dW_Y
            X1[0] = exp(-a_*dt)*X0[0] + sqrt(varianceY(t0,t0+dt))*dW[0];
            // Z1 = exp(-b dt) Z0 + sqrt{VarZ} dW_Z
            X1[1] = exp(-b_*dt)*X0[1] + sqrt(varianceZ(t0,t0+dt))*(rho_*dW[0] + sqrt(1-rho_*rho_)*dW[1]);
            return;
        }

        // stochastic process variables and payoffs

        inline virtual ActiveType asset(const DateType t, const VecA& X, const std::string& alias)                         { return futureAsset(t,t,X[0],X[1]);  }
        inline virtual ActiveType futureAsset(const DateType t, const DateType T, const VecA& X, const std::string& alias) { return futureAsset(t,T,X[0],X[1]);  }

    };

}

#undef _MIN_
#undef _MAX_

#endif  /* quantlib_template2fmeanreversionmodel_hpp */
