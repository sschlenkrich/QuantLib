/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2015 Sebastian Schlenkrich

*/



#ifndef quantlib_templatequadraticlvsvmodel_hpp
#define quantlib_templatequadraticlvsvmodel_hpp

#include <ql/shared_ptr.hpp>
#include <ql/errors.hpp>
#include <ql/experimental/templatemodels/stochasticprocessT.hpp>


namespace QuantLib {

    // Quadratic local-stochastic volatility model
    //
    //    dS(t) = phi(S) exp(z(t)) dW(t)
    //    dz(t) = -theta z(t) dt + nu dZ(t)
    //     z(0) = 0
    //    dW(t) dZ(t) = rho dt
    //    phi(S) = 0.5*curv*(S-S0)^2 + skew*(S-S0) + sigma0
    //
    template <class DateType, class PassiveType, class ActiveType>
    class QuadraticLVSVModelT : public StochasticProcessT<DateType,PassiveType,ActiveType> {
    private:
        ActiveType S0_, curv_, skew_, sigma0_, theta_, nu_, rho_;
    public:
        // from base class
        using typename StochasticProcessT<DateType, PassiveType, ActiveType>::VecP;
        using typename StochasticProcessT<DateType, PassiveType, ActiveType>::VecA;
        using typename StochasticProcessT<DateType, PassiveType, ActiveType>::MatA;


        // constructor
        QuadraticLVSVModelT( ActiveType S0, 
                             ActiveType curv, 
                             ActiveType skew, 
                             ActiveType sigma0,
                             ActiveType theta,
                             ActiveType nu,
                             ActiveType rho )
        : S0_(S0), curv_(curv), skew_(skew), sigma0_(sigma0), theta_(theta), nu_(nu), rho_(rho) {
            // check for valid parameter inputs
        }
        // stochastic process interface
        // dimension of X
        inline virtual size_t size() { return 2; }
        // stochastic factors of x and z (maybe distinguish if trivially eta=0)
        inline virtual size_t factors() { return 2; }
        // initial values for simulation
        inline virtual VecP initialValues() {
            VecP X(2);
            X[0] = S0_;
            X[1] = 0.0;
            return X;
        }
        // a[t,X(t)]
        inline virtual VecA drift( const DateType t, const VecA& X) {
            VecA a(2);
            // S-variable drift-less
            a[0] = 0.0;
            // z-variable -theta z(t)
            a[1] = -theta_* X[1];
            return a;
        }
        // b[t,X(t)]
        inline virtual MatA diffusion( const DateType t, const VecA& X) {
            MatA B(2);
            B[0].resize(2);
            B[1].resize(2);
            // S-variable phi(S) exp(z(t)) dW(t)
            ActiveType phi = 0.5*curv_*(X[0]-S0_)*(X[0]-S0_) + skew_*(X[0]-S0_) + sigma0_;
            if (phi < 0.0)          phi = 0.0;
            if (phi > 10.0*sigma0_) phi = 10.0*sigma0_;            
            B[0][0] = phi * exp(X[1]);
            B[0][1] = 0.0;
            // z-variable nu dZ(t)
            B[1][0] = rho_ * nu_;
            B[1][1] = sqrt(1-rho_*rho_) * nu_;
            // finished
            return B;
        }

        // the numeraire in the domestic currency used for discounting future payoffs
        inline virtual ActiveType numeraire(const DateType t, const VecA& X) { return 1.0; }

        // a domestic currency zero coupon bond
        inline virtual ActiveType zeroBond(const DateType t, const DateType T, const VecA& X) { return 1.0; }

        // default implementation for single-asset models
        inline virtual ActiveType asset(const DateType t, const VecA& X, const std::string& alias) { return X[0]; }

        // the expectation E^T in the domestic currency terminal meassure
        // this is required to calculate asset adjusters without knowing the implementation of the model
        inline virtual ActiveType forwardAsset(const DateType t, const DateType T, const VecA& X, const std::string& alias) { return asset(t, X, alias); }
        
    };

}

#undef _MIN_
#undef _MAX_

#endif  /* ifndef quantlib_templatequadraticlvsvmodel_hpp */
