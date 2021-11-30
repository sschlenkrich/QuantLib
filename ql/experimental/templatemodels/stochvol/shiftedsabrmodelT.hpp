/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2015 Sebastian Schlenkrich

*/



#ifndef quantlib_templateshiftedsabrmodel_hpp
#define quantlib_templateshiftedsabrmodel_hpp

#include <ql/shared_ptr.hpp>
#include <ql/errors.hpp>
#include <ql/experimental/templatemodels/stochasticprocessT.hpp>



#define _MIN_( a, b ) ( (a) < (b) ? (a) : (b) )
#define _MAX_( a, b ) ( (a) > (b) ? (a) : (b) )

namespace QuantLib {

    // SABR model with shift
    //
    //    dS(t) = z(t) (S(t) + lambda)^beta dW(t)
    //    dz(t) = nu z(t) dZ(t)
    //     z(0) = alpha
    //    dW(t) dZ(t) = rho dt
    //
    template <class DateType, class PassiveType, class ActiveType>
    class ShiftedSABRModelT : public StochasticProcessT<DateType,PassiveType,ActiveType> {
    private:
        ActiveType S0_, lambda_, alpha_, beta_, rho_, nu_;
    public:
        // from base class
        using typename StochasticProcessT<DateType, PassiveType, ActiveType>::VecD;
        using typename StochasticProcessT<DateType, PassiveType, ActiveType>::VecP;
        using typename StochasticProcessT<DateType, PassiveType, ActiveType>::VecA;
        using typename StochasticProcessT<DateType, PassiveType, ActiveType>::MatA;


        // constructor
        ShiftedSABRModelT( ActiveType S0, 
                           ActiveType lambda, 
                           ActiveType alpha, 
                           ActiveType beta,
                           ActiveType rho,
                           ActiveType nu )
        : S0_(S0), lambda_(lambda), alpha_(alpha), beta_(beta), rho_(rho), nu_(rho) {
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
            X[1] = alpha_;
            return X;
        }
        // a[t,X(t)]
        inline virtual VecA drift( const DateType t, const VecA& X) {
            VecA a(2);
            // S-variable drift-less
            a[0] = 0.0;
            // z-variable drift-less
            a[1] = 0.0;
            return a;
        }
        // b[t,X(t)]
        inline virtual MatA diffusion( const DateType t, const VecA& X) {
            MatA B(2);
            B[0].resize(2);
            B[1].resize(2);
            // S-variable z(t) (S(t) + lambda)^beta dW(t)
            B[0][0] = X[1] * pow(X[0] + lambda_, beta_);
            B[0][1] = 0.0;
            // z-variable nu z(t) dZ(t)
            B[1][0] = rho_ * nu_*X[1];
            B[1][1] = sqrt(1-rho_*rho_) * nu_*X[1];
            // finished
            return B;
        }

        inline virtual bool truncate( const DateType t, VecA& X ) { 
            // better check dimensions
            // S >= -lambda
            bool truncated = false;
            if ((beta_>0.0) & (X[0]<-lambda_)) {
                X[0] = -lambda_;
                truncated = true;
            }
            // z >= 0
            if (X[1]<0.0) {
                X[1] = 0.0;
                truncated = true;
            }
            return truncated;
        }

        // integrate X = ln(S + lambda) and v = ln(z) unless beta = 0
        inline virtual void evolve( const DateType t0, const VecA& X0, const DateType dt, const VecD& dW, VecA& X1 ) {
            // ensure X1 has size of X0
            // S-variable, need to distinguish normal, log-normal and in between
            if (beta_==0.0) {
                X1[0] = X0[0] + X0[1]*dW[0]*sqrt(dt);
            }
            if (beta_==1.0) {
                ActiveType vol = X0[1];
                X1[0] = (X0[0] + lambda_) * exp(-vol*vol/2.0*dt + vol*dW[0]*sqrt(dt)) - lambda_;
            }
            if ((beta_>0.0) & (beta_<1.0)) {
                X1[0] = -lambda_;  // fall back is absorption
                // check for absorbtion at -lambda
                if (X0[0]>-lambda_) {
                    ActiveType eps  = QL_EPSILON;
                    ActiveType SpLmin = pow( -X0[1]*X0[1]*dt/log(eps)/2.0, 1.0/2.0/(1.0-beta_) );
                    // only evolve if there is no absorbtion (yet)
                    if ((X0[0] + lambda_)>SpLmin) {
                        ActiveType vol = X0[1] * pow(X0[0] + lambda_, beta_-1.0);
                        X1[0] = (X0[0] + lambda_) * exp(-vol*vol/2.0*dt + vol*dW[0]*sqrt(dt)) - lambda_;
                    }
                }
            }
            // z-variable
            // X1[1] = X0[1] + (b[1][0]*dW[0]+b[1][1]*dW[1])*sqrt(dt);
            ActiveType dZ = rho_*dW[0] + sqrt(1-rho_*rho_)*dW[1];
            X1[1] = X0[1] * exp( -nu_*nu_/2.0*dt + nu_*dZ*sqrt(dt) );
            // truncate, should not be neccessary
            // truncate(t0+dt, X1);
            return;
        }

    };

}

#undef _MIN_
#undef _MAX_

#endif  /* ifndef quantlib_templateshiftedsabrmodel_hpp */
