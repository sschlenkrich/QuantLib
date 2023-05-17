/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2018 Sebastian Schlenkrich

*/



#ifndef quantlib_templateg2pp_hpp
#define quantlib_templateg2pp_hpp

#include <ql/shared_ptr.hpp>
#include <ql/errors.hpp>
#include <ql/termstructures/yieldtermstructure.hpp>
#include <ql/experimental/templatemodels/stochasticprocessT.hpp>

#define _MIN_( a, b ) ( (a) < (b) ? (a) : (b) )
#define _MAX_( a, b ) ( (a) > (b) ? (a) : (b) )

namespace QuantLib {

    // G2++ model
    //     r(t) = phi(t) + x(t) + y(t)
    //    dx(t) = -a x(t) dt + sigma dW_1(t)
    //    dy(t) = -b y(t) dt + eta dW_2(t)
    //    dB(t) = B(t) r(t) dt
    //    dW_1(t) dW_2(t) = rho dt
    //
    template <class DateType, class PassiveType, class ActiveType>
    class G2ppT : public StochasticProcessT<DateType,PassiveType,ActiveType> {
        // from base class
        using typename StochasticProcessT<DateType, PassiveType, ActiveType>::VecD;
        using typename StochasticProcessT<DateType, PassiveType, ActiveType>::VecP;
        using typename StochasticProcessT<DateType, PassiveType, ActiveType>::VecA;
        using typename StochasticProcessT<DateType, PassiveType, ActiveType>::MatA;


    private:
        ActiveType sigma_, eta_, a_, b_, rho_;
        Handle<YieldTermStructure> termStructure_;  // the yield curve is   

        inline ActiveType V(DateType t) const {
            ActiveType expat = exp(-a_*t);
            ActiveType expbt = exp(-b_*t);
            ActiveType cx = sigma_ / a_;
            ActiveType cy = eta_ / b_;
            ActiveType valuex = cx*cx*(t + (2.0*expat - 0.5*expat*expat - 1.5) / a_);
            ActiveType valuey = cy*cy*(t + (2.0*expbt - 0.5*expbt*expbt - 1.5) / b_);
            ActiveType value = 2.0*rho_*cx*cy* (t + (expat - 1.0) / a_
                + (expbt - 1.0) / b_
                - (expat*expbt - 1.0) / (a_ + b_));
            return valuex + valuey + value;
        }

        inline ActiveType A(DateType t, DateType T) const {
            return termStructure_->discount(T) / termStructure_->discount(t) * exp(0.5*(V(T - t) - V(T) + V(t)));
        }

        inline ActiveType B(PassiveType x, DateType t) const {
            return (1.0 - exp(-x*t)) / x;
        }

    public:
        // constructor
        G2ppT( const Handle<YieldTermStructure>& termStructure,
               const ActiveType sigma, 
               const ActiveType eta, 
               const ActiveType a, 
               const ActiveType b, 
               const ActiveType rho )
        : termStructure_(termStructure), sigma_(sigma), eta_(eta), a_(a), b_(b), rho_(rho) {
            // check for valid parameter inputs

        }         
        // stochastic process interface
        // dimension of X = [ x, y, B ]
        inline virtual size_t size() { return 3; } 
        // stochastic factors (x and y)
        inline virtual size_t factors() { return 2; }
        // initial values for simulation
        inline virtual VecP initialValues() {
            VecP X(3);
            X[0] = 0.0;  // x(t)
            X[1] = 0.0;  // y(t)
            X[2] = 1.0;  // B(t) bank account numeraire
            return X;
        }
        // a[t,X(t)]
        inline virtual VecA drift( const DateType t, const VecA& X) {
            QL_FAIL("Drift is not implemented");
            return VecA(0);
        }
        // b[t,X(t)]
        inline virtual MatA diffusion( const DateType t, const VecA& X) {
            QL_FAIL("Diffusion is not implemented");
            return MatA(0);
        }

        // integrate X1 = X0 + drift()*dt + diffusion()*dW*sqrt(dt)
        inline void evolve(const DateType t0, const VecA& X0, const DateType dt, const VecD& dW, VecA& X1) {
            // ensure X1 has size of X0
            ActiveType ExpX  = X0[0] * exp(-a_*dt);
            ActiveType ExpY  = X0[1] * exp(-b_*dt);
            ActiveType VarX  = sigma_*sigma_ / 2.0 / a_*(1.0 - exp(-2.0*a_*dt));
            ActiveType VarY  = eta_  *eta_   / 2.0 / b_*(1.0 - exp(-2.0*b_*dt));
            ActiveType CovXY = rho_*sigma_*eta_ / (a_ + b_)*(1.0 - exp(-(a_ + b_)*dt));
            ActiveType corr = CovXY / sqrt(VarX) / sqrt(VarY);
            ActiveType dZ = corr*dW[0] + sqrt(1 - corr*corr)*dW[1];
            X1[0] = ExpX + sqrt(VarX)*dW[0];
            X1[1] = ExpY + sqrt(VarY)*dZ;
            // Brigo, Corollary 4.2.1
            ActiveType expIntPhi = termStructure_->discount(t0) / termStructure_->discount(t0 + dt) * exp(0.5*(V(t0 + dt) - V(t0)));
            ActiveType expIntX = exp(0.5*(X0[0] + X1[0])*dt);  // numerical integration
            ActiveType expIntY = exp(0.5*(X0[1] + X1[1])*dt);  // numerical integration
            X1[2] = X0[2] * expIntPhi * expIntX * expIntY;  // using that r = phi + x + y
        }

        inline ActiveType zeroBond(const DateType t, const DateType T, const VecA& X) {
            QL_REQUIRE(t <= T, "G2++ Model ZeroBond t <= T required");
            if (t == T) return (ActiveType)1.0;
            return A(t, T) * exp(-B(a_, (T - t))*X[0] - B(b_, (T - t))*X[1]);  // assume X = [x, y, B]
        }

        inline ActiveType numeraire(const DateType t, const VecA& X) {
            return X[2];
        }

    };

}

#undef _MIN_
#undef _MAX_

#endif  /* ifndef quantlib_templateshiftedsabrmodel_hpp */
