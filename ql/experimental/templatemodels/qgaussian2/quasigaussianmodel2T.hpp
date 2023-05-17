/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2018 Sebastian Schlenkrich

*/

/*! \file quasigaussianmodel2T.hpp
    \brief (MC) pricing for multi-factor quasi-Gaussian model with stochastic vol
           Re-factoring initial implementation with alternative local volatility specification and clean-up
               
               r(t) = f(0,t) + 1^T*x(t)
               
               dx(t)     = [ y(t)*1 - a*x(t) ] dt                                        + sqrt[z(t)]*sigma_x^T(t,x,y) dW
               dy(t)     = [ z(t)*sigma_x^T(t,x,y)*sigma_x(t,x,y) - a*y(t) - y(t)*a ] dt
               dz(t)     = theta [ z0 - z(t) ] dt                                        + eta(t)*sqrt[z(t)]           dZ
               ds(t)     = r(t) dt  ( s(t) = int_0^t r(s) ds, for bank account numeraire)
               
           All methods are template based to allow incorporation of Automatic Differentiation
           tools
*/


#ifndef quantlib_templatequasigaussian2_hpp
#define quantlib_templatequasigaussian2_hpp

#include <ql/termstructures/yieldtermstructure.hpp>
#include <ql/experimental/templatemodels/stochasticprocessT.hpp>
#include <ql/experimental/templatemodels/auxilliaries/qrfactorisationT.hpp>
#include <ql/experimental/templatemodels/auxilliaries/choleskyfactorisationT.hpp>


namespace QuantLib {

    // Declaration of the quasi-Gaussian model class
    template <class DateType, class PassiveType, class ActiveType>
    class QuasiGaussianModel2T : public StochasticProcessT<DateType, PassiveType, ActiveType> {
    protected:

        // container class definitions
        typedef std::vector<DateType>                      VecD;
        typedef std::vector<PassiveType>                   VecP; 
        typedef std::vector<ActiveType>                    VecA;
        typedef std::vector< std::vector<DateType> >       MatD;
        typedef std::vector< std::vector<PassiveType> >    MatP;
        typedef std::vector< std::vector<ActiveType> >     MatA;

        // attributes defining the model
        Handle<YieldTermStructure> termStructure_;  // the yield curve is assumed to be passive
        // number of yield curve factors (excluding stoch. vol)
        size_t                     d_;       // (d+1)-dimensional Brownian motion for [x(t), z(t)]^T
        // unique grid for time-dependent parameters
        VecD                       times_;   // time-grid of left-constant model parameter values
        // time-dependent parameters, left-piecewise constant on times_-grid
        // we implement a local volatility function sigma_f = curve*|f-f0| + smile*(f-f0) + sigma
        MatA                       sigma_;  // volatility
        MatA                       slope_;  // skew
        MatA                       curve_;  // smile
        VecA                       eta_;    // vol-of-vol
        // time-homogeneous parameters
        VecP                       delta_;   // maturity of benchmark rates f(t,t+delta_i) 		
        VecP                       chi_;     // mean reversions
        MatP                       Gamma_;   // (benchmark rate) correlation matrix
        // stochastic volatility process parameters
        PassiveType                theta_;   // mean reversion speed
        PassiveType                z0_;      // mean reversion level z0=z(0)=1

        // additional parameters (calculated at initialisation)
        VecP                       f0_;      // f(0,delta_i)
        MatP                       DfT_;     // factorized correlation matrix Df^T with Df^T * Df = Gamma
        MatP                       HHfInv_;  // weighting matrix H*Hf^-1 = [ exp{-chi_j*delta_i} ]^-1
        PassiveType                condHHf_; // condition number of H*Hf; this is calculated from singular values and acts as control to avoid numerical instability

        inline virtual bool checkModelParameters(const bool throwException=true){
            bool ok = true;
            // check yield curve...
            // non-zero dimension
            if (d_<1)               { ok = false; if (throwException) QL_REQUIRE(false,"QuasiGaussianModel number of factors larger zero required."); }
            // non-empty time-grid
            size_t n=times_.size();
            if (n<1)                { ok = false; if (throwException) QL_REQUIRE(false,"QuasiGaussianModel non-empty time-grid required."); }
            // dimensions of time-dependent parameters
            if (sigma_.size()!=d_) { ok = false; if (throwException) QL_REQUIRE(false,"QuasiGaussianModel wrong sigma dimension."); }
            if (slope_.size()!=d_) { ok = false; if (throwException) QL_REQUIRE(false,"QuasiGaussianModel wrong slope dimension.");  }
            if (curve_.size()!=d_) { ok = false; if (throwException) QL_REQUIRE(false,"QuasiGaussianModel wrong curve dimension."); }
            for (size_t k=0; k<d_; ++k) {
                if (sigma_[k].size()!=n) { ok = false; if (throwException) QL_REQUIRE(false,"QuasiGaussianModel wrong sigma time dimension."); }
                if (slope_[k].size()!=n) { ok = false; if (throwException) QL_REQUIRE(false,"QuasiGaussianModel wrong slope time dimension."); }
                if (curve_[k].size()!=n) { ok = false; if (throwException) QL_REQUIRE(false,"QuasiGaussianModel wrong curve time dimension."); }
            }
            if (eta_.size() !=n) { ok = false; if (throwException) QL_REQUIRE(false,"QuasiGaussianModel wrong eta time dimension.");    }
            // dimensions of time-homogeneous parameters
            if (delta_.size()!=d_)   { ok = false; if (throwException) QL_REQUIRE(false,"QuasiGaussianModel wrong delta dimension."); }
            if (chi_.size()  !=d_)   { ok = false; if (throwException) QL_REQUIRE(false,"QuasiGaussianModel wrong chi dimension.");   }
            if (Gamma_.size()  !=d_) { ok = false; if (throwException) QL_REQUIRE(false,"QuasiGaussianModel wrong Gamma i-dimension."); }
            for (size_t k=0; k<d_; ++k) {
                if (Gamma_[k].size() !=d_) { ok = false; if (throwException) QL_REQUIRE(false,"QuasiGaussianModel wrong Gamma j-dimension."); }
            }
            // plausible parameter values
            // ascending time-grid
            for (size_t k=0; k<times_.size()-1; ++k) {
                if (times_[k]>=times_[k+1]) { ok = false; if (throwException) QL_REQUIRE(false,"QuasiGaussianModel ascending time-grid required."); }
            }
            // non-negative values
            for (size_t j=0; j<n; ++j) {
                for (size_t i=0; i<d_; ++i) {
                    if (sigma_[i][j]<0.0) { ok = false; if (throwException) QL_REQUIRE(false,"QuasiGaussianModel sigma>=0 required."); }
                }
                if (eta_[j]<0.0) { ok = false; if (throwException) QL_REQUIRE(false,"QuasiGaussianModel eta>=0 required.");      }
            }
            // positive/ascending values
            if (delta_[0]<=0.0) { ok = false; if (throwException) QL_REQUIRE(false,"QuasiGaussianModel delta>0 required."); }
            if (chi_[0]  <=0.0) { ok = false; if (throwException) QL_REQUIRE(false,"QuasiGaussianModel chi>0 required."); }  // we might allow negative mean reversion...
            for (size_t k=0; k<d_-1; ++k) {
                if (delta_[k]>=delta_[k+1]) { ok = false; if (throwException) QL_REQUIRE(false,"QuasiGaussianModel ascending delta values required."); }
                if (chi_[k]  >=chi_[k+1])   { ok = false; if (throwException) QL_REQUIRE(false,"QuasiGaussianModel ascending chi values required."); }
            }
            // plausible correlation values
            for (size_t i=0; i<d_; ++i) {
                for (size_t j=i; j<d_; ++j) {
                    //if (Gamma_[i][j]<0.0) {
                    //	ok = false;
                    //	if (throwException) QL_REQUIRE(false,"QuasiGaussianModel Gamma[i][j]>=0 required."); 
                    //}
                    if (i==j) {
                        if (Gamma_[i][j]!=1.0) {
                            ok = false;
                            if (throwException) QL_REQUIRE(false,"QuasiGaussianModel Gamma[i][i]=1 required."); 
                        }
                    }
                    if (Gamma_[i][j]!=Gamma_[j][i]) {
                        ok = false;
                        if (throwException) QL_REQUIRE(false,"QuasiGaussianModel Gamma[i][j]=Gamma[j][i] required."); 
                    }
                    //if (i<j) {
                    //	if((Gamma_[i][j-1]<=Gamma_[i][j])|(Gamma_[i+1][j]<=Gamma_[i][j])) {
                    //		ok = false;
                    //		if (throwException) QL_REQUIRE(false,"QuasiGaussianModel Gamma descending sub-diagonals required."); 
                    //	}
                    //}
                }
            }
            // stochastic vol parameters
            if (theta_<=0.0) { ok = false; if (throwException) QL_REQUIRE(false,"QuasiGaussianModel theta>0 required."); }
            if (z0_!=1.0)    { ok = false; if (throwException) QL_REQUIRE(false,"QuasiGaussianModel z0=1 required.");    }
            // finished
            return ok;
        }

        // evaluate Df^T with Df^T * Df = Gamma and H*Hf^-1 via singular value decomposition
        // return false (and throw exception) on error
        inline bool factorMatrices(const bool throwException=true){
            bool ok = true;
            DfT_ = TemplateAuxilliaries::cholesky(Gamma_);
            // [ Hf H^{-1} ] = [ exp{-chi_j*delta_i} ]
            MatA A(d_, VecA(d_, 0.0));
            for (size_t i = 0; i<d_; ++i) {
                for (size_t j = 0; j<d_; ++j) {
                    A[i][j] = exp(-chi_[j] * delta_[i]);
                }
            }
            HHfInv_ = TemplateAuxilliaries::qrinverse(A);
            condHHf_ = 0.0;  // we don't have this information with QR factorisation
            return ok;
        }

    protected:
        // analytic formulas based on internal state representation
        // these functions need to be consistent to corresponding public methods below!

        // interface for accessing the state of the yield curve
        class State {
            const VecA& X_;
            const size_t d_;
        public:
            State(const VecA& X, const size_t d) : X_(X), d_(d) {
                QL_REQUIRE(X.size() == d + d*d + 1 + 1, "TemplateQuasiGaussianModel::State Constructor: Dimensions mismatch.");
            }
            inline const ActiveType x(const size_t i) const { return X_[i]; }
            inline const ActiveType y(const size_t i, const size_t j) const {
                return X_[d_ + i*d_ + j];  // y row-wise
            }
            inline const ActiveType z() const { return X_[d_ + d_*d_]; }
            inline const ActiveType s() const { return X_[d_ + d_*d_ + 1]; }
        };

        inline ActiveType shortRate(const DateType t, const State& s) {
            ActiveType r = termStructure_->forwardRate(t, t, Continuous);
            for (size_t k = 0; k<d_; ++k) r += s.x(k);
            return r;
        }

        inline ActiveType forwardRate(const DateType t, const DateType T, const State& s) {
            ActiveType f = termStructure_->forwardRate(T, T, Continuous);  // check t,T
            for (size_t i = 0; i<d_; ++i) {
                ActiveType tmp = s.x(i);
                for (size_t j = 0; j<d_; ++j) tmp += s.y(i,j) * G(j, t, T);
                f += exp(-chi_[i] * (T - t)) * tmp;
            }
            return f;
        }

        // f - f0 for local volatility calculation
        inline ActiveType deltaForwardRate(const DateType t, const DateType T, const State& s) {
            ActiveType df = 0.0;
            for (size_t i = 0; i<d_; ++i) {
                ActiveType tmp = s.x(i);
                for (size_t j = 0; j<d_; ++j) tmp += s.y(i, j) * G(j, t, T);
                df += exp(-chi_[i] * (T - t)) * tmp;
            }
            return df;
        }

        inline ActiveType ZeroBond(const DateType t, const DateType T, const State& s) {
            QL_REQUIRE(t <= T, "QuasiGaussianModel ZeroBond t <= T required");
            if (t == T) return (ActiveType)1.0;
            PassiveType DF1 = termStructure_->discount(t);
            PassiveType DF2 = termStructure_->discount(T);
            ActiveType  Gx = 0;   // G^T * x
            for (size_t i = 0; i<d_; ++i) Gx += s.x(i) * G(i, t, T);
            ActiveType  GyG = 0;   // G^T * y * G
            for (size_t i = 0; i<d_; ++i) {
                ActiveType tmp = 0;
                for (size_t j = 0; j<d_; ++j) tmp += s.y(i,j) * G(j, t, T);
                GyG += G(i, t, T)*tmp;
            }
            ActiveType ZCB = DF2 / DF1 * exp(-Gx - 0.5*GyG);
            return ZCB;
        }

        inline VecA sigma_f(const DateType t, const State& s) { // diagonal vector
            VecA res(d_);
            PassiveType eps = 1.0e-4;  // 1bp lower cutoff for volatility
            for (size_t k = 0; k < d_; ++k) {
                ActiveType df = deltaForwardRate(t, t + delta_[k], s);
                //res[k] = sigma(k, t) + slope(k, t)*df + 0.5*curve(k, t)*df*df;
                res[k] = sigma(k, t) + slope(k, t)*df + curve(k, t)*((df<0.0)?(-df):(df));
                if (res[k] < eps) res[k] = eps;
                if (res[k] > 10.0*sigma(k, t)) res[k] = 10.0*sigma(k, t); // cap local volatility to avoid exploding states
            }
            return res;
        }

        inline virtual MatA sigma_xT(const DateType t, const State& s) {  // sigma_x^T
            MatA tmp(d_), res(d_);
            VecA sigmaf = sigma_f(t, s);
            // tmp = sigma_f * Df^T
            for (size_t i = 0; i<d_; ++i) {
                tmp[i].resize(d_);
                for (size_t j = 0; j<d_; ++j) {
                    tmp[i][j] = sigmaf[i] * DfT_[i][j];
                }
            }
            // res = H*Hf^-1 * tmp
            for (size_t i = 0; i<d_; ++i) {
                res[i].resize(d_);
                for (size_t j = 0; j<d_; ++j) {
                    res[i][j] = 0;
                    for (size_t k = 0; k<d_; ++k) res[i][j] += HHfInv_[i][k] * tmp[k][j];
                }
            }
            return res;
        }

    public:  

        QuasiGaussianModel2T(
            const Handle<YieldTermStructure>& termStructure,
            // number of yield curve factors (excluding stoch. vol)
            const size_t                d,       // (d+1)-dimensional Brownian motion for [x(t), z(t)]^T
            // unique grid for time-dependent parameters
            const VecD &                times,   // time-grid of left-constant model parameter values
            // time-dependent parameters, left-piecewise constant on times_-grid
            const MatA &                sigma,   // volatility
            const MatA &                slope,   // skew
            const MatA &                curve,   // smile
            const VecA &                eta,     // vol-of-vol
            // time-homogeneous parameters
            const VecP &                delta,   // maturity of benchmark rates f(t,t+delta_i) 		
            const VecP &                chi,     // mean reversions
            const MatP &                Gamma,   // (benchmark rate) correlation matrix
            const PassiveType           theta    // mean reversion speed
            ) : termStructure_(termStructure), d_(d), times_(times), sigma_(sigma), slope_(slope), curve_(curve), eta_(eta),
                delta_(delta), chi_(chi), Gamma_(Gamma), theta_(theta), z0_((PassiveType)1.0) {
                checkModelParameters();
                // calculate  DfT_
                // calculate  HHfInv_
                factorMatrices();
                f0_.resize(d_);
                for (size_t k = 0; k < d_; ++k) f0_[k] = termStructure_->forwardRate(delta_[k], delta_[k], Continuous);
            }

        virtual ~QuasiGaussianModel2T() = default;

        // helpers

        inline size_t maxidx( const size_t i ) { return (i<d_) ? i : d_-1; }

        // find index in ascending vector, evaluate n s.t. t[n-1] < t <= t[n]
        inline size_t idx( const DateType t ) { 
            if ((t <= times_[0]) | (times_.size()<2)) return 0;
            if (t >  times_[times_.size() - 2])       return times_.size() - 1;
            // bisection search
            size_t a = 0, b = times_.size() - 2;
            while (b-a>1) {
                size_t s = (a + b) / 2;
                if (t <= times_[s]) b = s;
                else                a = s;
            }
            return b;
        }

        // unsafe update
        inline void update(const MatA &                sigma,   // volatility
                           const MatA &                slope,   // skew
                           const VecA &                eta) {   // vol-of-vol
            sigma_ = sigma, slope_ = slope, eta_ = eta;
        }

        // unsafe update with curve
        inline void update(const MatA &                sigma,   // volatility
                           const MatA &                slope,   // skew
                           const MatA &                curve) { // curve
            sigma_ = sigma, slope_ = slope, curve_ = curve;
        }


        // clone the model
        virtual ext::shared_ptr<QuasiGaussianModel2T> clone() { return ext::shared_ptr<QuasiGaussianModel2T>(new QuasiGaussianModel2T(*this)); }

        // inspectors
        inline const Handle<YieldTermStructure> termStructure() { return termStructure_; }
        inline const VecD& times()   { return times_; }
        inline const MatA& sigma()   { return sigma_; }
        inline const MatA& slope()   { return slope_; }
        inline const MatA& curve()   { return curve_; }
        inline const VecA& eta()     { return  eta_;  }

        inline const MatP& DfT()    { return DfT_;    }
        inline const MatP& HHfInv() { return HHfInv_; }
        inline const VecP& delta()  { return delta_;  }
        inline const VecP& chi()    { return chi_;    }
        inline const ActiveType theta() { return theta_; } 
        inline const ActiveType z0()    { return z0_;    } 

        inline const PassiveType condHHf() { return condHHf_; }

        // parameter functions (no dimension checks)
        inline virtual ActiveType sigma ( const size_t i, const DateType t) { return sigma_[maxidx(i)][idx(t)]; }
        inline virtual ActiveType slope ( const size_t i, const DateType t) { return slope_[maxidx(i)][idx(t)]; }
        inline virtual ActiveType curve ( const size_t i, const DateType t) { return curve_[maxidx(i)][idx(t)]; }
        inline virtual ActiveType eta   ( const DateType t)                 { return eta_[idx(t)];              }

        // analytic formulas

        inline ActiveType G(const size_t i, const DateType t, const DateType T) { return (1.0-exp(-chi_[i]*(T-t)))/chi_[i]; }

        inline
        ActiveType shortRate ( const DateType t, const VecA& x ) {
            ActiveType r = termStructure_->forwardRate(t,t,Continuous);
            for (size_t k=0; k<d_; ++k) r += x[k];
            return r;
        }

        inline
        ActiveType forwardRate( const DateType t, const DateType T, const VecA& x, const MatA&  y) {
            ActiveType f = termStructure_->forwardRate(T,T,Continuous);  // check t,T
            for (size_t i=0; i<d_; ++i) {
                ActiveType tmp = x[i];
                for (size_t j=0; j<d_; ++j) tmp += y[i][j]*G(j,t,T);
                f += exp(-chi_[i]*(T-t)) * tmp;
            }
            return f;
        }

        inline  // f - f0 for local volatility calculation
        ActiveType deltaForwardRate( const DateType t, const DateType T, const VecA& x, const MatA&  y) {
            ActiveType df = 0.0;
            for (size_t i=0; i<d_; ++i) {
                ActiveType tmp = x[i];
                for (size_t j=0; j<d_; ++j) tmp += y[i][j]*G(j,t,T);
                df += exp(-chi_[i]*(T-t)) * tmp;
            }
            return df;
        }

        inline
        ActiveType ZeroBond( const DateType t, const DateType T, const VecA& x, const MatA&  y) {
            QL_REQUIRE(t<=T,"QuasiGaussianModel ZeroBond t <= T required");
            if (t==T) return (ActiveType)1.0;
            PassiveType DF1  = termStructure_->discount(t);
            PassiveType DF2  = termStructure_->discount(T);
            ActiveType  Gx   = 0;   // G^T * x
            for (size_t i=0; i<d_; ++i) Gx += x[i]*G(i,t,T);
            ActiveType  GyG  = 0;   // G^T * y * G
            for (size_t i=0; i<d_; ++i) {
                ActiveType tmp = 0;
                for (size_t j=0; j<d_; ++j) tmp += y[i][j]*G(j,t,T);
                GyG += G(i,t,T)*tmp;
            }
            ActiveType ZCB = DF2 / DF1 * exp(-Gx - 0.5*GyG);
            return ZCB;
        }

        inline  // diagonal vector
        VecA sigma_f( const DateType t, const VecA& x, const MatA&  y) {
            VecA res(d_);
            PassiveType eps = 1.0e-4;  // 1bp lower cutoff for volatility
            for (size_t k = 0; k < d_; ++k) {
                ActiveType df = deltaForwardRate(t, t + delta_[k], x, y);
                //res[k] = sigma(k, t) + slope(k, t)*df + 0.5*curve(k, t)*df*df;
                res[k] = sigma(k, t) + slope(k, t)*df + curve(k, t)*((df<0.0) ? (-df) : (df));
                if (res[k] < eps) res[k] = eps;
                if (res[k] > 10.0*sigma(k, t)) res[k] = 10.0*sigma(k, t); // cap local volatility to avoid exploding states
            }
            return res;
        }

        inline virtual // sigma_x^T
        MatA sigma_xT( const DateType t, const VecA& x, const MatA&  y) {
            MatA tmp(d_), res(d_);
            VecA sigmaf = sigma_f(t,x,y);
            // tmp = sigma_f * Df^T
            for (size_t i=0; i<d_; ++i) {
                tmp[i].resize(d_);
                for (size_t j=0; j<d_; ++j) {
                    tmp[i][j] = sigmaf[i] * DfT_[i][j];
                }
            }
            // res = H*Hf^-1 * tmp
            for (size_t i=0; i<d_; ++i) {
                res[i].resize(d_);
                for (size_t j=0; j<d_; ++j) {
                    res[i][j] = 0;
                    for (size_t k=0; k<d_; ++k) res[i][j] += HHfInv_[i][k]*tmp[k][j];
                }
            }
            return res;
        }

        // conditional moments of vol process used for z-integration, Piterbarg, 8.3.3.
        // E[ z(T) | z(t) ]
        inline virtual ActiveType expectationZ( DateType t, ActiveType zt, DateType dT ) {
            return z0_ + (zt - z0_)*exp(-theta_*dT);
        }
        // Var[ z(T) | z(t) ]
        inline virtual ActiveType varianceZ( DateType t, ActiveType zt, DateType dT ) {
            ActiveType expmThDT = exp(-theta_*dT);
            ActiveType onemETDT = 1 - expmThDT;
            ActiveType eta2oThe = eta(t+dT/2.0)*eta(t+dT/2.0)/theta_;  // approx eta(t)=eta for s \in [t, t+dT]
            return zt*eta2oThe*expmThDT*onemETDT + z0_*eta2oThe/2.0*onemETDT*onemETDT; 
        }

        // subset of QL's StochasticProcess interface for X = [ x, y, z, d ] (y row-wise)
        // with dX = a[t,X(t)] dt + b[t,X(t)] dW

        // dimension of X
        inline size_t size()    { return d_ + d_*d_ + 1 + 1; }
        // stochastic factors of x and z (maybe distinguish if trivially eta=0)
        inline size_t factors() { return d_ + 1; }
        // initial values for simulation
        inline VecP initialValues() {
            VecP X(size());
            for (size_t k=0; k<d_ + d_*d_; ++k) X[k] = 0.0;  // x(0), y(0)
            X[d_+d_*d_]                              = 1.0;  // z(0)
            X[d_+d_*d_+1]                            = 0.0;  // s(0)
            return X;
        }

        // a[t,X(t)]
        inline VecA drift( const DateType t, const VecA& X) {
            VecA a(size());
            State state(X,d_);
            // x-variable [ y(t)*1 - chi*x(t) ]
            for (size_t k=0; k<d_; ++k) {
                a[k] = -chi_[k]*state.x(k);
                for (size_t j=0; j<d_; ++j) a[k] += state.y(k,j);
            }
            // y-variable [ z(t)*sigma_x^T(t,x,y)*sigma_x(t,x,y) - chi*y(t) - y(t)*chi ]
            MatA sigmaxT = sigma_xT(t,state);
            for (size_t i=0; i<d_; ++i) {
                for (size_t j=0; j<d_; ++j) {
                    a[d_+i*d_+j] = 0.0;
                    for (size_t k=0; k<d_; ++k) a[d_+i*d_+j] += sigmaxT[i][k]*sigmaxT[k][j];
                    a[d_+i*d_+j] *= ((state.z()>0)?(state.z()):(0.0));  // full truncation
                    a[d_+i*d_+j] -= (chi_[i]+chi_[j])*state.y(i,j);
                }
            }
            // z-variable theta [ z0 - z(t)^+ ]  (full truncation)
            //a[d_+d_*d_] = theta_*(z0_ - ((state.z>0)?(state.z):(0.0)));
            // push to positive teritory
            a[d_+d_*d_] = theta_*(z0_ - state.z());
            // s-variable r(t)
            a[d_+d_*d_+1] = shortRate(t,state);
            // finished
            return a;
        }

        // b[t,X(t)]
        inline MatA diffusion( const DateType t, const VecA& X) {
            MatA b(size());
            for (size_t k=0; k<size(); ++k) b[k].resize(factors());
            State state(X,d_);
            ActiveType sqrtz = ((state.z()>0)?(sqrt(state.z())):(0.0));   // full truncation
            MatA sigmaxT = sigma_xT(t,state);
            // x-variable sqrt[z(t)]*sigma_x^T(t,x,y)
            for (size_t i=0; i<d_; ++i) {
                for (size_t j=0; j<d_; ++j) b[i][j] =  sqrtz * sigmaxT[i][j];
                b[i][d_] = 0.0;
            }
            // y-variable 0
            for (size_t i=d_; i<d_+d_*d_; ++i) {
                for (size_t j=0; j<d_+1; ++j) {
                    b[i][j] = 0.0;
                }
            }
            // z-variable eta(t)*sqrt[z(t)]
            for (size_t j=0; j<d_; ++j) b[d_+d_*d_][j] = 0.0;
            b[d_+d_*d_][d_] = eta(t)*sqrtz;
            // s-variable 0
            for (size_t j=0; j<d_+1; ++j) b[d_+d_*d_+1][j] = 0.0;
            // finished
            return b;
        }

        // simulate Quasi-Gaussian model as Gaussian model with frozen vol matrix
        inline void evolve(const DateType t0, const VecA& X0, const DateType dt, const VecD& dW, VecA& X1) {
            // first we simulate stochastic vol via lognormal approximation...
            ActiveType e = expectationZ(t0, X0[X0.size() - 2], dt);
            ActiveType v = varianceZ(t0, X0[X0.size() - 2], dt);
            ActiveType dZ = dW[dW.size() - 1];  // last risk factor is for vol process
            ActiveType si = sqrt(log(1.0 + v / e / e));
            ActiveType mu = log(e) - si*si / 2.0;
            X1[X1.size() - 2] = exp(mu + si*dZ);
            ActiveType averageZ = 0.5*(X0[X0.size() - 2] + X1[X1.size() - 2]);  // we freeze z for subsequent calculation
            // next we need V = z * sigmaxT sigmax
            State state(X0, d_);
            MatA sigmaxT = sigma_xT(t0+0.5*dt, state);  // assume constant vol on (t0, t0+dt)
            MatA V(d_, VecA(d_, 0.0));
            for (size_t i = 0; i < d_; ++i) {
                for (size_t j = 0; j <= i; ++j) {
                    for (size_t k = 0; k < d_; ++k) V[i][j] += sigmaxT[i][k] * sigmaxT[j][k];  // sigma_x^T * sigma_x (!)
                    V[i][j] *= averageZ;
                    if (j < i) V[j][i] = V[i][j];  // exploit symmetry
                }
            }
            // we also calculate intermediate variables; these do strictly depend on X and could be cached as well
            VecP expMinusChiDT(d_);
            for (size_t i = 0; i < d_; ++i) expMinusChiDT[i] = exp(-chi_[i] * dt);
            VecP OneMinusExpMinusChiDT(d_);
            for (size_t i = 0; i < d_; ++i) OneMinusExpMinusChiDT[i] = 1.0 - expMinusChiDT[i];
            // now we can calculate y
            MatA a(d_, VecA(d_, 0.0)), b(d_, VecA(d_,0.0));
            for (size_t i = 0; i < d_; ++i) {
                for (size_t j = 0; j <= i; ++j) {
                    b[i][j] = V[i][j] / (chi_[i] + chi_[j]);
                    a[i][j] = X0[d_ + i*d_ + j] - b[i][j];
                    if (j < i) {
                        b[j][i] = b[i][j];
                        a[j][i] = a[i][j];
                    }
                }
            }
            for (size_t i = 0; i < d_; ++i)
                for (size_t j = 0; j < d_; ++j) X1[d_ + i*d_ + j] = a[i][j] * expMinusChiDT[i] * expMinusChiDT[j] + b[i][j]; //= y[i][j]
            // finally we calculate x
            for (size_t i = 0; i < d_; ++i) {  // E[ x1 | x0 ]
                X1[i] = X0[i];
                for (size_t j = 0; j < d_; ++j) X1[i] += a[i][j] / chi_[j] * OneMinusExpMinusChiDT[j];
                X1[i] *= expMinusChiDT[i];
                ActiveType sumB = 0.0;
                for (size_t j = 0; j < d_; ++j) sumB += b[i][j];
                X1[i] += sumB / chi_[i] * OneMinusExpMinusChiDT[i];
            }
            // we overwrite V by the variance
            for (size_t i = 0; i < d_; ++i)
                for (size_t j = 0; j < d_; ++j)
                    V[i][j] *= (1.0 - expMinusChiDT[i] * expMinusChiDT[j]) / (chi_[i] + chi_[j]);
            MatA L = TemplateAuxilliaries::cholesky(V);
            for (size_t i = 0; i < d_; ++i)
                for (size_t j = 0; j < d_; ++j)
                    X1[i] += L[i][j] * dW[j];  // maybe also exploit that L is lower triangular
            // we don't truncate for the moment
            // finally, we need to update s as well
            ActiveType r0 = 0.0;
            for (size_t i = 0; i < d_; ++i) r0 += X0[i];
            ActiveType r1 = 0.0;
            for (size_t i = 0; i < d_; ++i) r1 += X1[i];
            X1[X1.size() - 1] = X0[X0.size() - 1] + 0.5*(r0 + r1)*dt;
        }

        inline ActiveType numeraire(const DateType t, const VecA& X) {
            State state(X,d_);
            return exp(state.s())/termStructure_->discount(t,true);
        }

        inline ActiveType zeroBond(const DateType t, const DateType T, const VecA& X) {
            State state(X,d_);
            return ZeroBond(t, T, state);
        }

        inline ActiveType zeroBond(const DateType t, const DateType T, const VecA& X, const std::string& alias) { 
            // we might want to introduce a domestic alias and check for consistency
            return zeroBond(t, T, X);  // delegate to model zcb
        }

        // the short rate over an integration time period
        // this is required for drift calculation in multi-asset and hybrid models
        inline ActiveType shortRate(const DateType t0, const DateType dt, const VecA& X0, const VecA& X1) {
            PassiveType B_d = termStructure_->discount(t0,true) / termStructure_->discount(t0 + dt);  // deterministic drift part for r_d
            State s0(X0, d_), s1(X1, d_);
            ActiveType x_av(0.0);
            for (size_t k = 0; k < d_; ++k) x_av += (s0.x(k) + s1.x(k)); // don't forget to divide by 2
            return std::log(B_d) / dt + 0.5 * x_av;
        }

        inline virtual VecA zeroBondVolatility(const DateType t, const DateType T, const VecA& X) { 
            // sigma_r^T = sigma_x^T sqrt(z) 
            // sigma_P^T = G(t,T)^T sigma_x^T sqrt(z)
            State state(X, d_);
            MatA sigmaxT = sigma_xT(0.5*(t+T), state);  // assume constant vol on (t, T)
            VecA sigmaPT(sigmaxT.size(),0.0);
            for (size_t j = 0; j < sigmaPT.size(); ++j) {
                for (size_t i = 0; i < sigmaxT.size(); ++i)
                    sigmaPT[j] += G(i, t, T) * sigmaxT[i][j];
                sigmaPT[j] *= sqrt(state.z());  // don't forget stoch vol scaling
            }
            return sigmaPT;
        }

        inline virtual VecA zeroBondVolatilityPrime(const DateType t, const DateType T, const VecA& X) {
            // sigma_r^T = sigma_x^T sqrt(z) 
            // sigma_P^T = G(t,T)^T sigma_x^T sqrt(z)
            State state(X, d_);
            MatA sigmaxT = sigma_xT(0.5*(t+T), state); // assume constant vol on (t, T)
            VecA sigmaPT(sigmaxT.size(), 0.0);
            for (size_t j = 0; j < sigmaPT.size(); ++j) {
                for (size_t i = 0; i < sigmaxT.size(); ++i) {
                    ActiveType GPrime = exp(-chi_[i] * (T - t));
                    sigmaPT[j] += GPrime * sigmaxT[i][j];
                }
                sigmaPT[j] *= sqrt(state.z());  // don't forget stoch vol scaling
            }
            return sigmaPT;
        }


        virtual std::vector< std::string > stateAliases() {
            std::vector< std::string > aliases(size());
            for (size_t i = 0; i < d_; ++i) aliases[i] = "x_" + std::to_string(i);
            for (size_t i = 0; i < d_; ++i)
                for (size_t j = 0; j < d_; ++j) aliases[d_ + (i*d_) + j] = "y_" + std::to_string(i)+ "_" + std::to_string(j);
            aliases[size() - 2] = "z";
            aliases[size() - 1] = "s";
            return aliases;
        }

        virtual std::vector< std::string > factorAliases() {
            std::vector< std::string > aliases(factors());
            for (size_t i = 0; i < d_; ++i) aliases[i] = "x_" + std::to_string(i);
            aliases[factors() - 1] = "z";
            return aliases;
        }

    };
}

#endif  /* ifndef quantlib_templatequasigaussian2_hpp */
