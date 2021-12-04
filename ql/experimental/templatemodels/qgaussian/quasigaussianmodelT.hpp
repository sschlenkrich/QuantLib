/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2015 Sebastian Schlenkrich

*/

/*! \file quasigaussianmodelT.hpp
    \brief (MC) pricing for multi-factor quasi-Gaussian model with stochastic vol
               
               r(t) = f(0,t) + 1^T*x(t)
               
               dx(t)     = [ y(t)*1 - a*x(t) ] dt                                        + sqrt[z(t)]*sigma_x^T(t,x,y) dW
               dy(t)     = [ z(t)*sigma_x^T(t,x,y)*sigma_x(t,x,y) - a*y(t) - y(t)*a ] dt
               dz(t)     = theta [ z0 - z(t) ] dt                                        + eta(t)*sqrt[z(t)]           dZ
               ds(t)     = r(t) dt  ( s(t) = int_0^t r(s) ds, for bank account numeraire)
               
           All methods are template based to allow incorporation of Automatic Differentiation
           tools
*/


#ifndef quantlib_templatequasigaussian_hpp
#define quantlib_templatequasigaussian_hpp

#include <ql/termstructures/yieldtermstructure.hpp>
#include <ql/option.hpp>
#include <ql/experimental/templatemodels/stochasticprocessT.hpp>
#include <ql/experimental/templatemodels/auxilliaries/auxilliariesT.hpp>
#include <ql/experimental/templatemodels/auxilliaries/integratorsT.hpp>
#include <ql/experimental/templatemodels/auxilliaries/svdT.hpp>
#include <ql/experimental/templatemodels/auxilliaries/choleskyfactorisationT.hpp>



namespace QuantLib {

    // Declaration of the quasi-Gaussian model class
    template <class DateType, class PassiveType, class ActiveType>
    class QuasiGaussianModelT : public StochasticProcessT<DateType, PassiveType, ActiveType> {
        // from base class
        using typename StochasticProcessT<DateType, PassiveType, ActiveType>::VolEvolv;


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
        MatA                       lambda_;  // volatility
        MatA                       alpha_;   // shift
        MatA                       b_;       // f-weighting
        VecA                       eta_;     // vol-of-vol
        // scaling parameters
        bool                       useSwapRateScaling_;
        MatA                       S0_;
        MatA                       D_;
        // time-homogeneous parameters
        VecP                       delta_;   // maturity of benchmark rates f(t,t+delta_i) 		
        VecP                       chi_;     // mean reversions
        MatP                       Gamma_;   // (benchmark rate) correlation matrix
        // stochastic volatility process parameters
        PassiveType                theta_;   // mean reversion speed
        PassiveType                z0_;      // mean reversion level z0=z(0)=1
        VolEvolv                   volEvolv_; // use FullTruncation or LogNormalApproximation for volatility process integration

        // truncate stochastic process to avoid numerical instabilities (fpn overflow)
        VecP                       procLimit_;  // [ z-limit, y-limit, x-limit ]; lower/upper limit for x, y; upper limit for z

        // additional parameters (calculated at initialisation via SVD)
        MatP                       DfT_;     // factorized correlation matrix Df^T with Df^T * Df = Gamma
        MatP                       HHfInv_;  // weighting matrix H*Hf^-1 = [ exp{-chi_j*delta_i} ]^-1

        // lightweight container holding the current state of the yield curve
        class State {
        public:
            VecA        x;
            MatA        y;
            ActiveType  z;
            ActiveType  s;
            // constructor
            State( const VecA& X, const size_t d) {
                QL_REQUIRE(X.size()==d+d*d+1+1,"TemplateQuasiGaussianModel::State Constructor: Dimensions mismatch.");
                x.resize(d);
                y.resize(d);
                for (size_t k=0; k<d; ++k) x[k] = X[k];
                for (size_t i=0; i<d; ++i) {
                    y[i].resize(d);
                    for (size_t j=0; j<d; ++j) y[i][j] = X[d+i*d+j];  // y row-wise
                }
                z = X[d+d*d  ];
                s = X[d+d*d+1];
            }
            inline void toVec(VecA& X) {
                size_t d = x.size();
                // check dimensions
                QL_REQUIRE(y.size()==d,"TemplateQuasiGaussianModel::State Assignment: y-row dimension mismatch.");
                for (size_t k=0; k<d; ++k) {
                    QL_REQUIRE(y[k].size()==d,"TemplateQuasiGaussianModel::State Assignment: y-column dimension mismatch.");
                }
                X.resize(d+d*d+1+1);
                for (size_t k=0; k<d; ++k) X[k]       = x[k];
                // this was wrong: for (size_t k=0; k<d; ++k) X[d+i*d+j] = y[i][j];  
                for (size_t i=0; i<d; ++i)
                    for (size_t j=0; j<d; ++j) X[d+i*d+j] = y[i][j];
                X[d+d*d  ]                            = z;
                X[d+d*d+1]                            = s;
            }
        };

        inline virtual bool checkModelParameters(const bool throwException=true){
            bool ok = true;
            // check yield curve...
            // non-zero dimension
            if (d_<1)               { ok = false; if (throwException) QL_REQUIRE(false,"QuasiGaussianModel number of factors larger zero required."); }
            // non-empty time-grid
            size_t n=times_.size();
            if (n<1)                { ok = false; if (throwException) QL_REQUIRE(false,"QuasiGaussianModel non-empty time-grid required."); }
            // dimensions of time-dependent parameters
            if (lambda_.size()!=d_) { ok = false; if (throwException) QL_REQUIRE(false,"QuasiGaussianModel wrong lambda dimension."); }
            if (alpha_.size() !=d_) { ok = false; if (throwException) QL_REQUIRE(false,"QuasiGaussianModel wrong alpha dimension.");  }
            if (b_.size()     !=d_) { ok = false; if (throwException) QL_REQUIRE(false,"QuasiGaussianModel wrong b dimension.");      }
            for (size_t k=0; k<d_; ++k) {
                if (lambda_[k].size()!=n) { ok = false; if (throwException) QL_REQUIRE(false,"QuasiGaussianModel wrong lambda time dimension."); }
                if (alpha_[k].size() !=n) { ok = false; if (throwException) QL_REQUIRE(false,"QuasiGaussianModel wrong alpha time dimension.");  }
                if (b_[k].size()!=n)      { ok = false; if (throwException) QL_REQUIRE(false,"QuasiGaussianModel wrong b time dimension.");      }
            }
            if (eta_.size() !=n) { ok = false; if (throwException) QL_REQUIRE(false,"QuasiGaussianModel wrong eta time dimension.");    }
            // dimensions of time-homogeneous parameters
            if (delta_.size()!=d_)   { ok = false; if (throwException) QL_REQUIRE(false,"QuasiGaussianModel wrong delta dimension."); }
            if (chi_.size()  !=d_)   { ok = false; if (throwException) QL_REQUIRE(false,"QuasiGaussianModel wrong chi i-dimension.");   }
            if (Gamma_.size()  !=d_) { ok = false; if (throwException) QL_REQUIRE(false,"QuasiGaussianModel wrong Gamma j-dimension."); }
            for (size_t k=0; k<d_; ++k) {
                if (Gamma_[k].size() !=d_) { ok = false; if (throwException) QL_REQUIRE(false,"QuasiGaussianModel wrong Gamma dimension."); }
            }
            // plausible parameter values
            // ascending time-grid
            for (size_t k=0; k<times_.size()-1; ++k) {
                if (times_[k]>=times_[k+1]) { ok = false; if (throwException) QL_REQUIRE(false,"QuasiGaussianModel ascending time-grid required."); }
            }
            // non-negative values
            for (size_t j=0; j<n; ++j) {
                for (size_t i=0; i<d_; ++i) {
                    if (lambda_[i][j]<0.0) { ok = false; if (throwException) QL_REQUIRE(false,"QuasiGaussianModel lambda>=0 required."); }
                    if (alpha_[i][j] <0.0) { ok = false; if (throwException) QL_REQUIRE(false,"QuasiGaussianModel alpha>=0 required.");  }
                    if (b_[i][j]     <0.0) { ok = false; if (throwException) QL_REQUIRE(false,"QuasiGaussianModel b>=0 required.");      }
                }
                if (eta_[j]<0.0) { ok = false; if (throwException) QL_REQUIRE(false,"QuasiGaussianModel eta>=0 required.");      }
            }
            // positive/ascending values
            if (delta_[0]<=0.0) { ok = false; if (throwException) QL_REQUIRE(false,"QuasiGaussianModel delta>0 required."); }
            if (chi_[0]  <=0.0) { ok = false; if (throwException) QL_REQUIRE(false,"QuasiGaussianModel chi>0 required."); }
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
            // adjust stochastic process limits to defaults
            VecP tmp(3,0.0); // [ z-limit, y-limit, x-limit ], default no limit
            for (size_t k=0; k<tmp.size(); ++k) if (k<procLimit_.size()) tmp[k] = (procLimit_[k]<0.0) ? (0.0) : (procLimit_[k]);
            procLimit_ = tmp;
            // finished
            return ok;
        }

        // evaluate Df^T with Df^T * Df = Gamma and H*Hf^-1 via singular value decomposition
        // return false (and throw exception) on error
        inline bool factorMatrices(const bool throwException=true){
            bool ok = true;
            // row-wise matrices
            size_t dim = d_;
            PassiveType *A  = new PassiveType[dim*dim];
            PassiveType *U  = new PassiveType[dim*dim];
            PassiveType *S  = new PassiveType[dim];
            PassiveType *VT = new PassiveType[dim*dim];
            PassiveType minS;
            // dummy auxilliary variables
            PassiveType work;
            int lwork, info;
            DfT_ = TemplateAuxilliaries::cholesky(Gamma_);

            // [ Hf H^{-1} ] = [ exp{-chi_j*delta_i} ] = V^T S U
            for (size_t i=0; i<dim; ++i) {
                for (size_t j=0; j<dim; ++j) {
                    A[i*dim+j] = exp(-chi_[j]*delta_[i]);
                }
            }
            TemplateAuxilliaries::svd("S","S",(int*)&dim,(int*)&dim,A,(int*)&dim,S,U,(int*)&dim,VT,(int*)&dim,&work,&lwork,&info);
            // check min(S)>0
            minS=S[0];
            for (size_t i=1; i<dim; ++i) if (S[i]<minS) minS = S[i];
            if (minS<=0) { ok = false; if (throwException) QL_REQUIRE(false,"QuasiGaussianModel non-singular Gamma required."); }
            // evaluate H*Hf^-1 = U^T S^{-1} V
            HHfInv_.resize(dim);
            for (size_t i=0; i<dim; ++i) {
                HHfInv_[i].resize(dim);
                for (size_t j=0; j<dim; ++j) {
                    HHfInv_[i][j] = 0.0;
                    for (size_t k=0; k<dim; ++k) HHfInv_[i][j] += U[k*dim+i] * VT[j*dim+k] / S[k];
                }
            }
            // finished
            delete[] A;
            delete[] U;
            delete[] S;
            delete[] VT;
            return ok;
        }

        // required for swap rate gradient calculation
        inline VecA zcbGradient( const DateType T) {
            VecA grad(d_);
            ActiveType DF = termStructure_->discount(T);
            for (size_t k=0; k<grad.size(); ++k) grad[k] = -DF * G(k,0.0,T);
            return grad;
        }

        // simplified swap rate and gradient evaluation for model parameter scaling
        inline VecA swapGrad(const DateType T0, const DateType TN, PassiveType& swapRate) {
            PassiveType num = termStructure_->discount(T0) - termStructure_->discount(TN);
            PassiveType den = 0.0;
            for (PassiveType Ti = T0; Ti<TN; Ti+=1.0) {
                PassiveType T = (Ti+1.0>TN) ? (TN) : (Ti+1.0);
                den += (T-Ti) * termStructure_->discount(T);
            }
            swapRate = num / den;
            // gradient evaluation
            VecA grad = zcbGradient(T0);
            VecA gZCB = zcbGradient(TN);
            for (size_t k=0; k<grad.size(); ++k) grad[k] = (grad[k] - gZCB[k]) / den;
            ActiveType dSdAnnuity = - swapRate / den;
            for (PassiveType Ti = T0; Ti<TN; Ti+=1.0) {
                PassiveType T = (Ti+1.0>TN) ? (TN) : (Ti+1.0);
                gZCB = zcbGradient(T);
                for (size_t k=0; k<grad.size(); ++k) grad[k] += dSdAnnuity * (T-Ti) * gZCB[k];
            }
            return grad;
        }

        inline void rescaleAlphaB() {
            S0_.resize(d_);
            D_.resize(d_);
            for (size_t k=0; k<d_; ++k) {
                S0_[k].resize(times_.size());
                D_[k].resize(times_.size());
                for (size_t i=0; i<times_.size(); ++i) {
                    VecA grad = swapGrad(times_[i],times_[i]+delta_[k],S0_[k][i]);
                    D_[k][i] = 0.0;
                    for (size_t j=0; j<d_; ++j) D_[k][i] += grad[j] * exp(-chi_[j]*delta_[i]);
                    // overwrite alpha, b
                    // alpha_[k][i] = (1.0 - b_[k][i])*S0_[k][i];
                    // b_[k][i] *= D_[k][i];
                }
            }
        }


        public:  


        // Constructor
        QuasiGaussianModelT() { } // do nothing; (unsafe)

        QuasiGaussianModelT(
            const Handle<YieldTermStructure>& termStructure,
            // number of yield curve factors (excluding stoch. vol)
            const size_t                d,       // (d+1)-dimensional Brownian motion for [x(t), z(t)]^T
            // unique grid for time-dependent parameters
            const VecD &                times,   // time-grid of left-constant model parameter values
            // time-dependent parameters, left-piecewise constant on times_-grid
            const MatA &                lambda,  // volatility
            const MatA &                alpha,   // shift
            const MatA &                b,       // f-weighting
            const VecA &                eta,     // vol-of-vol
            // time-homogeneous parameters
            const VecP &                delta,   // maturity of benchmark rates f(t,t+delta_i) 		
            const VecP &                chi,     // mean reversions
            const MatP &                Gamma,   // (benchmark rate) correlation matrix
            // stochastic volatility process parameters
            const PassiveType           theta,   // mean reversion speed
            const VolEvolv              volEvolv  = VolEvolv::FullTruncation,
            const VecP &                procLimit = VecP(0),     // stochastic process limits
            const bool                  useSwapRateScaling = true
            ) : termStructure_(termStructure), d_(d), times_(times), lambda_(lambda), alpha_(alpha), b_(b), eta_(eta),
                delta_(delta), chi_(chi), Gamma_(Gamma), theta_(theta), z0_((PassiveType)1.0), volEvolv_(volEvolv), procLimit_(procLimit), useSwapRateScaling_(useSwapRateScaling) {
                checkModelParameters();
                // calculate  DfT_
                // calculate  HHfInv_
                factorMatrices();
                // adjust alpha and beta to approximate swap dynamics
                if (useSwapRateScaling_) rescaleAlphaB();
            }

        // update model parameters (e.g. during calibration)
        void update( const MatA &                lambda,  // volatility
                     const MatA &                b,       // f-weighting
                     const VecA &                eta ) {  // vol-of-vol
            // perform some checks on inputs...
            lambda_ = lambda;
            b_      = b;
            eta_    = eta;
            if (useSwapRateScaling_) rescaleAlphaB();
        }

        // helpers

        inline size_t maxidx( const size_t i ) { return (i<d_) ? i : d_-1; }

        inline size_t idx( const DateType t ) { return TemplateAuxilliaries::idx(times_,t); }

        // clone the model
        virtual ext::shared_ptr<QuasiGaussianModelT> clone() { return ext::shared_ptr<QuasiGaussianModelT>(new QuasiGaussianModelT(*this)); }

        // inspectors
        inline const Handle<YieldTermStructure> termStructure() { return termStructure_; }
        inline const VecD& times()   { return times_;  }
        inline const MatA& lambda()  { return lambda_; }
        inline const MatA& alpha()   { return alpha_;  }
        inline const MatA& b()       { return b_;      }
        inline const VecA& eta()     { return  eta_;   }


        inline const MatP& DfT()    { return DfT_;    }
        inline const MatP& HHfInv() { return HHfInv_; }
        inline const VecP& delta()  { return delta_;  }
        inline const VecP& chi()    { return chi_;    }
        inline const ActiveType theta() { return theta_; } 
        inline const ActiveType z0()    { return z0_;    } 

        inline VolEvolv volEvolv()  { return volEvolv_; }  

        // parameter functions (no dimension checks)
        inline virtual ActiveType lambda( const size_t i, const DateType t) { return lambda_[maxidx(i)][idx(t)]; }
        inline virtual ActiveType alpha ( const size_t i, const DateType t) { return (useSwapRateScaling_) ? ((1.0-b_[maxidx(i)][idx(t)])*S0_[maxidx(i)][idx(t)]) : (alpha_[maxidx(i)][idx(t)]) ;  }
        inline virtual ActiveType b     ( const size_t i, const DateType t) { return (useSwapRateScaling_) ? (D_[maxidx(i)][idx(t)]*b_[maxidx(i)][idx(t)])        : (b_[maxidx(i)][idx(t)])     ;  }
        inline virtual ActiveType eta   ( const DateType t)                 { return eta_[idx(t)];       }

        // analytic formulas

        inline ActiveType G(const size_t i, const DateType t, const DateType T) { return (1.0-exp(-chi_[i]*(T-t)))/chi_[i]; }

        inline virtual
        ActiveType shortRate ( const DateType t, const VecA& x ) {
            ActiveType r = termStructure_->forwardRate(t,t,Continuous);
            for (size_t k=0; k<d_; ++k) r += x[k];
            return r;
        }

        inline virtual
        ActiveType forwardRate( const DateType t, const DateType T, const VecA& x, const MatA&  y) {
            ActiveType f = termStructure_->forwardRate(T,T,Continuous);  // check t,T
            for (size_t i=0; i<d_; ++i) {
                ActiveType tmp = x[i];
                for (size_t j=0; j<d_; ++j) tmp += y[i][j]*G(j,t,T);
                f += exp(-chi_[i]*(T-t)) * tmp;
            }
            return f;
        }

        inline virtual
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
            for (size_t k=0; k<d_; ++k) res[k] = lambda(k,t) * (alpha(k,t) + b(k,t)*forwardRate(t,t+delta_[k],x,y));
            return res;
        }

        inline  // sigma_x^T
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
        inline virtual size_t size()    { return d_ + d_*d_ + 1 + 1; }
        // stochastic factors of x and z (maybe distinguish if trivially eta=0)
        inline virtual size_t factors() { return d_ + 1; }
        // initial values for simulation
        inline virtual VecP initialValues() {
            VecP X(size());
            for (size_t k=0; k<d_ + d_*d_; ++k) X[k] = 0.0;  // x(0), y(0)
            X[d_+d_*d_]                              = 1.0;  // z(0)
            X[d_+d_*d_+1]                            = 0.0;  // s(0)
            return X;
        }

        // a[t,X(t)]
        inline virtual VecA drift( const DateType t, const VecA& X) {
            VecA a(size());
            State state(X,d_);
            // x-variable [ y(t)*1 - chi*x(t) ]
            for (size_t k=0; k<d_; ++k) {
                a[k] = -chi_[k]*state.x[k];
                for (size_t j=0; j<d_; ++j) a[k] += state.y[k][j];
            }
            // y-variable [ z(t)*sigma_x^T(t,x,y)*sigma_x(t,x,y) - chi*y(t) - y(t)*chi ]
            MatA sigmaxT = sigma_xT(t,state.x,state.y);
            for (size_t i=0; i<d_; ++i) {
                for (size_t j=0; j<d_; ++j) {
                    a[d_+i*d_+j] = 0.0;
                    for (size_t k=0; k<d_; ++k) a[d_+i*d_+j] += sigmaxT[i][k]*sigmaxT[k][j];
                    a[d_+i*d_+j] *= ((state.z>0)?(state.z):(0.0));  // full truncation
                    a[d_+i*d_+j] -= (chi_[i]+chi_[j])*state.y[i][j];
                }
            }
            // z-variable theta [ z0 - z(t)^+ ]  (full truncation)
            //a[d_+d_*d_] = theta_*(z0_ - ((state.z>0)?(state.z):(0.0)));
            // push to positive teritory
            a[d_+d_*d_] = theta_*(z0_ - state.z);
            // s-variable r(t)
            a[d_+d_*d_+1] = shortRate(t,state.x);
            // finished
            return a;
        }

        // b[t,X(t)]
        inline virtual MatA diffusion( const DateType t, const VecA& X) {
            MatA b(size());
            for (size_t k=0; k<size(); ++k) b[k].resize(factors());
            State state(X,d_);
            ActiveType sqrtz = ((state.z>0)?(sqrt(state.z)):(0.0));   // full truncation
            MatA sigmaxT = sigma_xT(t,state.x,state.y);
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

        // integrate X1 = X0 + drift()*dt + diffusion()*dW*sqrt(dt)
        inline virtual void evolve( const DateType t0, const VecA& X0, const DateType dt, const VecD& dW, VecA& X1 ) {
            if (volEvolv() == VolEvolv::LocalGaussian) {
                evolveAsLocalGaussian(t0, X0, dt, dW, X1);  // this is not really good coding style
                return;
            }
            // ensure X1 has size of X0
            VecA a = drift(t0, X0);
            MatA b = diffusion(t0, X0);
            // default via drift and diffusion
            for (size_t i=0; i<X1.size(); ++i) {
                X1[i] = 0.0;
                for (size_t j=0; j<dW.size(); ++j) X1[i] += b[i][j]*dW[j];
                X1[i] = X0[i] + a[i]*dt + X1[i]*sqrt(dt);
            }
            if (volEvolv()==VolEvolv::FullTruncation) {
                if (X1[X1.size()-2]<0.0) X1[X1.size()-2] = 0.0;
            }
            if (volEvolv()==VolEvolv::LogNormalApproximation) {
                ActiveType e = expectationZ(t0, X0[X0.size()-2], dt);
                ActiveType v = varianceZ(t0, X0[X0.size()-2], dt);
                //ActiveType dZ = dW[dW.size()-2];  // last risk factor is for vol process
                ActiveType dZ = dW[dW.size()-1];  // last risk factor is for vol process
                ActiveType si = sqrt(log(1.0 + v/e/e));
                ActiveType mu = log(e) - si*si/2.0;
                X1[X1.size()-2] = exp(mu + si*dZ);
            }
            truncate( t0+dt, X1 );
            return;
        }

        // simulate Quasi-Gaussian model as Gaussian model with frozen vol matrix
        inline void evolveAsLocalGaussian(const DateType t0, const VecA& X0, const DateType dt, const VecD& dW, VecA& X1) {
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
            MatA sigmaxT = sigma_xT(t0, state.x, state.y);
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
            ActiveType r0 = termStructure_->forwardRate(t0, t0, Continuous);
            for (size_t i = 0; i < d_; ++i) r0 += X0[i];
            ActiveType r1 = termStructure_->forwardRate(t0+dt, t0+dt, Continuous);
            for (size_t i = 0; i < d_; ++i) r1 += X1[i];
            X1[X1.size() - 1] = X0[X0.size() - 1] + 0.5*(r0 + r1)*dt;
        }


        // truncate process to its well-defined domain and return true (truncated) or false (not truncated)
        inline virtual bool truncate( const DateType t, VecA& X ) {
            // ensure X.size()==size()
            bool trunc = false;
            if (procLimit_[2]>0.0) {
                for (size_t k=0; k<d_; ++k) { // state x
                    if (X[k]<-procLimit_[2]) { X[k]=-procLimit_[2]; trunc = true; }
                    if (X[k]> procLimit_[2]) { X[k]= procLimit_[2]; trunc = true; }
                }
            }
            if (procLimit_[1]>0.0) {
                for (size_t k=d_; k<d_+d_*d_; ++k) { // state y
                    if (X[k]<-procLimit_[1]) { X[k]=-procLimit_[1]; trunc = true; }
                    if (X[k]> procLimit_[1]) { X[k]= procLimit_[1]; trunc = true; }
                }
            }			
            if (procLimit_[0]>0.0) { // state s
                if (X[d_+d_*d_] < 0.0)           { X[d_+d_*d_] = 0.0;    trunc = true; }
                if (X[d_+d_*d_] > procLimit_[0]) { X[d_+d_*d_] = procLimit_[0]; trunc = true; }
            }
            return trunc;
        }


        inline virtual ActiveType numeraire(const DateType t, const VecA& X) {
            State state(X,d_);
            return exp(state.s);
        }

        inline virtual ActiveType zeroBond(const DateType t, const DateType T, const VecA& X) {
            State state(X,d_);
            return ZeroBond(t, T, state.x, state.y);
        }

    };

}

#endif  /* ifndef quantlib_templatequasigaussian_hpp */
