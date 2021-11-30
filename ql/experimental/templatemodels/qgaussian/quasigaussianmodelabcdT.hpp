/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2015 Sebastian Schlenkrich

*/

/*! \file quasigaussianmodelabcdT.hpp
    \brief (MC) pricing for multi-factor quasi-Gaussian model with stochastic vol
           use abcd function for model parameters      
*/


#ifndef quantlib_templatequasigaussianabcd_hpp
#define quantlib_templatequasigaussianabcd_hpp

#include <ql/experimental/templatemodels/qgaussian/quasigaussianmodelT.hpp>

namespace QuantLib {

    // Declaration of the quasi-Gaussian model class
    template <class DateType, class PassiveType, class ActiveType>
    class QuasiGaussianModelAbcdT : public QuasiGaussianModelT<DateType, PassiveType, ActiveType> {
        // from base class
        using typename StochasticProcessT<DateType, PassiveType, ActiveType>::VolEvolv;
        using QuasiGaussianModelT<DateType, PassiveType, ActiveType>::d_;
        using QuasiGaussianModelT<DateType, PassiveType, ActiveType>::times_;
        using QuasiGaussianModelT<DateType, PassiveType, ActiveType>::lambda_;
        using QuasiGaussianModelT<DateType, PassiveType, ActiveType>::b_;
        using QuasiGaussianModelT<DateType, PassiveType, ActiveType>::alpha_;
        using QuasiGaussianModelT<DateType, PassiveType, ActiveType>::eta_;
        using QuasiGaussianModelT<DateType, PassiveType, ActiveType>::delta_;
        using QuasiGaussianModelT<DateType, PassiveType, ActiveType>::chi_;
        using QuasiGaussianModelT<DateType, PassiveType, ActiveType>::Gamma_;
        using QuasiGaussianModelT<DateType, PassiveType, ActiveType>::theta_;
        using QuasiGaussianModelT<DateType, PassiveType, ActiveType>::z0_;
        using QuasiGaussianModelT<DateType, PassiveType, ActiveType>::procLimit_;
        using QuasiGaussianModelT<DateType, PassiveType, ActiveType>::termStructure_;
        using QuasiGaussianModelT<DateType, PassiveType, ActiveType>::volEvolv_;
        using QuasiGaussianModelT<DateType, PassiveType, ActiveType>::useSwapRateScaling_;
        using QuasiGaussianModelT<DateType, PassiveType, ActiveType>::factorMatrices;
        using QuasiGaussianModelT<DateType, PassiveType, ActiveType>::maxidx;


    protected:

        // container class definitions
        typedef std::vector<DateType>                      VecD;
        typedef std::vector<PassiveType>                   VecP; 
        typedef std::vector<ActiveType>                    VecA;
        typedef std::vector< std::vector<DateType> >       MatD;
        typedef std::vector< std::vector<PassiveType> >    MatP;
        typedef std::vector< std::vector<ActiveType> >     MatA;


        inline virtual bool checkModelParameters(const bool throwException=true){
            bool ok = true;
            // check yield curve...
            // non-zero dimension
            if (d_<1)               { ok = false; if (throwException) QL_REQUIRE(false,"QuasiGaussianModelAbcd number of factors larger zero required."); }
            // non-empty time-grid
            size_t n=times_.size();
            if (n!=4)               { ok = false; if (throwException) QL_REQUIRE(false,"QuasiGaussianModelAbcd wrong time-grid dimension."); }
            // dimensions of time-dependent parameters
            if (lambda_.size()!=d_) { ok = false; if (throwException) QL_REQUIRE(false,"QuasiGaussianModelAbcd wrong lambda dimension."); }
            if (alpha_.size() !=d_) { ok = false; if (throwException) QL_REQUIRE(false,"QuasiGaussianModelAbcd wrong alpha dimension.");  }
            if (b_.size()     !=d_) { ok = false; if (throwException) QL_REQUIRE(false,"QuasiGaussianModelAbcd wrong b dimension.");      }
            for (size_t k=0; k<d_; ++k) {
                if (lambda_[k].size()!=n) { ok = false; if (throwException) QL_REQUIRE(false,"QuasiGaussianModelAbcd wrong lambda time dimension."); }
                if (alpha_[k].size() !=n) { ok = false; if (throwException) QL_REQUIRE(false,"QuasiGaussianModelAbcd wrong alpha time dimension.");  }
                if (b_[k].size()!=n)      { ok = false; if (throwException) QL_REQUIRE(false,"QuasiGaussianModelAbcd wrong b time dimension.");      }
            }
            if (eta_.size() !=n) { ok = false; if (throwException) QL_REQUIRE(false,"QuasiGaussianModelAbcd wrong eta time dimension.");    }
            // dimensions of time-homogeneous parameters
            if (delta_.size()!=d_)   { ok = false; if (throwException) QL_REQUIRE(false,"QuasiGaussianModelAbcd wrong delta dimension."); }
            if (chi_.size()  !=d_)   { ok = false; if (throwException) QL_REQUIRE(false,"QuasiGaussianModelAbcd wrong chi i-dimension.");   }
            if (Gamma_.size()  !=d_) { ok = false; if (throwException) QL_REQUIRE(false,"QuasiGaussianModelAbcd wrong Gamma j-dimension."); }
            for (size_t k=0; k<d_; ++k) {
                if (Gamma_[k].size() !=d_) { ok = false; if (throwException) QL_REQUIRE(false,"QuasiGaussianModelAbcd wrong Gamma dimension."); }
            }
            // plausible parameter values
            // zero time-grid - SHOULD NOT BE USED FOR abcd PARAMETRISATION
            for (size_t k=0; k<times_.size(); ++k) {
                if (times_[k]!=(DateType)0.0) { ok = false; if (throwException) QL_REQUIRE(false,"QuasiGaussianModelAbcd zero time-grid required."); }
            }
            // non-negative values
            //for (size_t j=0; j<n; ++j) {
            //	for (size_t i=0; i<d_; ++i) {
            //        if (lambda_[i][j]<0.0) { ok = false; if (throwException) QL_REQUIRE(false,"QuasiGaussianModelAbcd lambda>=0 required."); }
            //        if (alpha_[i][j] <0.0) { ok = false; if (throwException) QL_REQUIRE(false,"QuasiGaussianModelAbcd alpha>=0 required.");  }
            //        if (b_[i][j]     <0.0) { ok = false; if (throwException) QL_REQUIRE(false,"QuasiGaussianModelAbcd b>=0 required.");      }
            //	}
            //    if (eta_[j]<0.0) { ok = false; if (throwException) QL_REQUIRE(false,"QuasiGaussianModelAbcd eta>=0 required.");      }
            //}
            // positive/ascending values
            if (delta_[0]<=0.0) { ok = false; if (throwException) QL_REQUIRE(false,"QuasiGaussianModelAbcd delta>0 required."); }
            if (chi_[0]  <=0.0) { ok = false; if (throwException) QL_REQUIRE(false,"QuasiGaussianModelAbcd chi>0 required."); }
            for (size_t k=0; k<d_-1; ++k) {
                if (delta_[k]>=delta_[k+1]) { ok = false; if (throwException) QL_REQUIRE(false,"QuasiGaussianModelAbcd ascending delta values required."); }
                if (chi_[k]  >=chi_[k+1])   { ok = false; if (throwException) QL_REQUIRE(false,"QuasiGaussianModelAbcd ascending chi values required."); }
            }
            // plausible correlation values
            for (size_t i=0; i<d_; ++i) {
                for (size_t j=i; j<d_; ++j) {
                    //if (Gamma_[i][j]<0.0) {
                    //	ok = false;
                    //	if (throwException) QL_REQUIRE(false,"QuasiGaussianModelAbcd Gamma[i][j]>=0 required."); 
                    //}
                    if (i==j) {
                        if (Gamma_[i][j]!=1.0) {
                            ok = false;
                            if (throwException) QL_REQUIRE(false,"QuasiGaussianModelAbcd Gamma[i][i]=1 required."); 
                        }
                    }
                    if (Gamma_[i][j]!=Gamma_[j][i]) {
                        ok = false;
                        if (throwException) QL_REQUIRE(false,"QuasiGaussianModelAbcd Gamma[i][j]=Gamma[j][i] required."); 
                    }
                    //if (i<j) {
                    //	if((Gamma_[i][j-1]<=Gamma_[i][j])|(Gamma_[i+1][j]<=Gamma_[i][j])) {
                    //		ok = false;
                    //		if (throwException) QL_REQUIRE(false,"QuasiGaussianModelAbcd Gamma descending sub-diagonals required."); 
                    //	}
                    //}
                }
            }
            // stochastic vol parameters
            if (theta_<=0.0) { ok = false; if (throwException) QL_REQUIRE(false,"QuasiGaussianModelAbcd theta>0 required."); }
            if (z0_!=1.0)    { ok = false; if (throwException) QL_REQUIRE(false,"QuasiGaussianModelAbcd z0=1 required.");    }
            // adjust stochastic process limits to defaults
            VecP tmp(3,0.0); // [ z-limit, y-limit, x-limit ], default no limit
            for (size_t k=0; k<tmp.size(); ++k) if (k<procLimit_.size()) tmp[k] = (procLimit_[k]<0.0) ? (0.0) : (procLimit_[k]);
            procLimit_ = tmp;
            // finished
            return ok;
        }

        public:  

        // Constructor

        QuasiGaussianModelAbcdT(
            const Handle<YieldTermStructure>& termStructure,
            // number of yield curve factors (excluding stoch. vol)
            const size_t                d,       // (d+1)-dimensional Brownian motion for [x(t), z(t)]^T
            // abcd parameters for parametric curves
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
            const VecP &                procLimit = VecP(0)     // stochastic process limits
            ) { // call base model default constructor
                // initialise members manually...
                termStructure_ = termStructure;
                d_             = d;
                times_         = VecP(4,(DateType)0.0);
                lambda_        = lambda;
                alpha_         = alpha;
                b_             = b;
                eta_           = eta;
                delta_         = delta;
                chi_           = chi;
                Gamma_         = Gamma;
                theta_         = theta;
                z0_            =1.0;
                volEvolv_      = volEvolv;
                procLimit_     = procLimit;
                useSwapRateScaling_ = false;
                // check model parameters for abcd parametrisation
                checkModelParameters();
                // calculate  DfT_, HHfInv_
                factorMatrices();
            }

        // clone the model
        virtual ext::shared_ptr<QuasiGaussianModelT<DateType, PassiveType, ActiveType>> clone() { return ext::shared_ptr<QuasiGaussianModelT<DateType, PassiveType, ActiveType>>(new QuasiGaussianModelAbcdT(*this)); }

        // parameter functions based on abcd parametrisation
        // f(t) = [ a + b*t ] e^{-c*t} + d
        inline virtual ActiveType lambda( const size_t i, const DateType t) { return (lambda_[maxidx(i)][0] + lambda_[maxidx(i)][1]*t)*exp(-lambda_[maxidx(i)][2]*t) + lambda_[maxidx(i)][3]; }
        inline virtual ActiveType alpha ( const size_t i, const DateType t) { return (alpha_[maxidx(i)][0]  + alpha_[maxidx(i)][1] *t)*exp(-alpha_[maxidx(i)][2]*t)  + alpha_[maxidx(i)][3];  }
        inline virtual ActiveType b     ( const size_t i, const DateType t) { return (b_[maxidx(i)][0]      + b_[maxidx(i)][1]     *t)*exp(-b_[maxidx(i)][2]*t)      + b_[maxidx(i)][3];      }
        inline virtual ActiveType eta   ( const DateType t)                 { return (eta_[0]               + eta_[1]              *t)*exp(-eta_[2]*t)               + eta_[3];               }

    };

}

#endif  /* ifndef quantlib_templatequasigaussian_hpp */
