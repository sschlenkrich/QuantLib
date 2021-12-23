/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2015 Sebastian Schlenkrich

*/

/*! \file stochasticprocessT.hpp
    \brief define interface for general multi-dimensional stochastic process
               
*/


#ifndef quantlib_templatestochasticprocess_hpp
#define quantlib_templatestochasticprocess_hpp

#include <vector>
#include <ql/shared_ptr.hpp>
#include <ql/types.hpp>
#include <ql/errors.hpp>

namespace QuantLib {


    // Declaration of stochastic process class
    template <class DateType, class PassiveType, class ActiveType>
    class StochasticProcessT : public ext::enable_shared_from_this< StochasticProcessT<DateType,PassiveType,ActiveType> > {
    public:
        // container class definitions
        typedef std::vector<DateType>                      VecD;
        typedef std::vector<PassiveType>                   VecP; 
        typedef std::vector<ActiveType>                    VecA;
        typedef std::vector< std::vector<DateType> >       MatD;
        typedef std::vector< std::vector<PassiveType> >    MatP;
        typedef std::vector< std::vector<ActiveType> >     MatA;

        // subset of QL's StochasticProcess interface for X = [ x, y, z, d ] (y row-wise)
        // with dX = a[t,X(t)] dt + b[t,X(t)] dW

        // dimension of X
        virtual size_t size() = 0;
        // stochastic factors (underlying, volatilities and spreads)
        virtual size_t factors() = 0;
        // initial values for simulation
        virtual VecP initialValues() = 0;
        // a[t,X(t)]
        virtual VecA drift( const DateType t, const VecA& X) = 0;
        // b[t,X(t)]
        virtual MatA diffusion( const DateType t, const VecA& X) = 0;

        // truncate process to its well-defined domain and return true (truncated) or false (not truncated)
        inline virtual bool truncate( const DateType t, VecA& X ) { return false; } // default do nothing

        // integrate X1 = X0 + drift()*dt + diffusion()*dW*sqrt(dt)
        // default implementation
        inline virtual void evolve( const DateType t0, const VecA& X0, const DateType dt, const VecD& dW, VecA& X1 ) {
            // ensure X1 has size of X0
            VecA a = drift(t0, X0);
            MatA b = diffusion(t0, X0);
            for (size_t i=0; i<X1.size(); ++i) {
                X1[i] = 0.0;
                for (size_t j=0; j<dW.size(); ++j) X1[i] += b[i][j]*dW[j];
                X1[i] = X0[i] + a[i]*dt + X1[i]*sqrt(dt);
            }
            truncate( t0+dt, X1 );
            return;
        }

        // we set up a common interface such that models can easily be interchanged
        // the concrete model needs to make sure that it is fit for purpose in a
        // particular application, i.e. implement required methods.

        // the numeraire in the domestic currency used for discounting future payoffs
        inline virtual ActiveType numeraire(const DateType t, const VecA& X)                       { QL_FAIL("StochasticProcessT: numeraire not implemented"); return 0; }

        // a zero coupon bond for the model
        inline virtual ActiveType zeroBond(const DateType t, const DateType T, const VecA& X)      { QL_FAIL("StochasticProcessT: zeroBond not implemented"); return 0; }

        // a domestic/foreign currency zero coupon bond
        inline virtual ActiveType zeroBond(const DateType t, const DateType T, const VecA& X, const std::string& alias) { QL_FAIL("StochasticProcessT: zeroBond with alias not implemented"); return 0; }

        // an asset with (individual) drift and volatility
        inline virtual ActiveType asset(const DateType t, const VecA& X, const std::string& alias) { QL_FAIL("StochasticProcessT: (multi) asset not implemented"); return 0; }

        // the short rate over an integration period
        // this is required for drift calculation in multi-asset and hybrid models
        inline virtual ActiveType shortRate(const DateType t0, const DateType dt, const VecA& X0, const VecA& X1) { QL_FAIL("StochasticProcessT: shortRate not implemented"); return 0; }

        // the expectation E^T in the domestic currency terminal meassure
        // this is required to calculate asset adjusters without knowing the implementation of the model
        inline virtual ActiveType forwardAsset(const DateType t, const DateType T, const VecA& X, const std::string& alias) { QL_FAIL("StochasticProcessT: (multi) forwardAsset not implemented"); return 0; }

        // calculate the local volatility of the log-process of the asset
        // this is required continuous barrier estimation via Brownian Bridge
        inline virtual ActiveType assetVolatility(const DateType t, const VecA& X, const std::string& alias) { QL_FAIL("StochasticProcessT: (multi) assetVolatility not implemented"); return 0; }

        // a (domestic) zero coupon bond volatility for the model
        // this is required e.g. for hybrid model stochastic rates adjustment
        inline virtual VecA zeroBondVolatility(const DateType t, const DateType T, const VecA& X) { QL_FAIL("StochasticProcessT: zeroBondVolatility not implemented"); return VecA(0); }

        // a (domestic) zero coupon bond volatility derivative w.r.t. T for the model
        // this is required e.g. for hybrid model stochastic rates adjustment
        inline virtual VecA zeroBondVolatilityPrime(const DateType t, const DateType T, const VecA& X) { QL_FAIL("StochasticProcessT: zeroBondVolatilityPrime not implemented"); return VecA(0); }

        // the expectation E^Q in the domestic currency risk-neutral meassure
        // this is currently used for commodity payoffs
        inline virtual ActiveType futureAsset(const DateType t, const DateType T, const VecA& X, const std::string& alias) { QL_FAIL("StochasticProcessT: (multi) asset not implemented"); return 0; }

        // we want to keep track of the model details
        virtual std::vector< std::string > stateAliases() { QL_FAIL("StochasticProcessT: stateAliases not implemented"); return std::vector< std::string >(0); }

        virtual std::vector< std::string > factorAliases() { QL_FAIL("StochasticProcessT: factorAliases not implemented"); return std::vector< std::string >(0); }

        // options for integration
        enum VolEvolv {
            FullTruncation         = 0,
            LogNormalApproximation = 1,
            LocalGaussian          = 2,
            Other                  = -1
        };

        // default full truncation
        virtual inline VolEvolv volEvolv() { return FullTruncation; }

    };

    typedef StochasticProcessT<QuantLib::Time,QuantLib::Real,QuantLib::Real> RealStochasticProcess;


}

#endif  /* ifndef quantlib_templatestochasticprocess_hpp */
