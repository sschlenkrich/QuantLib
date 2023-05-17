/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2019 Sebastian Schlenkrich

*/

/*! \file assetmodelT.hpp
    \brief (MC) pricing for asset model with stochastic interest rates
               
               X(t) = X0 * exp{x(t)}

            We use extended state variables Y = [ x, (z), mu, volAdj ],
            
               dx(t)     = [mu - 0.5*sigma^2]dt + sigma dW   
               dr_d(t)   = 0 dt (domestic rates)
               dr_f(t)   = 0 dt (foreign rates)
               dz(t)     = 0 dt (stochastic volatility, currently not implemented)
               mu        = r_d - r_f (rates differential, provided exogenously)
               sigma     = sigmaLV + volAdj (hybrid volatility adjuster, provided exogenously)
               
           All methods are template based to allow incorporation of Automatic Differentiation
           tools
           
           This is essentially a component of a hybrid model
*/


#ifndef quantlib_templateassetmodel_hpp
#define quantlib_templateassetmodel_hpp

#include <ql/experimental/templatemodels/stochasticprocessT.hpp>

namespace QuantLib {

    // Declaration of the asset model class
    template <class DateType, class PassiveType, class ActiveType>
    class AssetModelT : public StochasticProcessT<DateType, PassiveType, ActiveType> {
    protected:

        // container class definitions
        typedef std::vector<DateType>                      VecD;
        typedef std::vector<PassiveType>                   VecP; 
        typedef std::vector<ActiveType>                    VecA;
        typedef std::vector< std::vector<DateType> >       MatD;
        typedef std::vector< std::vector<PassiveType> >    MatP;
        typedef std::vector< std::vector<ActiveType> >     MatA;

        // as an initial setup we work with a constant lognormal volatility
        ActiveType X0_;
        ActiveType sigma_;

    public:  

        AssetModelT( const ActiveType X0,
                     const ActiveType sigma )
                     : X0_(X0), sigma_(sigma) {}

        // inspectors
        inline const ActiveType X0()    { return X0_;    }
        inline const ActiveType sigma() { return sigma_; }

        // subset of QL's StochasticProcess interface for X = [ x, y, z, d ] (y row-wise)
        // with dX = a[t,X(t)] dt + b[t,X(t)] dW

        // dimension of Y
        inline size_t size()    { return 1; }
        // stochastic factors of x (and z)
        inline size_t factors() { return 1; }
        // initial values for simulation
        inline VecP initialValues() {
            VecP Y(size(),0.0);
            return Y;
        }

        // a[t,X(t)]
        inline VecA drift( const DateType t, const VecA& Y) {
            VecA a(size(),0.0);
            //a[0] = Y[1] - 0.5*sigma_*sigma_;
            QL_FAIL("AssetModel: drift not implemented. Use evolve.");
            return a;
        }

        // b[t,X(t)]
        inline MatA diffusion( const DateType t, const VecA& Y) {
            MatA b(size(), VecA(factors(),0.0));
            //b[0][0] = sigma_;
            QL_FAIL("AssetModel: diffusion not implemented. Use evolve.");
            // finished
            return b;
        }

        // for quanto adjustment we also need volatility
        inline virtual ActiveType volatility(const DateType t, const VecA& Y) { 
            // this is trivial for lognormal model, but it can be more subtile for LV or SLV model
            QL_REQUIRE(Y.size()==3, "AssetModel: Y.size()==3 required");
            return sigma_ + Y[2];
        }

        // simple Euler step
        inline virtual void evolve(const DateType t0, const VecA& Y0, const DateType dt, const VecD& dW, VecA& Y1) {
            // input state Y0 also carries additional parameters, Y0 = [ x0, mu, volAdj ]
            // output state Y1 only requires x1
            QL_REQUIRE(Y0.size()==3, "AssetModel: Y0.size()==3 required");
            ActiveType sigma = volatility(t0,Y0);
            Y1[0] = Y0[0] + (Y0[1] - 0.5*sigma*sigma)*dt + sigma*dW[0]*std::sqrt(dt);
        }

        // asset calculation is the main purpose of this model
        inline virtual ActiveType asset(const DateType t, const VecA& Y, const std::string& alias) { 
            return X0_ * exp(Y[0]);
        }

        virtual std::vector< std::string > stateAliases() {
            std::vector< std::string > aliases(size());
            aliases[0] = "logS";
            // aliases[1] = "mu";
            return aliases;
        }

        virtual std::vector< std::string > factorAliases() {
            std::vector< std::string > aliases(factors());
            aliases[0] = "logS";
            return aliases;
        }

    };

}

#endif  /* ifndef quantlib_templateassetmodel_hpp */
