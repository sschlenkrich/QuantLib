/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2020 Sebastian Schlenkrich

*/

/*! \file localvolassetmodel.hpp
    \brief (MC) pricing for local volatility asset model with stochastic interest rates

           This is a specialisation of the simple asset model
               
               X(t) = X0 * exp{x(t)}

            We use extended state variables Y = [ x, mu ]
            
               dx(t)     = [mu - 0.5*sigma^2]dt + sigma dW   
               mu        = r_d - r_f (rates differential, provided exogenously)
                       
           This is a component of a hybrid models
*/


#ifndef quantlib_templatelocalvolassetmodel_hpp
#define quantlib_templatelocalvolassetmodel_hpp

#include <ql/termstructures/volatility/equityfx/localvolsurface.hpp>
#include <ql/experimental/templatemodels/hybrid/hybridmodels.hpp>

namespace QuantLib {

    class LocalvolAssetModel : public AssetModel {
    protected:

        // we replace constant lognormal volatility by local vol term structure
        Handle<LocalVolTermStructure> localvol_;

    public:  

        LocalvolAssetModel(
            const Real X0, 
            const Handle<LocalVolTermStructure> localvol )
            : AssetModel(X0, 0.0), localvol_(localvol) {}

        // inspectors
        inline const Handle<LocalVolTermStructure>& localVolatility() const { return localvol_; }

        // for quanto adjustment we also need volatility
        inline virtual Real volatility(const Time t, const std::vector<Real>& Y) { 
            // this is trivial for lognormal model, but it can be more subtile for LV or SLV model
            QL_REQUIRE(Y.size() == 3, "AssetModel: Y.size()==3 required");
            return localvol_->localVol(t, asset(t, Y, ""), true) + Y[2];
        }

    };

}

#endif  /* ifndef quantlib_templatelocalvolassetmodel_hpp */
