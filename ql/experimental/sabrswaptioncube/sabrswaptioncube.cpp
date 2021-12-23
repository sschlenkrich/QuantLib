/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2010 Sebastian Schlenkrich
*/

/*! \file bondoptionengine.cpp
    \brief engine for bond options with Hull White model
*/



#include <ql/experimental/sabrswaptioncube/sabrswaptioncube.hpp>


namespace QuantLib {

    SabrSwaptionCube::SabrSwaptionCube (
                           const std::vector<Time>&                   optionTimes,
                           const std::vector<Time>&                   swapTimes,
                           const std::vector< std::vector< Real > >&  alpha,
                           const std::vector< std::vector< Real > >&  beta,
                           const std::vector< std::vector< Real > >&  rho,
                           const std::vector< std::vector< Real > >&  nu,
                           const std::vector< std::vector< Real > >&  fwd,
                           BusinessDayConvention                      bdc,
                           const DayCounter&                          dc,
                           const bool                                 useNormalVols)
    : SwaptionVolatilityStructure(bdc, dc), 
      optionTimes_(optionTimes), swapTimes_(swapTimes),
      alpha_(optionTimes.size(),swapTimes.size(),0.0), 
      beta_(optionTimes.size(),swapTimes.size(),0.0), 
      rho_(optionTimes.size(),swapTimes.size(),0.0),  
      nu_(optionTimes.size(),swapTimes.size(),0.0), 
      fwd_(optionTimes.size(),swapTimes.size(),0.0), 
      maxSwapTenor_(100*Years),
      referenceDate_(Settings::instance().evaluationDate()),
      useNormalVols_(useNormalVols) {
        // transform input data
        for (Size i=0; i<optionTimes.size(); ++i) {
            for (Size j=0; j<swapTimes.size(); ++j) {
                if (i<alpha.size() && j<alpha[i].size()) alpha_[i][j] = alpha[i][j];
                if (i<beta.size()  && j<beta[i].size())  beta_[i][j]  = beta[i][j];
                if (i<rho.size()   && j<rho[i].size())   rho_[i][j]   = rho[i][j];
                if (i<nu.size()    && j<nu[i].size())    nu_[i][j]    = nu[i][j];
                if (i<fwd.size()   && j<fwd[i].size())   fwd_[i][j]   = fwd[i][j];
            }
        }
        // set up bilinear interpolators
        ext::shared_ptr<Interpolation2D> interpolation;
        alphaInterp_ = ext::shared_ptr<Interpolation2D>( new FlatExtrapolator2D( ext::shared_ptr<Interpolation2D>(new BilinearInterpolation (swapTimes_.begin(), swapTimes_.end(), optionTimes_.begin(), optionTimes_.end(), alpha_))));
        betaInterp_  = ext::shared_ptr<Interpolation2D>( new FlatExtrapolator2D( ext::shared_ptr<Interpolation2D>(new BilinearInterpolation (swapTimes_.begin(), swapTimes_.end(), optionTimes_.begin(), optionTimes_.end(), beta_))));
        rhoInterp_   = ext::shared_ptr<Interpolation2D>( new FlatExtrapolator2D( ext::shared_ptr<Interpolation2D>(new BilinearInterpolation (swapTimes_.begin(), swapTimes_.end(), optionTimes_.begin(), optionTimes_.end(), rho_))));
        nuInterp_    = ext::shared_ptr<Interpolation2D>( new FlatExtrapolator2D( ext::shared_ptr<Interpolation2D>(new BilinearInterpolation (swapTimes_.begin(), swapTimes_.end(), optionTimes_.begin(), optionTimes_.end(), nu_))));
        fwdInterp_   = ext::shared_ptr<Interpolation2D>( new FlatExtrapolator2D( ext::shared_ptr<Interpolation2D>(new BilinearInterpolation (swapTimes_.begin(), swapTimes_.end(), optionTimes_.begin(), optionTimes_.end(), fwd_))));
        alphaInterp_->enableExtrapolation();
        betaInterp_->enableExtrapolation();
        rhoInterp_->enableExtrapolation();
        nuInterp_->enableExtrapolation();
        fwdInterp_->enableExtrapolation();
    }


    ext::shared_ptr<SmileSection> SabrSwaptionCube::smileSectionImpl(
                           Time optionTime,
                           Time swapLength) const {
        // interpolate and set up smile section
        std::vector<Real> sabrParameters(4);
        sabrParameters[0] = (*alphaInterp_)(swapLength, optionTime);
        sabrParameters[1] = (*betaInterp_)(swapLength, optionTime);
        sabrParameters[2] = (*nuInterp_)(swapLength, optionTime);
        sabrParameters[3] = (*rhoInterp_)(swapLength, optionTime);
        Real fwd = (*fwdInterp_)(swapLength, optionTime);
        return ext::shared_ptr<SmileSection>( new SabrSmileSection(optionTime, fwd, sabrParameters,0.0,(useNormalVols_)?(VolatilityType::Normal):(VolatilityType::Lognormal)));
    }

}

