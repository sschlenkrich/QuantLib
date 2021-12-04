/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2017 Cord Harms

*/



#ifndef quantlib_localcorrelationbsmodel_hpp
#define quantlib_localcorrelationbsmodel_hpp


#include <ql/experimental/templatemodels/multiasset/multiassetbsmodel.hpp>
#include <ql/experimental/termstructures/localcorrtermstructure.hpp>
#include <ql/processes/blackscholesprocess.hpp>

namespace QuantLib {


    /* We model a multi-asset local volatility model by means of the normalized log-processes X_i = log[S_i/S_(0)].
       Correlation between assets is modelled state- and time-dependent using a LocalCorrTermStructure: p(t,S_t^1,...,S_t^n).
       cf. J. Guyon, 2013, A new Class of local correlation models
       
       d(ln(S_t^i)) = (r_t-q_t-0.5*\sigma_i^2(t,S_t^i)) dt + \sigma_ i(t,S_t^i) dW^i 
       dW^i dW^j    = p(t,S_t^1,...,S_t^n)              dt
    
    */
    class LocalCorrelationBSModel : public MultiAssetBSModel {
    private:
        Handle<LocalCorrTermStructure> localCorrTermStructure_;
        RealStochasticProcess::MatA corrMatrix_;
    public:
        LocalCorrelationBSModel(const Handle<YieldTermStructure>&                                         termStructure,
                          const std::vector<std::string>&                                                 aliases,
                          const std::vector<ext::shared_ptr<QuantLib::GeneralizedBlackScholesProcess>>&   processes,
                          const Handle<LocalCorrTermStructure>&											  localCorrTermStructure);
        LocalCorrelationBSModel(const Handle<YieldTermStructure>&                                         termStructure,
            const std::vector<std::string>&																  aliases,
            const std::vector<ext::shared_ptr<QuantLib::LocalVolSurface>>&				                  localVolSurfaces,
            const Handle<LocalCorrTermStructure>&											              localCorrTermStructure);

        virtual void evolve(const QuantLib::Time t0, const VecA& X0, const QuantLib::Time dt, const VecD& dW, VecA& X1);
    };

}



#endif  /* ifndef quantlib_localcorrelationbsmodel_hpp */
