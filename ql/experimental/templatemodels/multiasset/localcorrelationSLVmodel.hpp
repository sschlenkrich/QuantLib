/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2019 Cord Harms

*/



#ifndef quantlib_localcorrelationslvmodel_hpp
#define quantlib_localcorrelationslvmodel_hpp


#include <ql/experimental/templatemodels/multiasset/multiassetslvmodel.hpp>
#include <ql/experimental/termstructures/localcorrtermstructure.hpp>
#include <ql/experimental/processes/hestonslvprocess.hpp>

namespace QuantLib {


    /* We model a multi-asset local stochastic volatility model by means of the normalized processes log(S_i/S0) 
           
       d(ln(S_t^i)) = (r_t-q_t-0.5*\L_i^2(t,S_t^i)v_i)  dt + \L_ i(t,S_t^i)sqrt(v_i) dW^(i,S) 
       d(v_i)       = \kappa*(\theta-v_i)               dt + \sigma sqrt(v_i)        dW_t^(i,v)
       dW^i dW^j    = p(t,S_t^1,...,S_t^n,v_1,....,v_n) dt

       cf. J. Guyon, 2013, A new Class of local correlation models, with a_i:=sqrt(v_i)
    
    */
    class LocalCorrelationSLVModel : public MultiAssetSLVModel {
    private:
        Handle<LocalCorrTermStructure> localCorrTermStructure_;
        RealStochasticProcess::MatA corrMatrix_;
    public:
        LocalCorrelationSLVModel(const Handle<YieldTermStructure>&                                        termStructure,
                          const std::vector<std::string>&                                                 aliases,
                          const std::vector<ext::shared_ptr<QuantLib::HestonSLVProcess>>&				  processes,
                          const Handle<LocalCorrTermStructure>&											  localCorrTermStructureAsset);
        
        virtual ~LocalCorrelationSLVModel() = default;

        virtual void evolve(const QuantLib::Time t0, const VecA& X0, const QuantLib::Time dt, const VecD& dW, VecA& X1);
    };

}



#endif  /* ifndef quantlib_localcorrelationslvmodel_hpp */
