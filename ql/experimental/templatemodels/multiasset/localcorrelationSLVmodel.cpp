/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2019 Cord Harms

*/


#include <vector>
#include <ql/errors.hpp>
#include <ql/experimental/templatemodels/auxilliaries/choleskyfactorisationT.hpp>

#include <ql/experimental/templatemodels/multiasset/localcorrelationSLVmodel.hpp>


namespace QuantLib {

    LocalCorrelationSLVModel::LocalCorrelationSLVModel(
        const Handle<YieldTermStructure>&                                               termStructure,
        const std::vector<std::string>&                                                 aliases,
        const std::vector<ext::shared_ptr<QuantLib::HestonSLVProcess>>&				    processes,
        const Handle<LocalCorrTermStructure>&											localCorrTermStructure)
        : MultiAssetSLVModel(termStructure, aliases, processes),
        localCorrTermStructure_(localCorrTermStructure) {
        corrMatrix_ = getPureHestonImpliedCorrelationMatrix(); //init corrMatrix, it is overwritten by term structure in evolve-function.
    }

    void LocalCorrelationSLVModel::evolve(const QuantLib::Time t0, const VecA& X0, const QuantLib::Time dt, const VecD& dW, VecA& X1) {
        // need to calculate the time- and state-dependent correlation matrix
        localCorrTermStructure_->localCorr(corrMatrix_, t0, X0 , true); //true because first X0 are always 0 in simulation.
        //DT_ = TemplateAuxilliaries::svdSqrt(corrMatrix_);
        for (size_t i = 0; i < corrMatrix_.size(); i++)
        {
            for (size_t j = i; j < corrMatrix_.size(); j++)
            {
                DT_[i][j] = corrMatrix_[i][j];
                DT_[j][i] = corrMatrix_[j][i];
            }
        }
        
        TemplateAuxilliaries::performCholesky(DT_, corrMatrix_.size(),true);
        
        //now we call the function of multiassetmodel
        MultiAssetSLVModel::evolve(t0, X0, dt, dW, X1);
    }

}



