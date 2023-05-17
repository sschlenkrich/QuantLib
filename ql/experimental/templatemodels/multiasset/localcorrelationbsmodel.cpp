/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2017 Cord Harms

*/

#include <vector>
#include <ql/errors.hpp>
#include <ql/experimental/templatemodels/auxilliaries/choleskyfactorisationT.hpp>

#include <ql/experimental/templatemodels/multiasset/localcorrelationbsmodel.hpp>



namespace QuantLib {

    LocalCorrelationBSModel::LocalCorrelationBSModel(
        const Handle<YieldTermStructure>&                                               termStructure,
        const std::vector<std::string>&                                                 aliases,
        const std::vector<ext::shared_ptr<QuantLib::GeneralizedBlackScholesProcess>>&   processes,
        const Handle<LocalCorrTermStructure>& localCorrTermStructure)
    : MultiAssetBSModel(termStructure,aliases,processes), 
        localCorrTermStructure_(localCorrTermStructure) {
        corrMatrix_ = RealStochasticProcess::MatA(processes.size());

        for (size_t k = 0; k<processes.size(); ++k) corrMatrix_[k].resize(processes.size());
    }

    LocalCorrelationBSModel::LocalCorrelationBSModel(const Handle<YieldTermStructure>&                termStructure,
        const std::vector<std::string>&																  aliases,
        const std::vector<ext::shared_ptr<QuantLib::LocalVolSurface>>&				                  localVolSurfaces,
        const Handle<LocalCorrTermStructure>&											              localCorrTermStructure) 
        : MultiAssetBSModel(termStructure, aliases, localVolSurfaces),
        localCorrTermStructure_(localCorrTermStructure)
    {
        corrMatrix_ = RealStochasticProcess::MatA(processes_.size());

        for (size_t k = 0; k<processes_.size(); ++k) corrMatrix_[k].resize(processes_.size());
    }


    void LocalCorrelationBSModel::evolve(const QuantLib::Time t0, const VecA& X0, const QuantLib::Time dt, const VecD& dW, VecA& X1) {
        // we merely need to calculate the correlation matrix
        localCorrTermStructure_->localCorr(corrMatrix_, t0, X0 , true); //true because first X0 are always 0 in simulation.
        //DT_ = TemplateAuxilliaries::svdSqrt(corrMatrix_);
        TemplateAuxilliaries::performCholesky(corrMatrix_, corrMatrix_.size(),true);
        DT_ = corrMatrix_;
        //now we call the function of multiassetmodel
        MultiAssetBSModel::evolve(t0, X0, dt, dW, X1);
    }

}



