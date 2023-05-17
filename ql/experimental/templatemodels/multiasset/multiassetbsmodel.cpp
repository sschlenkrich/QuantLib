/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2017 Sebastian Schlenkrich

*/


#include <vector>
#include <ql/errors.hpp>
#include <ql/experimental/templatemodels/auxilliaries/choleskyfactorisationT.hpp>

#include <ql/experimental/templatemodels/multiasset/multiassetbsmodel.hpp>


namespace QuantLib {

    MultiAssetBSModel::MultiAssetBSModel(
        const Handle<YieldTermStructure>&                                               termStructure,
        const std::vector<std::string>&                                                 aliases,
        const std::vector<ext::shared_ptr<QuantLib::GeneralizedBlackScholesProcess>>&   processes,
        const RealStochasticProcess::MatA&                                              correlations)
    : termStructure_(termStructure), processes_(processes) {
        QL_REQUIRE(processes_.size() > 0, "No BS processes supplied");
        QL_REQUIRE(processes_.size() == aliases.size(), "Number of processes doesn't match aliases");
        for (size_t k = 0; k < aliases.size(); ++k) index_[aliases[k]] = k;
        QL_REQUIRE(processes_.size() == correlations.size(), "Number of processes doesn't match correlation");
        for (size_t k=0; k< correlations.size(); ++k)
            QL_REQUIRE(processes_.size() == correlations[k].size(), "Number of processes doesn't match correlation");
        DT_ = TemplateAuxilliaries::cholesky(correlations);
        //check whether it is a diagonal matrix
        bool isDiagonal = true;
        for (size_t k = 0; k < correlations.size(); ++k) {
            for (size_t l = k + 1; l < correlations.size(); ++l) {
                if (correlations[k][l] != 0) isDiagonal = false;
            }
        }		
        DT_ = RealStochasticProcess::MatA(processes.size());
        for (size_t k = 0; k<DT_.size(); ++k) DT_[k].resize(processes.size());

        for (size_t i = 0; i < processes.size(); i++)
        {
            for (size_t j = i; j < processes.size(); j++)
            {
                DT_[i][j] = correlations[i][j];
                DT_[j][i] = correlations[i][j];
            }
        }
        
        if (! isDiagonal) {
            TemplateAuxilliaries::performCholesky(DT_, DT_.size(),true);
            //DT_ = TemplateAuxilliaries::svdSqrt(correlations);
        }
        else {
            //due to ones on diagonal simply copy the matrix.
        }
    }

    MultiAssetBSModel::MultiAssetBSModel(
        const Handle<YieldTermStructure>&                                               termStructure,
        const std::vector<std::string>&                                                 aliases,
        const std::vector<ext::shared_ptr<QuantLib::GeneralizedBlackScholesProcess>>&   processes)
        : termStructure_(termStructure), processes_(processes) {
        QL_REQUIRE(processes_.size() > 0, "No BS processes supplied");
        //no correlation matrix means, we simply assume independence
        RealStochasticProcess::MatA corrM = RealStochasticProcess::MatA(processes.size());
        for (size_t k = 0; k<corrM.size(); ++k) corrM[k].resize(processes.size());

        for (size_t i = 0; i < processes.size(); i++)
        {
            for (size_t j = 0; j < processes.size(); j++)
            {
                if (i == j) {
                    corrM[i][j] = 1;
                }
                else corrM[i][j] = 0;
            }
        }

        DT_ = RealStochasticProcess::MatA(MultiAssetBSModel(termStructure, aliases, processes, corrM).DT_);
        for (size_t k = 0; k < aliases.size(); ++k) index_[aliases[k]] = k; // not transferable from other constructor.
    }

    MultiAssetBSModel::MultiAssetBSModel(const Handle<YieldTermStructure>&              termStructure,
        const std::vector<std::string>&                                                 aliases,
        const std::vector<ext::shared_ptr<QuantLib::LocalVolSurface>>&                  localVolSurfaces,
        const RealStochasticProcess::MatA&                                              correlations)
        : termStructure_(termStructure), localVolSurfaces_(localVolSurfaces){
    
        initProcessesFromSurface();
        DT_ = RealStochasticProcess::MatA(MultiAssetBSModel(termStructure, aliases, processes_, correlations).DT_);
        for (size_t k = 0; k < aliases.size(); ++k) index_[aliases[k]] = k; 
    }

    MultiAssetBSModel::MultiAssetBSModel(const Handle<YieldTermStructure>&              termStructure,
        const std::vector<std::string>&                                                 aliases,
        const std::vector<ext::shared_ptr<QuantLib::LocalVolSurface>>&				    localVolSurfaces)
    {
        RealStochasticProcess::MatA corrM = RealStochasticProcess::MatA(localVolSurfaces.size());
        for (size_t k = 0; k<corrM.size(); ++k) corrM[k].resize(localVolSurfaces.size());

        for (size_t i = 0; i < localVolSurfaces.size(); i++)
        {
            for (size_t j = 0; j < localVolSurfaces.size(); j++)
            {
                if (i == j) {
                    corrM[i][j] = 1;
                }
                else corrM[i][j] = 0;
            }
        }
        *this = MultiAssetBSModel(termStructure, aliases, localVolSurfaces, corrM);
    }


    // initial values for simulation
    RealStochasticProcess::VecP MultiAssetBSModel::initialValues() {
        return RealStochasticProcess::VecP(size(), 0.0);
    }
    // a[t,X(t)]
    RealStochasticProcess::VecA MultiAssetBSModel::drift(const QuantLib::Time t, const VecA& X) {
        RealStochasticProcess::VecA nu(processes_.size());
        // todo: make sure all processes use same domestic/riskfree rate...
        for (Size k = 0; k < processes_.size(); ++k) {
            QuantLib::Real S = processes_[k]->x0() * std::exp(X[k]);
            nu[k] = processes_[k]->drift(t, S);
        }
        return nu;
    }
    // b[t,X(t)]
    RealStochasticProcess::MatA MultiAssetBSModel::diffusion(const QuantLib::Time t, const VecA& X) {
        RealStochasticProcess::MatA b(DT_);
        for (Size i = 0; i < size(); ++i) {
            QuantLib::Real S = processes_[i]->x0() * std::exp(X[i]);
            QuantLib::Real sigma = processes_[i]->diffusion(t, S);
            for (Size j = 0; j < factors(); ++j) {
                b[i][j] *= sigma;
            }
        }
        return b;
    }

    void MultiAssetBSModel::evolve(const QuantLib::Time t0, const VecA& X0, const QuantLib::Time dt, const VecD& dW, VecA& X1) {
        // we approximate the local vol diffusion by a Black-Scholes process on [t, t+dt]
        for (Size i = 0; i < X1.size(); ++i) {
            X1[i] = 0.0;
            for (Size j = 0; j < dW.size(); ++j) X1[i] += DT_[i][j] * dW[j];
            // sigma represents the average volatility on [t, t+dt]
            // here we use a first very simple approximation
            Real S = processes_[i]->x0() * std::exp(X0[i]);
            Real sigma;
            if (localVolSurfaces_.size() > 0)
                //localVolSurfaces_ may be an InterpolatedLocalVolSurface which has better performance.
                sigma = localVolSurfaces_[i]->localVol(t0,S,true);
            else
                sigma= processes_[i]->diffusion(t0, S);
            // We may integrate the drift exactly (given the approximate volatility)
            Real B_d = processes_[i]->riskFreeRate()->discount(t0) / processes_[i]->riskFreeRate()->discount(t0 + dt);
            Real B_f = processes_[i]->dividendYield()->discount(t0) / processes_[i]->dividendYield()->discount(t0 + dt);
            X1[i] = X0[i] + std::log(B_d / B_f) - 0.5 * sigma * sigma * dt + sigma * X1[i] * std::sqrt(dt);
        }
    }


    void MultiAssetBSModel::initProcessesFromSurface() {
        if (localVolSurfaces_.size() > 0) {
            processes_.resize(localVolSurfaces_.size());
            for (size_t i = 0; i < localVolSurfaces_.size(); i++)
            {
                processes_[i] = ext::shared_ptr<QuantLib::GeneralizedBlackScholesProcess>(new GeneralizedBlackScholesProcess(
                    localVolSurfaces_[i]->getUnderlying(),
                    localVolSurfaces_[i]->getDividendTS(),
                    localVolSurfaces_[i]->getInterestRateTS(),
                    localVolSurfaces_[i]->getBlackSurface()
                ));
            }
        }
    }
}



