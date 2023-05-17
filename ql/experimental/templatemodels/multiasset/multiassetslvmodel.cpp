/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2019 Cord Harms

*/


#include <vector>
#include <ql/errors.hpp>
#include <ql/experimental/templatemodels/auxilliaries/choleskyfactorisationT.hpp>

#include <ql/experimental/templatemodels/multiasset/multiassetslvmodel.hpp>


namespace QuantLib {

    MultiAssetSLVModel::MultiAssetSLVModel(
        const Handle<YieldTermStructure>&                                               termStructure,
        const std::vector<std::string>&                                                 aliases,
        const std::vector<ext::shared_ptr<QuantLib::HestonSLVProcess>>&				    processes,
        const RealStochasticProcess::MatA&                                              correlations)
    : termStructure_(termStructure), processes_(processes) {
        QL_REQUIRE(processes_.size() > 0, "No SLV processes supplied");
        QL_REQUIRE(processes_.size() == aliases.size(), "Number of processes doesn't match aliases");
        for (size_t k = 0; k < aliases.size(); ++k) index_[aliases[k]] = k;
        QL_REQUIRE(2*processes_.size() == correlations.size(), "Number of processes*2 doesn't match correlation");
        for (size_t k=0; k< correlations.size(); ++k)
            QL_REQUIRE(processes_.size()*2 == correlations[k].size(), "Number of processes*2 doesn't match correlation");
        //check whether it is a diagonal matrix
        bool isDiagonal = true;
        for (size_t k = 0; k < correlations.size(); ++k) {
            for (size_t l = k + 1; l < correlations.size(); ++l) {
                if (correlations[k][l] != 0) isDiagonal = false;
            }
        }
        
        
        DT_ = RealStochasticProcess::MatA(2*processes.size());
        for (size_t k = 0; k<DT_.size(); ++k) DT_[k].resize(2*processes.size());

        for (size_t i = 0; i < 2 * processes.size(); i++)
        {
            for (size_t j = i; j < 2 * processes.size(); j++)
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
    MultiAssetSLVModel::MultiAssetSLVModel(
        const Handle<YieldTermStructure>&                                               termStructure,
        const std::vector<std::string>&                                                 aliases,
        const std::vector<ext::shared_ptr<QuantLib::HestonSLVProcess>>&				    processes)
        : termStructure_(termStructure), processes_(processes) {
        QL_REQUIRE(processes_.size() > 0, "No SLV processes supplied");
        //no correlation matrix provided, consequently we simply infer asset-volvol-correlation from heston and assume
        //remaining correlations to be zero.
        RealStochasticProcess::MatA corrM = getPureHestonImpliedCorrelationMatrix();
        DT_ = RealStochasticProcess::MatA(MultiAssetSLVModel(termStructure, aliases, processes, corrM).DT_);
        for (size_t k = 0; k < aliases.size(); ++k) index_[aliases[k]] = k; // not transferable from other constructor.
    }
    // initial values for simulation
    RealStochasticProcess::VecP MultiAssetSLVModel::initialValues() {
        RealStochasticProcess::VecP init =  RealStochasticProcess::VecP(size(), 0.0);
        for (size_t i = 0; i < processes_.size(); i++)
        {
            init[i] = 0;
            init[i+ processes_.size()] = processes_[i]->v0();
        }
        return init;
    }
    // a[t,X(t)]
    RealStochasticProcess::VecA MultiAssetSLVModel::drift(const QuantLib::Time t, const VecA& X) {
        RealStochasticProcess::VecA nu(processes_.size());
        // todo: make sure all processes use same domestic/riskfree rate...
        // remember: X=[x1,x2,....,v1, v2,....]
        // HestonSLVProcess returns drift of d(lnSt)
        size_t amtProcesses = processes_.size();


        for (Size k = 0; k < processes_.size(); ++k) {
            QuantLib::Array tmp(2);
            tmp[0] = X[k];
            tmp[1] = X[k+ amtProcesses];
            tmp = processes_[k]->drift(t,tmp);
            nu[k] = tmp[0];
            nu[k+amtProcesses] = tmp[1];
        }
        return nu;
    }	
    // b[t,X(t)]
    RealStochasticProcess::MatA MultiAssetSLVModel::diffusion(const QuantLib::Time t, const VecA& X) {
        RealStochasticProcess::MatA b(DT_);
        size_t amtProcesses = processes_.size();
        for (Size i = 0; i < processes_.size(); ++i) {
            QuantLib::Array tmp(2);
            tmp[0] = X[i];
            tmp[1] = X[i + amtProcesses];
            QuantLib::Matrix sigma = processes_[i]->diffusion(t, tmp);
            for (Size j = 0; j < factors(); ++j) {
                b[i][j] *= sigma[0][0];
                b[i + amtProcesses][j] *= sigma[1][0];
                b[i + amtProcesses][j + amtProcesses] *= sigma[1][1];
            }
        }
        return b;
    }

    void MultiAssetSLVModel::evolve(const QuantLib::Time t0, const VecA& X0, const QuantLib::Time dt, const VecD& dW, VecA& X1) {
        
        size_t amtProcesses = processes_.size();
        QuantLib::Array tmp0(2);
        QuantLib::Array tmp1(2);
        QuantLib::Real s0;
        size_t amtDW;

        for (Size i = 0; i < processes_.size(); ++i) {
            
            s0 = processes_[i]->s0()->value();
            
            X1[i] = 0.0;
            X1[i + amtProcesses] = 0.0;

            amtDW = processes_[i]->isLocalVolProcess() ?  processes_.size() : dW.size();

            for (Size j = 0; j < amtDW; ++j)
            {
                X1[i] += DT_[i][j] * dW[j];
                if(!processes_[i]->isLocalVolProcess()) X1[i + amtProcesses] += DT_[i + amtProcesses][j] * dW[j];
            }
            tmp0[0] = s0*std::exp(X0[i]);
            tmp0[1] = X0[i + amtProcesses];
            if (tmp0[1] < 0) tmp0[1] = 0.00001; // may happen due to richardson extrapolation that vola becomes negative

            tmp1[0] = X1[i];
            tmp1[1] = X1[i + amtProcesses];

            //calculate an artificial random number with a correlation rho for the 1dim Heston model
            if (!processes_[i]->isLocalVolProcess()) tmp1[1] = (tmp1[1] - tmp1[0] * processes_[i]->rho()) / std::sqrt(1 - processes_[i]->rho()*processes_[i]->rho());
            
            tmp1 = processes_[i]->evolve(t0,tmp0,dt,tmp1); //evolve of HestonSLVProcess returns S instead of ln(S/S0)
            //transfer (S,v) to (log(S/S0),v)
            
            X1[i] = std::log(tmp1[0]/s0);


            if (tmp1[0] == 0 || tmp1[0] != tmp1[0]) {
                tmp1[0] = 2;
            }

            X1[i + amtProcesses] = tmp1[1];
        }
    }

    RealStochasticProcess::MatA MultiAssetSLVModel::getPureHestonImpliedCorrelationMatrix()
    {
        RealStochasticProcess::MatA corrM = RealStochasticProcess::MatA(2 * processes_.size());

        for (size_t k = 0; k<corrM.size(); ++k) corrM[k].resize(2 * processes_.size());

        for (size_t i = 0; i < 2 * processes_.size(); i++)
        {
            for (size_t j = 0; j < 2 * processes_.size(); j++)
            {
                if (i == j) {
                    corrM[i][j] = 1;
                }
                else if (i == j + processes_.size() || i + processes_.size() == j) {
                    int assetIndex = std::min(i, j);
                    corrM[i][j] = processes_[assetIndex]->rho();
                }
                else {
                    corrM[i][j] = 0;
                }
            }
        }

        return corrM;
    }
}



