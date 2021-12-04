/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2019 Cord Harms

*/



#ifndef quantlib_multiassetslvmodel_hpp
#define quantlib_multiassetslvmodel_hpp


#include <ql/experimental/templatemodels/stochasticprocessT.hpp>

#include <ql/experimental/processes/hestonslvprocess.hpp>

namespace QuantLib {


    // We model a multi-asset local volatility model by means of the normalized log-processes X_i = log[S_i/S_(0)]
    // Basis for the model, is local volatility model with integrated Heston component for stochastic volatility: cf. class HestonSLVModel.

    class MultiAssetSLVModel : public RealStochasticProcess {
    protected:
        /* information in list of HestonSLVProcesses is sufficient to imply a correlation matrix using asset-volvol-correlations.
           Caution: Remaining correlations are set to zero. Therefore, only call if no bespoke correlation matrix is provided.
        */
        RealStochasticProcess::MatA getPureHestonImpliedCorrelationMatrix();
        Handle<YieldTermStructure>                                               termStructure_;  // domestic discounting term structure
        std::map<std::string, size_t>                                            index_;
        std::vector<ext::shared_ptr<QuantLib::HestonSLVProcess>>				 processes_;
        RealStochasticProcess::MatA                                              DT_;  // D^T D = Correlations
    public:
        MultiAssetSLVModel(const Handle<YieldTermStructure>&                                              termStructure,
                          const std::vector<std::string>&                                                 aliases,
                          const std::vector<ext::shared_ptr<QuantLib::HestonSLVProcess>>&                 processes,
                          const RealStochasticProcess::MatA&                                              correlations);

        MultiAssetSLVModel(const Handle<YieldTermStructure>&                                              termStructure,
            const std::vector<std::string>&                                                               aliases,
            const std::vector<ext::shared_ptr<QuantLib::HestonSLVProcess>>&	                              processes);

        // dimension of X -> [x1,x2,....,v1, v2,....]
        virtual size_t size() { return processes_.size()*2; }
        // stochastic factors of x and z 
        virtual size_t factors() { return processes_.size()*2; } //UL and vol
        // initial values for simulation
        virtual RealStochasticProcess::VecP initialValues();
        // a[t,X(t)], not tested
        virtual RealStochasticProcess::VecA drift(const QuantLib::Time t, const VecA& X);
        // b[t,X(t)], not tested
        virtual RealStochasticProcess::MatA diffusion(const QuantLib::Time t, const VecA& X);

        virtual void evolve(const QuantLib::Time t0, const VecA& X0, const QuantLib::Time dt, const VecD& dW, VecA& X1);

        // default implementation, zero interest rates
        inline virtual QuantLib::Real numeraire(const QuantLib::Time t, const VecA& X) { return 1.0/termStructure_->discount(t); }

        // default implementation, zero interest rates
        inline virtual QuantLib::Real zeroBond(const QuantLib::Time t, const QuantLib::Time T, const VecA& X) { return termStructure_->discount(T) / termStructure_->discount(t); }

        // default implementation for single-asset models
        inline virtual QuantLib::Real asset(const QuantLib::Time t, const VecA& X, const std::string& alias) { return processes_[index_.at(alias)]->s0()->value()*std::exp(X[index_.at(alias)]); }

        // default implementation for single-asset models
        inline virtual QuantLib::Array assetAndVol(const QuantLib::Time t, const VecA& X, const std::string& alias) { 
            QuantLib::Array x0(2);
            x0[0] = processes_[index_.at(alias)]->s0()->value();
            x0[1] = processes_[index_.at(alias)]->v0();
            QuantLib::Array dx(2);
            dx[0] = X[0];
            dx[1] = X[1];
            return processes_[index_.at(alias)]->apply(x0,dx);
        }

        inline virtual QuantLib::Real forwardAsset(const QuantLib::Time t, const QuantLib::Time T, const VecA& X, const std::string& alias) {
            return asset(t, X, alias) *
                (processes_[index_.at(alias)]->dividendYield()->discount(T) / processes_[index_.at(alias)]->dividendYield()->discount(t)) /
                (processes_[index_.at(alias)]->riskFreeRate()->discount(T) / processes_[index_.at(alias)]->riskFreeRate()->discount(t));
        }

        // calculate the local volatility of the log-process of the asset
        // this is required continuous barrier estimation via Brownian Bridge
        // NOT VERIFIED THAT RETURN VALUE IS CORRECT
        inline virtual QuantLib::Real assetVolatility(const QuantLib::Time t, const VecA& X, const std::string& alias) { return processes_[index_.at(alias)]->diffusion(t, assetAndVol(t,X,alias))[0][0]; } 


    };

}



#endif  /* ifndef quantlib_multiassetslvmodel_hpp */
