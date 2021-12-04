/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2017 Sebastian Schlenkrich

*/



#ifndef quantlib_multiassetbsmodel_hpp
#define quantlib_multiassetbsmodel_hpp


#include <ql/experimental/templatemodels/stochasticprocessT.hpp>

#include <ql/processes/blackscholesprocess.hpp>
#include <ql/termstructures/volatility/equityfx/localvolsurface.hpp>

namespace QuantLib {


    // We model a multi-asset local stochastic volatility model by means of the normalized log-processes X_i = log[S_i/S_(0)]

    class MultiAssetBSModel : public RealStochasticProcess {
    protected:
        Handle<YieldTermStructure>                                               termStructure_;  // domestic discounting term structure
        std::map<std::string, size_t>                                            index_;
        std::vector<ext::shared_ptr<QuantLib::GeneralizedBlackScholesProcess>>   processes_;
        std::vector<ext::shared_ptr<QuantLib::LocalVolSurface>>				     localVolSurfaces_;
        RealStochasticProcess::MatA                                              DT_;  // D^T D = Correlations

        void initProcessesFromSurface();
    public:
        MultiAssetBSModel(const Handle<YieldTermStructure>&                                               termStructure,
                          const std::vector<std::string>&                                                 aliases,
                          const std::vector<ext::shared_ptr<QuantLib::GeneralizedBlackScholesProcess>>&   processes,
                          const RealStochasticProcess::MatA&                                              correlations);
        MultiAssetBSModel(const Handle<YieldTermStructure>&                                               termStructure,
            const std::vector<std::string>&                                                 aliases,
            const std::vector<ext::shared_ptr<QuantLib::GeneralizedBlackScholesProcess>>&   processes);
        
        //This constructor enables to directly pass a term structure, i.e. an InterpolatedLocalVolSurface
        MultiAssetBSModel(const Handle<YieldTermStructure>&                                 termStructure,
            const std::vector<std::string>&                                                 aliases,
            const std::vector<ext::shared_ptr<QuantLib::LocalVolSurface>>&				    localVolSurfaces,
            const RealStochasticProcess::MatA&                                              correlations);
        MultiAssetBSModel(const Handle<YieldTermStructure>&                                 termStructure,
            const std::vector<std::string>&                                                 aliases,
            const std::vector<ext::shared_ptr<QuantLib::LocalVolSurface>>&				    localVolSurfaces);

        // dimension of X
        inline virtual size_t size() { return processes_.size(); }
        // stochastic factors of x and z (maybe distinguish if trivially eta=0)
        inline virtual size_t factors() { return processes_.size(); }
        // initial values for simulation
        virtual RealStochasticProcess::VecP initialValues();
        // a[t,X(t)]
        virtual RealStochasticProcess::VecA drift(const QuantLib::Time t, const VecA& X);
        // b[t,X(t)]
        virtual RealStochasticProcess::MatA diffusion(const QuantLib::Time t, const VecA& X);

        virtual void evolve(const QuantLib::Time t0, const VecA& X0, const QuantLib::Time dt, const VecD& dW, VecA& X1);

        // default implementation, zero interest rates
        inline virtual QuantLib::Real numeraire(const QuantLib::Time t, const VecA& X) { return 1.0/termStructure_->discount(t); }

        // default implementation, zero interest rates
        inline virtual QuantLib::Real zeroBond(const QuantLib::Time t, const QuantLib::Time T, const VecA& X) { return termStructure_->discount(T) / termStructure_->discount(t); }

        // default implementation for single-asset models
        inline virtual QuantLib::Real asset(const QuantLib::Time t, const VecA& X, const std::string& alias) { return processes_[index_.at(alias)]->x0() * std::exp(X[index_.at(alias)]); }

        inline virtual QuantLib::Real forwardAsset(const QuantLib::Time t, const QuantLib::Time T, const VecA& X, const std::string& alias) {
            return asset(t, X, alias) *
                (processes_[index_.at(alias)]->dividendYield()->discount(T) / processes_[index_.at(alias)]->dividendYield()->discount(t)) /
                (processes_[index_.at(alias)]->riskFreeRate()->discount(T) / processes_[index_.at(alias)]->riskFreeRate()->discount(t));
        }

        // calculate the local volatility of the log-process of the asset
        // this is required continuous barrier estimation via Brownian Bridge
        inline virtual QuantLib::Real assetVolatility(const QuantLib::Time t, const VecA& X, const std::string& alias) { return processes_[index_.at(alias)]->diffusion(t, asset(t, X, alias)); }
        
    };

}



#endif  /* ifndef quantlib_multiassetbsmodel_hpp */
