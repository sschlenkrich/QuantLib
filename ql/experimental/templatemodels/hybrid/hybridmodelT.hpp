/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2019 Sebastian Schlenkrich

*/

/*! \file hybridmodelT.hpp
    \brief (MC) pricing for multi-asset model with stochastic interest rates
               
*/


#ifndef quantlib_templatehybridmodel_hpp
#define quantlib_templatehybridmodel_hpp

#include <map>

#include <ql/errors.hpp>

#include <ql/experimental/templatemodels/stochasticprocessT.hpp>
#include <ql/experimental/templatemodels/hybrid/assetmodelT.hpp>
#include <ql/experimental/templatemodels/qgaussian2/quasigaussianmodel2T.hpp>

//#include <ql/experimental/templatemodels/auxilliaries/qrfactorisationT.hpp>
//#include <ql/experimental/templatemodels/auxilliaries/choleskyfactorisationT.hpp>


namespace QuantLib {

    // Declaration of the asset model class
    template <class DateType, class PassiveType, class ActiveType>
    class HybridModelT : public StochasticProcessT<DateType, PassiveType, ActiveType> {
    protected:

        // container class definitions
        typedef std::vector<DateType>                      VecD;
        typedef std::vector<PassiveType>                   VecP; 
        typedef std::vector<ActiveType>                    VecA;
        typedef std::vector< std::vector<DateType> >       MatD;
        typedef std::vector< std::vector<PassiveType> >    MatP;
        typedef std::vector< std::vector<ActiveType> >     MatA;

        typedef AssetModelT<DateType, PassiveType, ActiveType>          AssetModelType;
        typedef StochasticProcessT<DateType, PassiveType, ActiveType>   RatesModelType;

        // class members
        std::string                      domAlias_;
        ext::shared_ptr<RatesModelType>  domRatesModel_;
        std::vector<std::string>         forAliases_;
        std::vector< ext::shared_ptr<AssetModelType> >  forAssetModels_;
        std::vector< ext::shared_ptr<RatesModelType> >  forRatesModels_;
        MatA                             correlations_;
        MatA                             L_;   // factorised correlation matrix

        size_t size_;     // the number of state variables
        size_t factors_;  // the number of random factors

        // keep the first index in state Y for the foreign (!) models
        // [ asset1, rates1, asset2, rates2, ...]
        // domestic rates model state starts at Y[0] and is not stored explicitely
        std::vector<size_t>  modelsStartIdx_;

        // aliases are keys
        std::map<std::string, size_t>  index_;

        // we incorporate hybrid volatility adjustments
        VecD  hybAdjTimes_;   // [ 0=T0, T1, ..., Tn ]
        MatA  localVol_;      // [ [ v0, ..., vn ] for k in forAliases_ ]
        MatA  hybrdVol_;      // [ [ v0, ..., vn ] for k in forAliases_ ]
        MatA  hybVolAdj_;     // [ [ v0, ..., vn ] for k in forAliases_ ]

    public:  
        HybridModelT(const std::string                        domAlias,
                     const ext::shared_ptr<RatesModelType>    domRatesModel,
                     const std::vector<std::string>&          forAliases,
                     const std::vector< ext::shared_ptr<AssetModelType> >& forAssetModels,
                     const std::vector< ext::shared_ptr<RatesModelType> >& forRatesModels,
                     const MatA&                                           correlations,
                     const VecD&                              hybAdjTimes = VecD(),
                     const MatA&                              hybVolAdj   = MatA()
        )
            : domAlias_(domAlias), domRatesModel_(domRatesModel), forAliases_(forAliases),
              forAssetModels_(forAssetModels), forRatesModels_(forRatesModels), correlations_(correlations),
              L_(0), modelsStartIdx_(2*forAliases.size()),
              hybAdjTimes_(hybAdjTimes), hybVolAdj_(hybVolAdj), localVol_(0), hybrdVol_(0) {
            // we perform a couple of consistency checks
            QL_REQUIRE(domAlias_.compare("") != 0, "HybridModel: Domestic alias required.");
            QL_REQUIRE(domRatesModel_, "HybridModel: Domestic rates model required.");
            // QL_REQUIRE(forAliases_.size()>0, "HybridModel: forAliases.size()>0 required.");
            QL_REQUIRE(forAliases_.size()==forAssetModels_.size(), "HybridModel: forAliases.size()==forAssetModels.size() required.");
            QL_REQUIRE(forAliases_.size()==forRatesModels_.size(), "HybridModel: forAliases.size()==forRatesModels.size() required.");
            for (size_t k = 0; k < forAliases_.size(); ++k) index_[forAliases_[k]] = k;
            // we also check the individual components to avoid surprises at later stages...
            for (size_t k = 0; k < forAliases.size(); ++k) {
                QL_REQUIRE(forAliases_[k].compare("") != 0, "HybridModel: Foreign alias required.");
                QL_REQUIRE(forAssetModels_[k], "HybridModel: Foreign asset model required.");
                QL_REQUIRE(forRatesModels_[k], "HybridModel: Foreign rates model required.");
            }
            // adjuster times be consistent
            if (hybAdjTimes_.size()>0) QL_REQUIRE(hybAdjTimes_[0]==0.0, "HybridModel: hybAdjTimes_[0]==0.0 required.");
            for (size_t k=1; k<hybAdjTimes_.size(); ++k)
                QL_REQUIRE(hybAdjTimes_[k]>hybAdjTimes_[k-1], "HybridModel: hybAdjTimes_[k]>hybAdjTimes_[k-1] required.");
            // adjuster values be consistent
            if (hybAdjTimes_.size()>0) {
                QL_REQUIRE(hybVolAdj_.size()==forAliases_.size(), "HybridModel: hybVolAdj_.size()==forAliases_.size() required.");
                for (size_t k=0; k<hybVolAdj_.size(); ++k)
                    QL_REQUIRE(hybAdjTimes_.size()==hybVolAdj_[k].size(), "HybridModel: hybAdjTimes_.size()==hybVolAdj_[k].size() required.");
            }
            // manage model indices in state variable
            size_ = domRatesModel_->size();
            factors_ = domRatesModel_->factors();
            size_t lastModelIdx = domRatesModel_->size();
            for (size_t k = 0; k < forAliases.size(); ++k) {
                size_ += (forAssetModels_[k]->size() + forRatesModels_[k]->size());
                factors_ += (forAssetModels_[k]->factors() + forRatesModels_[k]->factors());
                modelsStartIdx_[2*k]   = lastModelIdx;					
                modelsStartIdx_[2*k+1] = lastModelIdx + forAssetModels_[k]->size();
                lastModelIdx += (forAssetModels_[k]->size() + forRatesModels_[k]->size());
            }
            // we need to check and factor the hybrid correlation matrix
            if (correlations_.size() > 0) {  // an empty correlation matrix represents the identity matrix
                QL_REQUIRE(correlations_.size() == factors(), "HybridModel: correlations_.size()==factors() required.");
                for (size_t k = 0; k < correlations_.size(); ++k) {
                    QL_REQUIRE(correlations_[k].size() == factors(), "HybridModel: correlations_[k].size()==factors() required.");
                    QL_REQUIRE(correlations_[k][k] == 1.0, "HybridModel: correlations_[k][K] == 1.0 required.");
                }
                for (size_t k = 0; k < correlations_.size(); ++k) {
                    for (size_t j = k + 1; j < correlations_.size(); ++j)
                        QL_REQUIRE(correlations_[k][j] == correlations_[j][k], "HybridModel: correlations_[k][j] == correlations_[j][k] required.");
                }
                // We should also ensure that rates model correlations are identity here.
                // Rates model correlations are incorporated via individual rates models.
                L_ = TemplateAuxilliaries::cholesky(correlations_);
            }
        }

        // inspectors
        const std::string domAlias()                           { return domAlias_;       }
        const ext::shared_ptr<RatesModelType>& domRatesModel() { return domRatesModel_;  }
        const std::vector<std::string>&  forAliases()          { return forAliases_;     }
        const std::vector< ext::shared_ptr<AssetModelType> >& forAssetModels() { return forAssetModels_; }
        const std::vector< ext::shared_ptr<RatesModelType> >& forRatesModels() { return forRatesModels_; }
        const MatA& correlations()                             { return correlations_;   }
        const MatA& L()                                        { return L_;              }
        inline const std::vector<size_t>& modelsStartIdx()     { return modelsStartIdx_; }
        const VecD& hybAdjTimes() { return hybAdjTimes_; }
        const MatA& localVol()    { return localVol_; }
        const MatA& hybrdVol()    { return hybrdVol_; }
        const MatA& hybVolAdj()   { return hybVolAdj_; }

        // stochastic process interface
        // dimension of combined state variable
        inline virtual size_t size()    { return size_; }
        // stochastic factors
        inline virtual size_t factors() { return factors_; }
        // initial values for simulation
        inline virtual VecP initialValues() {
            VecP X(size(),0.0);
            // domestic model
            VecP x(domRatesModel_->initialValues());
            for (size_t i = 0; i < domRatesModel_->size(); ++i) X[i] = x[i];
            for (size_t k = 0; k < forAliases_.size(); ++k) {
                x = forAssetModels_[k]->initialValues();
                for (size_t i = 0; i < forAssetModels_[k]->size(); ++i) X[modelsStartIdx_[2 * k] + i] = x[i];
                x = forRatesModels_[k]->initialValues();
                for (size_t i = 0; i < forRatesModels_[k]->size(); ++i) X[modelsStartIdx_[2 * k + 1] + i] = x[i];
            }
            return X;
        }

        // a[t,X(t)]
        inline virtual VecA drift( const DateType t, const VecA& X) {
            VecA a(size(),0.0);
            QL_FAIL("HybridModel: drift not implemented. Use evolve.");
            // finished
            return a;
        }

        // b[t,X(t)]
        inline virtual MatA diffusion( const DateType t, const VecA& X) {
            MatA b(size(), VecA(factors(),0.0));
            QL_FAIL("HybridModel: diffusion not implemented. Use evolve.");
            // finished
            return b;
        }

        // simple Euler step
        inline virtual void evolve(const DateType t0, const VecA& X0, const DateType dt, const VecD& dW, VecA& X1) {
            // correlate Brownian increments
            VecD dZ(dW.size(),0.0);
            if (L_.size() > 0) {
                for (size_t i = 0; i < dW.size(); ++i) // exploit lower triangular structure
                    for (size_t j = 0; j <= i; ++j) dZ[i] += L_[i][j] * dW[j];
            }
            else {
                dZ = dW;  // we would rather avoid copying here...
            }
            if (forAliases_.size() == 0) { // we may take a short cut
                domRatesModel_->evolve(t0, X0, dt, dZ, X1);
                return;
            }
            // evolve domestic rates
            VecD dw(dZ.begin(), dZ.begin() + domRatesModel_->factors());
            VecA x0(X0.begin(), X0.begin() + domRatesModel_->size());
            VecA x1(domRatesModel_->size(),0.0);
            domRatesModel_->evolve(t0, x0, dt, dw, x1);
            for (size_t i = 0; i < domRatesModel_->size(); ++i) X1[i] = x1[i];
            // we need the domestic drift r_d
            ActiveType r_d = domRatesModel_->shortRate(t0, dt, x0, x1);
            // foreign models...
            size_t corrStartIdx = domRatesModel_->factors();
            for (size_t k = 0; k < forAliases_.size(); ++k) {
                // carefully collect Brownian increments
                VecD dw_asset(dZ.begin() + corrStartIdx, dZ.begin() + corrStartIdx + forAssetModels_[k]->factors());
                VecD dw_rates(dZ.begin() + corrStartIdx + forAssetModels_[k]->factors(),
                              dZ.begin() + corrStartIdx + forAssetModels_[k]->factors() + forRatesModels_[k]->factors());
                // we need the starting point states for evolution
                VecA y0(X0.begin() + modelsStartIdx_[2 * k], X0.begin() + modelsStartIdx_[2 * k] + forAssetModels_[k]->size());
                VecA x0(X0.begin() + modelsStartIdx_[2 * k + 1], X0.begin() + modelsStartIdx_[2 * k + 1] + forRatesModels_[k]->size());
                // Quanto adjustment
                // NOTE: this should only apply to random factors for x (not z); this has to be ensured by user via correlation
                // implementation ASSUMES correlation between FX and stoch vol is zero
                // we use the model-independent implementation to allow for credit hybrid components 
                VecA qAdj(correlations_[corrStartIdx].begin() + corrStartIdx + forAssetModels_[k]->factors(),
                          correlations_[corrStartIdx].begin() + corrStartIdx + forAssetModels_[k]->factors() + forRatesModels_[k]->factors());  // keep in mind: last random factor in qG model is for dz!
                // we need to extend the input state for our asset mode to account for drift and adjuster
                y0.resize(y0.size() + 2, 0.0);  // NOTE: this uses knowledge of the structure of asset model mc state!
                // y0[y0.size() - 2] = r_d - r_f ... is set below
                y0[y0.size() - 1] = hybridVolAdjuster(k, t0);  // hybrid vol adjuster, to be implemented...
                ActiveType assetVol = forAssetModels_[k]->volatility(t0, y0);
                for (size_t i = 0; i < qAdj.size(); ++i) qAdj[i] *= (assetVol*std::sqrt(dt));
                for (size_t i = 0; i < qAdj.size(); ++i) dw_rates[i] -= qAdj[i];  // only adjust d factors for x!
                // evolve foreign rates
                VecA x1(forRatesModels_[k]->size(),0.0);
                forRatesModels_[k]->evolve(t0, x0, dt, dw_rates, x1);
                // calculate FX drift, volAdjuster and extend input state
                ActiveType r_f = forRatesModels_[k]->shortRate(t0, dt, x0, x1);
                y0[y0.size()-2] = r_d - r_f;  // FX drift
                // finally we can evolve FX
                VecA y1(forAssetModels_[k]->size(),0.0);
                forAssetModels_[k]->evolve(t0, y0, dt, dw_asset, y1);
                // almost done; we copy MC states back in hybrid MC state
                for (size_t i = 0; i < forAssetModels_[k]->size(); ++i) X1[modelsStartIdx_[2 * k] + i]     = y1[i];
                for (size_t i = 0; i < forRatesModels_[k]->size(); ++i) X1[modelsStartIdx_[2 * k + 1] + i] = x1[i];
                // done.

                // update stochastic factor index
                corrStartIdx += (forAssetModels_[k]->factors() + forRatesModels_[k]->factors());
            }
            //QL_FAIL("HybridModel: evolve not yet implemented.");
        }

        // the numeraire in the domestic currency used for discounting future payoffs
        inline virtual ActiveType numeraire(const DateType t, const VecA& X) {
            VecA x(X.begin(), X.begin() + domRatesModel_->size());
            return domRatesModel_->numeraire(t, x); 
        }

        // asset calculation is the main purpose of this model
        inline virtual ActiveType asset(const DateType t, const VecA& Y, const std::string& alias) {
            size_t k = index_.at(alias);  // this should throw an exception if alias is unknown
            VecA y(Y.begin() + modelsStartIdx_[2 * k], Y.begin() + modelsStartIdx_[2 * k] + forAssetModels_[k]->size());
            return forAssetModels_[k]->asset(t, y, alias);
        }

        // a domestic currency zero coupon bond
        inline virtual ActiveType zeroBond(const DateType t, const DateType T, const VecA& X) { 
            VecA x(X.begin(), X.begin() + domRatesModel_->size());
            return domRatesModel_->zeroBond(t, T, x);
        }

        // a foreign currency zero coupon bond
        inline virtual ActiveType zeroBond(const DateType t, const DateType T, const VecA& X, const std::string& alias) { 
            if (alias.compare(domAlias_) == 0) return zeroBond(t, T, X);
            size_t k = index_.at(alias);  // this should throw an exception if alias is unknown
            VecA x(X.begin() + modelsStartIdx_[2 * k + 1], X.begin() + modelsStartIdx_[2 * k + 1] + forRatesModels_[k]->size());
            return forRatesModels_[k]->zeroBond(t, T, x);
        }

        virtual std::vector< std::string > stateAliases() {
            std::vector< std::string > aliases(size());
            std::vector< std::string > domAliases = domRatesModel_->stateAliases();
            for (size_t i = 0; i < domAliases.size(); ++i) aliases[i] = domAlias_ + "_" + domAliases[i];
            size_t startIdx = domRatesModel_->size();
            for (size_t k = 0; k < forAliases_.size(); ++k) {
                std::vector< std::string > forAssetAliases = forAssetModels_[k]->stateAliases();
                for (size_t i = 0; i < forAssetAliases.size(); ++i) aliases[startIdx + i] = forAliases_[k] + "_" + forAssetAliases[i];
                startIdx += forAssetModels_[k]->size();
                std::vector< std::string > forRatesAliases = forRatesModels_[k]->stateAliases();
                for (size_t i = 0; i < forRatesAliases.size(); ++i) aliases[startIdx + i] = forAliases_[k] + "_" + forRatesAliases[i];
                startIdx += forRatesModels_[k]->size();
            }
            return aliases;
        }

        virtual std::vector< std::string > factorAliases() {
            std::vector< std::string > aliases(factors());
            std::vector< std::string > domAliases = domRatesModel_->factorAliases();
            for (size_t i = 0; i < domAliases.size(); ++i) aliases[i] = domAlias_ + "_" + domAliases[i];
            size_t startIdx = domRatesModel_->factors();
            for (size_t k = 0; k < forAliases_.size(); ++k) {
                std::vector< std::string > forAssetAliases = forAssetModels_[k]->factorAliases();
                for (size_t i = 0; i < forAssetAliases.size(); ++i) aliases[startIdx + i] = forAliases_[k] + "_" + forAssetAliases[i];
                startIdx += forAssetModels_[k]->factors();
                std::vector< std::string > forRatesAliases = forRatesModels_[k]->factorAliases();
                for (size_t i = 0; i < forRatesAliases.size(); ++i) aliases[startIdx + i] = forAliases_[k] + "_" + forRatesAliases[i];
                startIdx += forRatesModels_[k]->factors();
            }
            return aliases;
        }

        // incorporate hybrid volatility adjuster methodology
        // this method is unsafe, we do not check index and time range
        inline ActiveType hybridVolAdjuster(size_t forIdx, DateType t) {
            if (hybAdjTimes_.size() == 0) return 0.0;  // default
            if (t >= hybAdjTimes_.back()) return hybVolAdj_[forIdx].back();  // constant extrapolation or single adjuster time
            // if we end up here we have at least two adjuster times; recall T0=0
            // find index in ascending vector, evaluate n s.t. t[n-1] < t <= t[n]
            // bisection search
            size_t a=0, b=hybAdjTimes_.size()-1;
            while (b-a>1) {
                size_t s = (a + b) / 2;
                if (t<= hybAdjTimes_[s]) b = s;
                else                     a = s;
            }
            DateType rho = (t - hybAdjTimes_[b-1]) / (hybAdjTimes_[b] - hybAdjTimes_[b-1]);
            return (1.0-rho)*hybVolAdj_[forIdx][b-1] + rho*hybVolAdj_[forIdx][b];
        }

        // this method does the actual work...
        inline void recalculateHybridVolAdjuster(const VecD& hybAdjTimes = VecD()) {
            if (hybAdjTimes.size() > 0) hybAdjTimes_ = hybAdjTimes;  // if we don't supply time grid we want to keep the current grid
            if (hybAdjTimes_.size() == 0)  // patological case, do nothing
                return;
            // we need to check for consistent times again
            QL_REQUIRE(hybAdjTimes_[0] == 0.0, "HybridModel: hybAdjTimes_[0]==0.0 required.");
            for (size_t k = 1; k < hybAdjTimes_.size(); ++k)
                QL_REQUIRE(hybAdjTimes_[k] > hybAdjTimes_[k - 1], "HybridModel: hybAdjTimes_[k]>hybAdjTimes_[k-1] required.");
            localVol_ = MatA(forAliases_.size(), VecA(hybAdjTimes_.size(), 0.0));
            hybrdVol_ = MatA(forAliases_.size(), VecA(hybAdjTimes_.size(), 1.0));  // 1.0 is required for p calculation
            hybVolAdj_ = MatA(forAliases_.size(), VecA(hybAdjTimes_.size(), 0.0));
            VecA S0(forAliases_.size());  // we need initial asset values for forwards and ATM vol
            for (size_t i = 0; i < S0.size(); ++i) {
                S0[i] = forAssetModels_[i]->asset(0.0, forAssetModels_[i]->initialValues(), "");
                VecA y0 = forAssetModels_[i]->initialValues();
                y0.resize(y0.size() + 2, 0.0);         // this is asset model-dependent
                localVol_[i][0] = forAssetModels_[i]->volatility(0.0, y0);
                hybrdVol_[i][0] = localVol_[i][0];  // initial condition
                hybVolAdj_[i][0] = hybrdVol_[i][0] - localVol_[i][0];
            }
            if (hybAdjTimes_.size() == 1)  // nothing else to do				
                return;
            // now we start with the actual methodology...
            size_t corrStartIdx = domRatesModel_->factors();
            for (size_t i = 0; i < S0.size(); ++i) {
                // we collect all relevant correlations
                VecA rhoY0X1(domRatesModel_->factors());      // ASSUME vol-FX correlation is zero
                VecA rhoX1Y1(forRatesModels_[i]->factors());  // ASSUME vol-FX correlation is zero
                MatA rhoY0Y1(domRatesModel_->factors(), VecA(forRatesModels_[i]->factors())); // ASSUME all vol-... correlation are zero
                for (size_t row = 0; row < rhoY0X1.size(); ++row)
                    rhoY0X1[row] = correlations_[row][corrStartIdx];
                for (size_t col = 0; col < rhoX1Y1.size(); ++col)
                    rhoX1Y1[col] = correlations_[corrStartIdx][corrStartIdx + forAssetModels_[i]->factors() + col];
                for (size_t row = 0; row < rhoY0X1.size(); ++row)
                    for (size_t col = 0; col < rhoX1Y1.size(); ++col)
                        rhoY0Y1[row][col] = correlations_[row][corrStartIdx + forAssetModels_[i]->factors() + col];
                // update stochastic factor index
                corrStartIdx += (forAssetModels_[i]->factors() + forRatesModels_[i]->factors());
                // 
                for (size_t k = 1; k < hybAdjTimes_.size(); ++k) {
                    // we will use current time repeatedly
                    DateType T = hybAdjTimes_[k];
                    // ATM forward and effective local volatility
                    ActiveType dfDom = domRatesModel_->zeroBond(0.0, T, domRatesModel_->initialValues());
                    ActiveType dfFor = forRatesModels_[i]->zeroBond(0.0, T, forRatesModels_[i]->initialValues());
                    ActiveType S = S0[i] * dfFor / dfDom;  // maybe it's worth to save S for debugging
                    VecA y(forAssetModels_[i]->size() + 2, 0.0);
                    y[0] = log(S / S0[i]);
                    localVol_[i][k] = forAssetModels_[i]->volatility(hybAdjTimes_[k], y);
                    VecA hPrime(k+1);
                    for (size_t j = 0; j < hPrime.size(); ++j) {
                        VecA sigmaP0      = domRatesModel_->zeroBondVolatility(hybAdjTimes_[j], T, domRatesModel_->initialValues());
                        VecA sigmaP0Prime = domRatesModel_->zeroBondVolatilityPrime(hybAdjTimes_[j], T, domRatesModel_->initialValues());
                        VecA sigmaP1      = forRatesModels_[i]->zeroBondVolatility(hybAdjTimes_[j], T, forRatesModels_[i]->initialValues());
                        VecA sigmaP1Prime = forRatesModels_[i]->zeroBondVolatilityPrime(hybAdjTimes_[j], T, forRatesModels_[i]->initialValues());

                        // we start with sigmaP0Prime * sigma0
                        VecA sigma0 = sigmaP0;
                        for (size_t l = 0; l < sigma0.size(); ++l)
                            for (size_t m = 0; m < sigmaP1.size(); ++m)
                                sigma0[l] -= rhoY0Y1[l][m] * sigmaP1[m];
                        for (size_t l = 0; l < sigma0.size(); ++l)
                            sigma0[l] += hybrdVol_[i][j] * rhoY0X1[l];  // bootstrapping enters here
                        ActiveType sum0 = 0.0;
                        for (size_t l = 0; l < sigma0.size(); ++l)
                            sum0 += sigmaP0Prime[l] * sigma0[l];

                        // now we look at sigma1 * sigmaP1Prime
                        VecA sigma1 = sigmaP1;
                        for (size_t l = 0; l < sigma1.size(); ++l)
                            for (size_t m = 0; m < sigmaP0.size(); ++m)
                                sigma1[l] -= rhoY0Y1[m][l] * sigmaP0[m];
                        for (size_t l = 0; l < sigma1.size(); ++l)
                            sigma1[l] -= hybrdVol_[i][j] * rhoX1Y1[l];  // mnus sign (!), bootstrapping enters here
                        ActiveType sum1 = 0.0;
                        for (size_t l = 0; l < sigma1.size(); ++l)
                            sum1 += sigma1[l] * sigmaP1Prime[l];

                        // collect terms and finish
                        hPrime[j] = 2.0*(sum0 + sum1);
                    }
                    ActiveType p = 0.5 * hPrime[k] * (hybAdjTimes_[k] - hybAdjTimes_[k - 1]);
                    ActiveType q = 0.5 * hPrime[k-1] * (hybAdjTimes_[k] - hybAdjTimes_[k - 1]);
                    for (size_t j=1; j<k; ++j)
                        q += 0.5 * (hPrime[j - 1] + hPrime[j]) * (hybAdjTimes_[j] - hybAdjTimes_[j - 1]);
                    // let's see if this works...
                    ActiveType root2 = p*p / 4.0 - q + localVol_[i][k] * localVol_[i][k];
                    QL_REQUIRE(root2 >= 0.0, "HybridModel: root2>=0.0 required.");
                    hybrdVol_[i][k] = -p / 2.0 + sqrt(root2);
                    QL_REQUIRE(hybrdVol_[i][k]>0.0, "HybridModel: hybrdVol_[i][k]>0.0 required.");
                    // maybe we should add some more safety checks here...
                    hybVolAdj_[i][k] = hybrdVol_[i][k] - localVol_[i][k];
                }
            }
        }
    };

}

#endif  /* ifndef quantlib_templatehybridmodel_hpp */
