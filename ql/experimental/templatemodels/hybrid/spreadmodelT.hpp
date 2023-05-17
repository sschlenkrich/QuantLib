/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2019 Sebastian Schlenkrich

*/

/*! \file spreadmodelT.hpp
    \brief (MC) pricing for multi-asset model with stochastic interest rates and stochastic credit spread
               
*/


#ifndef quantlib_templatespreadmodel_hpp
#define quantlib_templatespreadmodel_hpp

#include <ql/errors.hpp>

#include <ql/experimental/templatemodels/stochasticprocessT.hpp>
#include <ql/experimental/templatemodels/hybrid/assetmodelT.hpp>

//#include <ql/experimental/templatemodels/qgaussian2/quasigaussianmodel2T.hpp>

//#include <ql/experimental/templatemodels/auxilliaries/qrfactorisationT.hpp>
//#include <ql/experimental/templatemodels/auxilliaries/choleskyfactorisationT.hpp>


namespace QuantLib {

    // Declaration of the asset model class
    template <class DateType, class PassiveType, class ActiveType>
    class SpreadModelT : public StochasticProcessT<DateType, PassiveType, ActiveType> {
    protected:

        // container class definitions
        typedef std::vector<DateType>                      VecD;
        typedef std::vector<PassiveType>                   VecP; 
        typedef std::vector<ActiveType>                    VecA;
        typedef std::vector< std::vector<DateType> >       MatD;
        typedef std::vector< std::vector<PassiveType> >    MatP;
        typedef std::vector< std::vector<ActiveType> >     MatA;

        typedef StochasticProcessT<DateType, PassiveType, ActiveType> ModelType;

        // class members
        ext::shared_ptr<ModelType>  baseModel_;
        ext::shared_ptr<ModelType>  sprdModel_;
        MatA                        correlations_;
        MatA                        L_;   // factorised correlation matrix

        // [ base_state, credit_state ]
        size_t size_;     // the number of state variables
        size_t factors_;  // the number of random factors

    public:  
        SpreadModelT(
            const ext::shared_ptr<ModelType>    baseModel,
            const ext::shared_ptr<ModelType>    sprdModel,
            const MatA&                         correlations)
            : baseModel_(baseModel), sprdModel_(sprdModel), correlations_(correlations), L_(0) {
            // we perform a couple of consistency checks
            QL_REQUIRE(baseModel_, "SpreadModel: Base model required.");
            QL_REQUIRE(sprdModel_, "SpreadModel: Credit model required.");
            // manage model indices in state variable
            size_ = baseModel_->size() + sprdModel_->size();
            factors_ = baseModel_->factors() + sprdModel_->factors();
            // we need to check and factor the hybrid correlation matrix
            if (correlations_.size() > 0) {  // an empty correlation matrix represents the identity matrix
                QL_REQUIRE(correlations_.size() == factors(), "SpreadModel: correlations_.size()==factors() required.");
                for (size_t k = 0; k < correlations_.size(); ++k) {
                    QL_REQUIRE(correlations_[k].size() == factors(), "SpreadModel: correlations_[k].size()==factors() required.");
                    QL_REQUIRE(correlations_[k][k] == 1.0, "SpreadModel: correlations_[k][K] == 1.0 required.");
                }
                for (size_t k = 0; k < correlations_.size(); ++k) {
                    for (size_t j = k + 1; j < correlations_.size(); ++j)
                        QL_REQUIRE(correlations_[k][j] == correlations_[j][k], "SpreadModel: correlations_[k][j] == correlations_[j][k] required.");
                }
                // We assume that base model correlation is identity.
                // User must ensure that correlation is not applied multiple times to input Brownians.
                // Rates model correlations are incorporated via individual rates models.
                L_ = TemplateAuxilliaries::cholesky(correlations_);
            }
        }

        // inspectors
        const ext::shared_ptr<ModelType>& baseModel() { return baseModel_;      }
        const ext::shared_ptr<ModelType>& sprdModel() { return sprdModel_;      }
        const MatA& correlations()                    { return correlations_;   }
        const MatA& L()                               { return L_;              }


        // stochastic process interface
        // dimension of combined state variable
        inline virtual size_t size()    { return size_; }
        // stochastic factors
        inline virtual size_t factors() { return factors_; }
        // initial values for simulation
        inline virtual VecP initialValues() {
            VecP X(size(),0.0);
            // base model
            VecP x(baseModel_->initialValues());
            for (size_t i = 0; i < baseModel_->size(); ++i) X[i] = x[i];
            VecP y(sprdModel_->initialValues());
            for (size_t i = 0; i < sprdModel_->size(); ++i) X[baseModel_->size() + i] = y[i];
            return X;
        }

        // a[t,X(t)]
        inline virtual VecA drift( const DateType t, const VecA& X) {
            VecA a(size(),0.0);
            QL_FAIL("SpreadModel: drift not implemented. Use evolve.");
            // finished
            return a;
        }

        // b[t,X(t)]
        inline virtual MatA diffusion( const DateType t, const VecA& X) {
            MatA b(size(), VecA(factors(),0.0));
            QL_FAIL("SpreadModel: diffusion not implemented. Use evolve.");
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
            // evolve base model
            {
                VecD dw(dZ.begin(), dZ.begin() + baseModel_->factors());
                VecA x0(X0.begin(), X0.begin() + baseModel_->size());
                VecA x1(baseModel_->size(), 0.0);
                baseModel_->evolve(t0, x0, dt, dw, x1);
                for (size_t i = 0; i < baseModel_->size(); ++i) X1[i] = x1[i];
            }
            // evolve credit model
            {
                VecD dw(dZ.begin() + baseModel_->factors(), dZ.begin() + baseModel_->factors() + sprdModel_->factors());
                VecA x0(X0.begin() + baseModel_->size(), X0.begin() + baseModel_->size() + sprdModel_->size());
                VecA x1(sprdModel_->size(), 0.0);
                sprdModel_->evolve(t0, x0, dt, dw, x1);
                for (size_t i = 0; i < sprdModel_->size(); ++i) X1[baseModel_->size() + i] = x1[i];
            }
        }

        // the numeraire is credit-risky bank account
        inline virtual ActiveType numeraire(const DateType t, const VecA& X) {
            VecA x(X.begin(), X.begin() + baseModel_->size());
            VecA y(X.begin() + baseModel_->size(), X.begin() + baseModel_->size() + sprdModel_->size());
            return baseModel_->numeraire(t, x) * sprdModel_->numeraire(t, y);
        }

        // asset calculation is delegated to base model
        inline virtual ActiveType asset(const DateType t, const VecA& X, const std::string& alias) {
            VecA x(X.begin(), X.begin() + baseModel_->size());
            return baseModel_->asset(t, x, alias);
        }

        // a credit-risky zero coupon bond
        inline virtual ActiveType zeroBond(const DateType t, const DateType T, const VecA& X) { 
            VecA x(X.begin(), X.begin() + baseModel_->size());
            VecA y(X.begin() + baseModel_->size(), X.begin() + baseModel_->size() + sprdModel_->size());
            return baseModel_->zeroBond(t, T, x) * sprdModel_->zeroBond(t, T, y);
        }

        // a foreign or domestic currency zero coupon bond is delegated to base model
        inline virtual ActiveType zeroBond(const DateType t, const DateType T, const VecA& X, const std::string& alias) { 
            VecA x(X.begin(), X.begin() + baseModel_->size());
            return baseModel_->zeroBond(t, T, x, alias);
        }

        // the short rate over an integration time period
        // this is required for drift calculation in multi-asset and hybrid models
        inline virtual ActiveType shortRate(const DateType t0, const DateType dt, const VecA& X0, const VecA& X1) {
            VecA x0(X0.begin(), X0.begin() + baseModel_->size());
            VecA x1(X1.begin(), X1.begin() + baseModel_->size());
            VecA y0(X0.begin() + baseModel_->size(), X0.begin() + baseModel_->size() + sprdModel_->size());
            VecA y1(X1.begin() + baseModel_->size(), X1.begin() + baseModel_->size() + sprdModel_->size());
            return baseModel_->shortRate(t0, dt, x0, x1) + sprdModel_->shortRate(t0, dt, y0, y1);
        }

        virtual std::vector< std::string > stateAliases() {
            std::vector< std::string > aliases(size());
            std::vector< std::string > baseAliases = baseModel_->stateAliases();
            for (size_t i = 0; i < baseAliases.size(); ++i) aliases[i] = "base_" + baseAliases[i];
            std::vector< std::string > spreadAliases = sprdModel_->stateAliases();
            for (size_t i = 0; i < spreadAliases.size(); ++i) aliases[baseModel_->size() + i] = "sprd_" + spreadAliases[i];
            return aliases;
        }

        virtual std::vector< std::string > factorAliases() {
            std::vector< std::string > aliases(factors());
            std::vector< std::string > baseAliases = baseModel_->factorAliases();
            for (size_t i = 0; i < baseAliases.size(); ++i) aliases[i] = "base_" + baseAliases[i];
            std::vector< std::string > spreadAliases = sprdModel_->factorAliases();
            for (size_t i = 0; i < spreadAliases.size(); ++i) aliases[baseModel_->factors() + i] = "sprd_" + spreadAliases[i];
            return aliases;
        }


    };

}

#endif  /* ifndef quantlib_templatespreadmodel_hpp */
