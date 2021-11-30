/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2010 Sebastian Schlenkrich

*/



#ifndef quantlib_hestonmodels_hpp
#define quantlib_hestonmodels_hpp

#include <ql/types.hpp>
#include <ql/experimental/templatemodels/auxilliaries/minimADVariable2T.hpp>

#include <ql/experimental/templatemodels/stochasticprocessT.hpp>
#include <ql/experimental/templatemodels/stochvol/hestonmodelT.hpp>
#include <ql/experimental/templatemodels/stochvol/tdstochvolmodelT.hpp>
#include <ql/experimental/templatemodels/stochvol/shiftedsabrmodelT.hpp>
#include <ql/experimental/templatemodels/stochvol/quadraticlvsvmodelT.hpp>



#define _MIN_( a, b ) ( (a) < (b) ? (a) : (b) )
#define _MAX_( a, b ) ( (a) > (b) ? (a) : (b) )

namespace QuantLib {
    
    typedef HestonModelT<QuantLib::Time,QuantLib::Real,QuantLib::Real> RealHestonModel;
    typedef HestonModelT<QuantLib::Time,QuantLib::Real, MinimAD::Variable<Real> > MinimADHestonModel;

    typedef StochVolModelT<QuantLib::Time,QuantLib::Real,QuantLib::Real> RealStochVolModel;

    typedef TimeDependentStochVolModelT<QuantLib::Time,QuantLib::Real,QuantLib::Real> RealTDStochVolModel;
    typedef PWCAnalytical<QuantLib::Time,QuantLib::Real,QuantLib::Real> RealPWCStochVolModel;

    typedef ShiftedSABRModelT<QuantLib::Time,QuantLib::Real,QuantLib::Real> RealShiftedSABRModel;

    typedef QuadraticLVSVModelT<QuantLib::Time, QuantLib::Real, QuantLib::Real> RealQuadraticLVSVModel;

    class ActiveHestonModel : public RealHestonModel {
    private:
        ext::shared_ptr< MinimADHestonModel > amodel_;
    public:
        ActiveHestonModel( Real kappa,
                           Real theta,
                           Real sigma,
                           Real rho,
                           Real v0 )
        : RealHestonModel(kappa, theta, sigma, rho, v0) {
            amodel_ = ext::shared_ptr< MinimADHestonModel >( new MinimADHestonModel(kappa, theta, sigma, rho, v0) );
        }
        
        std::vector<Real> vanillaOption(
                const Real   forwardPrice,
                const Real   strikePrice,
                const Time   term,
                const int    callOrPut,
                const Real   accuracy,
                const size_t maxEvaluations) {
            std::vector<Real> results(6);
                MinimAD::Variable<Real> res = amodel_->vanillaOption(forwardPrice, strikePrice, term, callOrPut, accuracy, maxEvaluations);
                results[0] = res.value();
                results[1] = res % amodel_->kappa();
                results[2] = res % amodel_->theta();
                results[3] = res % amodel_->sigma();
                results[4] = res % amodel_->rho();
                results[5] = res % amodel_->v0();        
                return results;
        }
        

    };

}

#undef _MIN_
#undef _MAX_

#endif  /* ifndef quantlib_hestonmodels_hpp */
