/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2015 Sebastian Schlenkrich

*/



#ifndef quantlib_template2flognormalmodel_hpp
#define quantlib_template2flognormalmodel_hpp

#include <ql/experimental/templatemodels/commodity/twofactormeanreversionmodelT.hpp>



#define _MIN_( a, b ) ( (a) < (b) ? (a) : (b) )
#define _MAX_( a, b ) ( (a) > (b) ? (a) : (b) )

namespace QuantLib {

    // 2-factor mean reverting normal model
    //
    //     X(t) = phi(t) + Y(t) + Z(t), X(0)=Y(0)=0
    //    dY(t) = -a Y(t) dt  +  sigma(t) dW_Y(t)
    //    dZ(t) = -b Z(t) dt  +    eta(t) dW_Z(t)
    //    dW_Y(t) dW_Z(t) = rho dt
    //
    template <class DateType, class PassiveType, class ActiveType>
    class TwoFactorNormalModelT : public TwoFactorMeanReversionModelT<DateType,PassiveType,ActiveType> {
        // from base class
        using typename TwoFactorMeanReversionModelT<DateType,PassiveType,ActiveType>::VecD;
        using typename TwoFactorMeanReversionModelT<DateType,PassiveType,ActiveType>::VecP;
        using typename TwoFactorMeanReversionModelT<DateType,PassiveType,ActiveType>::VecA;
        using TwoFactorMeanReversionModelT<DateType,PassiveType,ActiveType>::rho_;
        using TwoFactorMeanReversionModelT<DateType,PassiveType,ActiveType>::a_;
        using TwoFactorMeanReversionModelT<DateType,PassiveType,ActiveType>::b_;


    public:
        // constructor
        
        TwoFactorNormalModelT( const Handle<IndexTermStructure>&    futureTS,
                               const VecD&                          times,
                               const VecA&                          sigma,
                               const VecA&                          eta,
                               const PassiveType                     a,
                               const PassiveType                     b,
                               const PassiveType                     rho 
                               )
        : TwoFactorMeanReversionModelT<DateType,PassiveType,ActiveType>(futureTS,times,sigma,eta,a,b,rho) {
            // check for valid parameter inputs
        }
    
        virtual ~TwoFactorNormalModelT() = default;

        // analytic formulas in base class

        // basic instruments

        // future expectation, E[ S(T) | F(t) ]
        inline virtual ActiveType futureAsset(const DateType t, const DateType T, const ActiveType Y, const ActiveType Z) {
            return phi(T) + exp(-a_*(T-t))*Y + exp(-b_*(T-t))*Z;
        }

        inline virtual ActiveType averageFuture ( const VecD& settlementTimes, const VecP& settlementWeights) {
            ActiveType fut = 0.0;
            size_t N = _MIN_(settlementTimes.size(),settlementWeights.size());
            for (size_t k=0; k<N; ++k) fut += settlementWeights[k]*phi(settlementTimes[k]);
            return fut;
        }

        inline virtual ActiveType varianceAverageFuture ( const DateType expiryTime, const VecD& settlementTimes, const VecP& settlementWeights) {
            PassiveType B = 0, C = 0;
            size_t N = _MIN_(settlementTimes.size(),settlementWeights.size());
            for (size_t k=0; k<N; ++k) {
                B += settlementWeights[k]*exp(-a_*(settlementTimes[k]-expiryTime));
                C += settlementWeights[k]*exp(-b_*(settlementTimes[k]-expiryTime));
            }
            ActiveType var = B*B*varianceY(0.0,expiryTime) + C*C*varianceZ(0.0,expiryTime) + 2.0*rho_*B*C*covarianceYZ(0.0,expiryTime);
            return var;
        }

        inline virtual ActiveType vanillaOption ( const DateType expiryTime, const VecD& settlementTimes, const VecP& settlementWeights, PassiveType strike, int callOrPut) {
            ActiveType fut = averageFuture(settlementTimes,settlementWeights);
            ActiveType var = varianceAverageFuture(expiryTime,settlementTimes,settlementWeights);
            ActiveType pv  = TemplateAuxilliaries::Bachelier(fut,strike,sqrt(var),1.0,callOrPut);
            return pv;
        }

    };

}

#undef _MIN_
#undef _MAX_

#endif  /* quantlib_template2flognormalmodel_hpp */
