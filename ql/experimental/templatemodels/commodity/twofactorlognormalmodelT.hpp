/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2015 Sebastian Schlenkrich

*/



#ifndef quantlib_template2fnormalmodel_hpp
#define quantlib_template2fnormalmodel_hpp

#include <ql/experimental/templatemodels/commodity/twofactormeanreversionmodelT.hpp>



#define _MIN_( a, b ) ( (a) < (b) ? (a) : (b) )
#define _MAX_( a, b ) ( (a) > (b) ? (a) : (b) )

namespace QuantLib {

    // 2-factor mean reverting normal model
    //
    //  ln X(t) = phi(t) + Y(t) + Z(t), X(0)=Y(0)=0
    //    dY(t) = -a Y(t) dt  +  sigma(t) dW_Y(t)
    //    dZ(t) = -b Z(t) dt  +    eta(t) dW_Z(t)
    //    dW_Y(t) dW_Z(t) = rho dt
    //
    template <class DateType, class PassiveType, class ActiveType>
    class TwoFactorLognormalModelT : public TwoFactorMeanReversionModelT<DateType,PassiveType,ActiveType> {
        // from base class
        using typename TwoFactorMeanReversionModelT<DateType,PassiveType,ActiveType>::VecD;
        using typename TwoFactorMeanReversionModelT<DateType,PassiveType,ActiveType>::VecP;
        using typename TwoFactorMeanReversionModelT<DateType,PassiveType,ActiveType>::VecA;
        using TwoFactorMeanReversionModelT<DateType,PassiveType,ActiveType>::rho_;
        using TwoFactorMeanReversionModelT<DateType,PassiveType,ActiveType>::futureTS_;
        using TwoFactorMeanReversionModelT<DateType,PassiveType,ActiveType>::a_;
        using TwoFactorMeanReversionModelT<DateType,PassiveType,ActiveType>::b_;


    public:
        // constructor
        
        TwoFactorLognormalModelT( const Handle<IndexTermStructure>& futureTS,
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

        virtual ~TwoFactorLognormalModelT() = default;

        // analytic formulas in base class

        inline virtual const ActiveType phi(const DateType t) const {
            ActiveType var = varianceY(0,t) + varianceZ(0,t) + 2.0*rho_*covarianceYZ(0,t);
            return log(futureTS_->value(t)) - var/2.0;
        }


        // basic instruments

        // future expectation
        inline virtual ActiveType futureAsset(const DateType t, const DateType T, const ActiveType Y, const ActiveType Z) {
            ActiveType mu  = phi(T) + exp(-a_*(T-t))*Y + exp(-b_*(T-t))*Z;
            ActiveType var = varianceY(t,T) + varianceZ(t,T) + 2.0*rho_*covarianceYZ(t,T);
            return exp(mu + var/2);
        }

        inline virtual ActiveType averageFuture ( const VecD& settlementTimes, const VecP& settlementWeights) {
            ActiveType fut = 0.0;
            size_t N = _MIN_(settlementTimes.size(),settlementWeights.size());
            for (size_t k=0; k<N; ++k) fut += settlementWeights[k] * futureAsset(0.0,settlementTimes[k],0.0,0.0);
            return fut;
        }

        // calculate approximate log-variance as input to Black formula below
        // various approaches are distinguished by approxType
        // NOTE: approximation quality needs to be verified
        inline virtual ActiveType averageFutureStDev ( const DateType expiryTime, const VecD& settlementTimes, const VecP& settlementWeights, int approxType = 0) {
            PassiveType B = 0, C = 0, B2 = 0, C2 = 0, BC = 0;
            PassiveType weight=0.0;
            size_t N = _MIN_(settlementTimes.size(),settlementWeights.size());
            if (approxType==0) { // average log-variance
                ActiveType  var=0.0;
                for (size_t k=0; k<N; ++k) {
                    B = exp(-a_*(settlementTimes[k]-expiryTime));
                    C = exp(-b_*(settlementTimes[k]-expiryTime));
                    B2 += settlementWeights[k]*B*B;
                    C2 += settlementWeights[k]*C*C;
                    BC += settlementWeights[k]*B*C;
                    weight += settlementWeights[k];
                }
                var = B2*varianceY(0.0,expiryTime) + C2*varianceZ(0.0,expiryTime) + 2.0*rho_*BC*covarianceYZ(0.0,expiryTime);
                var /= weight; // not sure this makes sense
                return sqrt(var);
            }
            if (approxType==1) { // average log-volatility
                ActiveType varY  = varianceY(0.0,expiryTime);
                ActiveType varZ  = varianceZ(0.0,expiryTime);
                ActiveType covYZ = covarianceYZ(0.0,expiryTime);
                ActiveType vol   = 0.0;
                for (size_t k=0; k<N; ++k) {
                    B = exp(-a_*(settlementTimes[k]-expiryTime));
                    C = exp(-b_*(settlementTimes[k]-expiryTime));
                    vol += settlementWeights[k]*sqrt(B*B*varY + C*C*varZ + 2*rho_*B*C*covYZ);
                    weight += settlementWeights[k];
                }
                vol /= weight;
                return vol;
            }
            if (approxType==2) { // moment matching via with approximate variances
                // apply lognormal variance approximation exp{2mu + sigma^2}[exp{sigma^2}-1] \approx exp{2mu}sigma^2
                // assume perfect correlation of individual futures then vol of sum equals the sum of vols
                ActiveType varY  = varianceY(0.0,expiryTime);
                ActiveType varZ  = varianceZ(0.0,expiryTime);
                ActiveType covYZ = covarianceYZ(0.0,expiryTime);
                ActiveType vol   = 0.0;
                for (size_t k=0; k<N; ++k) {
                    // sum of futures
                    B = exp(-a_*(settlementTimes[k]-expiryTime));
                    C = exp(-b_*(settlementTimes[k]-expiryTime));
                    vol += settlementWeights[k]*exp(phi(settlementTimes[k]))*sqrt(B*B*varY + C*C*varZ + 2*rho_*B*C*covYZ);
                }
                // we have for the (average) future vol = exp{mu}sigma = E[F]exp{-sigma^2/2}sigma
                // setting vols equal to above and solving for sigma via Newton iteration
                ActiveType expectFuture = averageFuture(settlementTimes,settlementWeights);
                ActiveType volFuture0 = 0;
                ActiveType volFuture1 = vol/expectFuture;
                long cnt=0;
                while (fabs(volFuture1-volFuture0)>1.0e-8) {
                    volFuture0 = volFuture1;
                    ActiveType tmp = exp(-volFuture0*volFuture0/2.0);
                    ActiveType f   = tmp*volFuture0 - vol/expectFuture;
                    ActiveType fpr = tmp*(1.0-volFuture0*volFuture0);
                    volFuture1 = volFuture0 - f/fpr;
                    if ((++cnt)>10) break;  // we do not want to get stuck in the iteration if something goes wrong
                } // now we should have the volatility
                return volFuture1;
            }
            return 0.0; // default value
        }

        inline virtual ActiveType vanillaOption ( const DateType expiryTime, const VecD& settlementTimes, const VecP& settlementWeights, PassiveType strike, int callOrPut) {
            ActiveType fut = averageFuture(settlementTimes,settlementWeights);
            ActiveType vol = averageFutureStDev(expiryTime,settlementTimes,settlementWeights);
            ActiveType pv  = TemplateAuxilliaries::Black76(fut,strike,vol,1.0,callOrPut);
            return pv;
        }

    };

}

#undef _MIN_
#undef _MAX_

#endif  /* quantlib_template2fnormalmodel_hpp */
