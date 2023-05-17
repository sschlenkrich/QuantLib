/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2016 Sebastian Schlenkrich

*/

/*! \file commoditypayoffT.hpp
    \brief MC payoffs 
    
*/


#ifndef quantlib_templatecommoditypayoffs_hpp
#define quantlib_templatecommoditypayoffs_hpp


#include <ql/experimental/templatemodels/montecarlo/mcpayoffT.hpp>


namespace QuantLib {

    template <class DateType, class PassiveType, class ActiveType>
    class CommodityPayoffT {
    protected:
        typedef MCPayoffT<DateType, PassiveType, ActiveType>  PayoffType;
        typedef MCSimulationT<DateType, PassiveType, ActiveType>                 SimulationType;
        typedef typename MCSimulationT<DateType, PassiveType, ActiveType>::Path  PathType;

    public:

        // call or put on a strip of futures
        class AverageFutureOption : public MCPayoffT<DateType,PassiveType,ActiveType> {
            std::vector<DateType>    settlementTimes_;
            std::vector<PassiveType> settlementWeights_;
            PassiveType              strike_;        // option strike
            PassiveType              callOrPut_;     // call (+1) or put (-1) option on swap rate
        public:
            AverageFutureOption( DateType                        obsTime, 
                                 const std::vector<DateType>&    settlementTimes, 
                                 const std::vector<PassiveType>& settlementWeights,
                                 const PassiveType               strike,
                                 const PassiveType               callOrPut )
                : PayoffType(obsTime), settlementTimes_(settlementTimes), settlementWeights_(settlementWeights), strike_(strike), callOrPut_(callOrPut) { }
            inline virtual ActiveType at(const ext::shared_ptr<PathType>& p) {
                ActiveType fut=0.0;
                for (size_t k=0; k<settlementTimes_.size(); ++k)
                    if (settlementTimes_[k]>=PayoffType::observationTime()) fut += settlementWeights_[k] * p->futureAsset(PayoffType::observationTime(),settlementTimes_[k],"");
                ActiveType V  = callOrPut_ * (fut - strike_);
                return (V>0.0) ? (V) : ((ActiveType)0.0);
            }
        };

        // covariance between quoted futures with roll-over
        class AverageFutureCovariance : public MCPayoffT<DateType,PassiveType,ActiveType> {
        protected:
            std::vector<DateType> obsTimes_;
            // future A definition
            std::vector<DateType>    settlementTimesA_;     // delivery period
            std::vector<PassiveType> settlementWeightsA_;   // usually equal weights for averaging
            PassiveType              obsLagA_;              // roll over delivery period if first date is less than observation lag
            // future B definition
            std::vector<DateType>    settlementTimesB_;     // delivery period
            std::vector<PassiveType> settlementWeightsB_;   // usually equal weights for averaging
            PassiveType              obsLagB_;              // roll over delivery period if first date is less than observation lag
            // further flags
            bool                     useLogReturns_;        // specify (normal vs lognormal) type of returns considered
            long                     calcType_;             // flag to distinguish what to do (quick and dirty)
                                                            // 0 - calculate covariance
                                                            // 1 - calculate correlation instead of covariance
                                                            // 2 - calculate spread variance 
            // helper methods
            void rollOver( std::vector<DateType>& settlementTimes ) {
                PassiveType rollPeriod = settlementTimes[settlementTimes.size()-1] - settlementTimes[0] + 1.0/365.0;
                for (size_t k=0; k<settlementTimes.size(); ++k) settlementTimes[k] += rollPeriod;
            }
            ActiveType averageFuture(const ext::shared_ptr<PathType>& p, const DateType t, const std::vector<DateType>& settlementTimes, const std::vector<PassiveType>& settlementWeights ) {
                ActiveType fut=0.0;
                for (size_t k=0; k<settlementTimes.size(); ++k) if (settlementTimes[k]>=t) fut += settlementWeights[k] * p->futureAsset(t,settlementTimes[k],"");
                return fut;
            }
        public:
            AverageFutureCovariance( const std::vector<DateType>&    obsTimes,
                                     const std::vector<DateType>&    settlementTimesA,     
                                     const std::vector<PassiveType>& settlementWeightsA,
                                     const PassiveType               obsLagA,
                                     const std::vector<DateType>&    settlementTimesB,
                                     const std::vector<PassiveType>& settlementWeightsB,
                                     const PassiveType               obsLagB,
                                     const bool                      useLogReturns,       
                                     const long                      calcType )
            : PayoffType(0.0), obsTimes_(obsTimes), settlementTimesA_(settlementTimesA), settlementWeightsA_(settlementWeightsA), obsLagA_(obsLagA),
            settlementTimesB_(settlementTimesB), settlementWeightsB_(settlementWeightsB), obsLagB_(obsLagB), useLogReturns_(useLogReturns), calcType_(calcType) {}
            // payoff should NOT be discounted
            inline virtual ActiveType discountedAt(const ext::shared_ptr<PathType>& p) { return at(p); }
            inline virtual ActiveType at(const ext::shared_ptr<PathType>& p) {
                // if only one observation time then return [FutA * FutB] for autocorrelation estimation
                if (obsTimes_.size()==1) return averageFuture(p,settlementTimesA_[0],settlementTimesA_,settlementWeightsA_) * averageFuture(p,settlementTimesB_[0],settlementTimesB_,settlementWeightsB_);
                std::vector<ActiveType> dFutA(obsTimes_.size()-1), dFutB(obsTimes_.size()-1), dSprd(obsTimes_.size()-1);
                ActiveType EdFutA = 0.0, EdFutB = 0.0, EdSprd = 0.0; 
                for (size_t i=1; i<obsTimes_.size(); ++i) {
                    // check rollover A and B
                    if (settlementTimesA_[0]-obsTimes_[i]<obsLagA_) rollOver(settlementTimesA_);
                    if (settlementTimesB_[0]-obsTimes_[i]<obsLagB_) rollOver(settlementTimesB_);
                    // calculate returns
                    if (useLogReturns_) {
                        dFutA[i-1] = log(averageFuture(p,obsTimes_[i],settlementTimesA_,settlementWeightsA_)) - log(averageFuture(p,obsTimes_[i-1],settlementTimesA_,settlementWeightsA_));
                        dFutB[i-1] = log(averageFuture(p,obsTimes_[i],settlementTimesB_,settlementWeightsB_)) - log(averageFuture(p,obsTimes_[i-1],settlementTimesB_,settlementWeightsB_));
                    } else {
                        dFutA[i-1] = averageFuture(p,obsTimes_[i],settlementTimesA_,settlementWeightsA_) - averageFuture(p,obsTimes_[i-1],settlementTimesA_,settlementWeightsA_);
                        dFutB[i-1] = averageFuture(p,obsTimes_[i],settlementTimesB_,settlementWeightsB_) - averageFuture(p,obsTimes_[i-1],settlementTimesB_,settlementWeightsB_);
                        if (calcType_==2) dSprd[i-1] = dFutB[i-1] - dFutA[i-1];
                    }
                    EdFutA += dFutA[i-1];
                    EdFutB += dFutB[i-1];
                }
                EdFutA /= dFutA.size();
                EdFutB /= dFutB.size();
                if (calcType_==2) EdSprd = EdFutB - EdFutA;
                ActiveType VarA=0.0, VarB=0.0, Cov=0.0, SprdVar=0.0;
                for (size_t i=0; i<obsTimes_.size()-1; ++i) {
                    if (calcType_<2) {
                        VarA += (dFutA[i] - EdFutA)*(dFutA[i] - EdFutA);
                        VarB += (dFutB[i] - EdFutB)*(dFutB[i] - EdFutB);
                        Cov  += (dFutA[i] - EdFutA)*(dFutB[i] - EdFutB);
                    } else {
                        SprdVar += (dSprd[i]-EdSprd)*(dSprd[i]-EdSprd);
                    }
                }
                if (calcType_==0) return Cov     / (dFutA.size()-1);    // covariance
                if (calcType_==1) return Cov     / sqrt(VarA*VarB);     // correlation
                if (calcType_==2) return SprdVar / (dSprd.size()-1);  // correlation
                return 0.0; // fall back
            }

        };

    }; // class CommodityPayoffT

}

#endif  /* ifndef quantlib_templatecommoditypayoffs_hpp */ 
