/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2015 Sebastian Schlenkrich

*/

/*! \file ratespayoffT.hpp
    \brief generic payoff interface for MC simulation
    
*/


#ifndef quantlib_templatemcswap_hpp
#define quantlib_templatemcswap_hpp


#include <ql/experimental/templatemodels/montecarlo/mcpayoffT.hpp>
#include <ql/experimental/basismodels/swaptioncfs.hpp>

#include <ql/indexes/iborindex.hpp>
#include <ql/indexes/swapindex.hpp>
#include <ql/termstructures/yieldtermstructure.hpp>


namespace QuantLib {

    template <class DateType, class PassiveType, class ActiveType>
    class RatesPayoffT {
    protected:
        typedef MCPayoffT<DateType, PassiveType, ActiveType>  PayoffType;
        typedef MCSimulationT<DateType, PassiveType, ActiveType>  SimulationType;
        typedef typename MCSimulationT<DateType, PassiveType, ActiveType>::Path  PathType;

    public:

        // payoffs if we want to specify all the cash flows

        // generalised physically settled European swaption
        class GeneralSwaption : public MCPayoffT<DateType,PassiveType,ActiveType> {
        protected:
            std::vector<DateType>    floatTimes_;     // T_1, .., T_M
            std::vector<PassiveType> floatWeights_;   // u_1, .., u_M
            std::vector<DateType>    fixedTimes_;     // T_1, .., T_N
            std::vector<PassiveType> fixedWeights_;   // w_1, .., w_N
            PassiveType              strikeRate_;     // option strike
            PassiveType              payOrRec_;       // call (+1) or put (-1) option on swap rate
            inline void checkForConsistency() {
                // check consistency of swap
                // float leg
                QL_REQUIRE(floatWeights_.size()>0,"GeneralSwaption: empty float weights.");
                QL_REQUIRE(floatTimes_.size()==floatWeights_.size(),"GeneralSwaption: float sizes mismatch.");
                QL_REQUIRE(floatTimes_[0]>0,"GeneralSwaption: future float times required");
                for (size_t k=1; k<floatTimes_.size(); ++k) QL_REQUIRE(floatTimes_[k]>=floatTimes_[k-1],"GeneralSwaption: ascending float times required");
                // fixed leg
                QL_REQUIRE(fixedWeights_.size()>0,"GeneralSwaption: empty fixed weights.");
                QL_REQUIRE(fixedTimes_.size()==fixedWeights_.size(),"GeneralSwaption: fixed sizes mismatch.");
                QL_REQUIRE(fixedTimes_[0]>0,"GeneralSwaption: future fixed times required");
                for (size_t k=1; k<fixedTimes_.size(); ++k) QL_REQUIRE(fixedTimes_[k]>=fixedTimes_[k-1],"GeneralSwaption: ascending fixed times required");
                // finished
            }
        public:
            GeneralSwaption( DateType                        obsTime,    // observation equals fixing time
                             const std::vector<DateType>&    floatTimes,
                             const std::vector<PassiveType>& floatWeights,
                             const std::vector<DateType>&    fixedTimes,
                             const std::vector<PassiveType>& fixedWeights,
                             PassiveType                     strikeRate,
                             PassiveType                     payOrRec      )
                : PayoffType(obsTime),  floatTimes_(floatTimes), floatWeights_(floatWeights),
                  fixedTimes_(fixedTimes), fixedWeights_(fixedWeights), strikeRate_(strikeRate), payOrRec_(payOrRec) { 
                checkForConsistency();
            }
            GeneralSwaption( DateType                              obsTime,    // observation equals fixing time
                             const ext::shared_ptr<SwapIndex>&   swapIndex,
                             const Handle<YieldTermStructure>&     discYTS,
                             PassiveType                           strikeRate,
                             PassiveType                           payOrRec      )
                : PayoffType(obsTime), strikeRate_(strikeRate), payOrRec_(payOrRec) {
                Date today      = discYTS->referenceDate(); // check if this is the correct date...
                Date fixingDate = today + ((BigInteger)ClosestRounding(0)(obsTime*365.0)); // assuming act/365 day counting
                SwapCashFlows scf(swapIndex->underlyingSwap(fixingDate),discYTS,true);        // assume continuous tenor spreads
                // set attributes
                floatTimes_   = scf.floatTimes();
                floatWeights_ = scf.floatWeights();
                fixedTimes_   = scf.fixedTimes();
                fixedWeights_ = scf.annuityWeights();
                checkForConsistency();
            }

            virtual ~GeneralSwaption() = default;

            inline virtual ActiveType at(const ext::shared_ptr<PathType>& p) {
                ActiveType floatleg = 0.0;
                ActiveType annuity  = 0.0;
                // float leg
                for (size_t k=0; k<floatTimes_.size(); ++k) floatleg += floatWeights_[k] * p->zeroBond(PayoffType::observationTime(),floatTimes_[k]);
                // annuity
                for (size_t k=0; k<fixedTimes_.size(); ++k) annuity  += fixedWeights_[k] * p->zeroBond(PayoffType::observationTime(),fixedTimes_[k]);
                // floatleg - fixedleg...
                ActiveType res = floatleg - strikeRate_*annuity;
                // payer or receiver swap...
                res *= payOrRec_;
                // exercise option...
                res = (res>0) ? res : 0.0;
                return res;
            }
        };

        // future swap rate
        class SwapRate : public GeneralSwaption {
            // we save this to be able to clone the swap rate
            ext::shared_ptr<SwapIndex>     swapIndex_;
            Handle<YieldTermStructure>     discYTS_;
        public:
            SwapRate( DateType                        obsTime,    // observation equals fixing time
                      const std::vector<DateType>&    floatTimes,
                      const std::vector<PassiveType>& floatWeights,
                      const std::vector<DateType>&    fixedTimes,
                      const std::vector<PassiveType>& annuityWeights )
                      : GeneralSwaption(obsTime,floatTimes,floatWeights,fixedTimes,annuityWeights,0.0,0.0) {}

            SwapRate( const DateType                        fixingTime,
                      const ext::shared_ptr<SwapIndex>&     swapIndex,
                      const Handle<YieldTermStructure>&     discYTS    )
                      : GeneralSwaption(fixingTime,swapIndex,discYTS,0.0,0.0), swapIndex_(swapIndex), discYTS_(discYTS) { }

            virtual ~SwapRate() = default;

            inline virtual ActiveType at(const ext::shared_ptr<PathType>& p) {
                ActiveType floatleg = 0.0;
                ActiveType annuity  = 0.0;
                // float leg
                for (size_t k=0; k<GeneralSwaption::floatTimes_.size(); ++k) floatleg += GeneralSwaption::floatWeights_[k] * p->zeroBond(PayoffType::observationTime(),GeneralSwaption::floatTimes_[k]);
                // annuity
                for (size_t k=0; k<GeneralSwaption::fixedTimes_.size(); ++k) annuity  += GeneralSwaption::fixedWeights_[k] * p->zeroBond(PayoffType::observationTime(),GeneralSwaption::fixedTimes_[k]);
                return floatleg / annuity;
            }
            inline virtual ext::shared_ptr<PayoffType> at(const DateType t) {
                if ((swapIndex_!=0)&&(!discYTS_.empty())) return ext::shared_ptr<PayoffType>(new SwapRate(t, swapIndex_, discYTS_));
                QL_FAIL("Can not clone swap rate");
            }
            // payoff should NOT be discounted
            inline virtual ActiveType discountedAt(const ext::shared_ptr<PathType>& p) { return at(p); }
        };


        // Libor rate based on index
        class LiborRate : public MCPayoffT<DateType, PassiveType, ActiveType> {
        protected:
            // we save this to be able to clone the Libor rate
            ext::shared_ptr<IborIndex>     iborIndex_;
            Handle<YieldTermStructure>     discYTS_;
            DateType    fixingTime_, startTime_, endTime_;
            PassiveType oneOverDaycount_;
            PassiveType D_; // tenor basis
        public:
            LiborRate(const DateType                        fixingTime,
                      const ext::shared_ptr<IborIndex>&     iborIndex,
                      const Handle<YieldTermStructure>&     discYTS)
                      : PayoffType(fixingTime), fixingTime_(fixingTime), iborIndex_(iborIndex), discYTS_(discYTS) {
                Date today = discYTS->referenceDate(); // check if this is the correct date...
                Date fixingDate = today + ((BigInteger)ClosestRounding(0)(fixingTime_*365.0)); // assuming act/365 day counting
                Date startDate = iborIndex->valueDate(fixingDate);
                Date endDate = iborIndex->maturityDate(startDate);
                startTime_ = (startDate - today) / 365.0;
                endTime_ = (endDate - today) / 365.0;
                PassiveType liborForward = iborIndex->fixing(fixingDate, true);
                PassiveType daycount = iborIndex->dayCounter().yearFraction(startDate, endDate);
                oneOverDaycount_ = 1.0 / daycount;
                // tenor basis calculation
                D_ = (1.0 + daycount*liborForward) * discYTS->discount(endDate) / discYTS->discount(startDate);
            }
            LiborRate( const DateType                        fixingTime,
                       const DateType                        startTime,  // we don't get start and end date from the index
                       const DateType                        endTime,    // therefore we need to supply it explicitely
                       const ext::shared_ptr<IborIndex>&     iborIndex,
                       const Handle<YieldTermStructure>&     discYTS )
                       : PayoffType(fixingTime), fixingTime_(fixingTime), startTime_(startTime), endTime_(endTime) {
                Date today      = discYTS->referenceDate(); // check if this is the correct date...
                Date fixingDate = today + ((BigInteger)ClosestRounding(0)(fixingTime_*365.0)); // assuming act/365 day counting
                Date startDate  = today + ((BigInteger)ClosestRounding(0)(startTime_ *365.0)); // assuming act/365 day counting
                Date endDate    = today + ((BigInteger)ClosestRounding(0)(endTime_   *365.0)); // assuming act/365 day counting
                PassiveType liborForward = iborIndex->fixing(fixingDate,true);
                PassiveType daycount     = iborIndex->dayCounter().yearFraction(startDate,endDate);
                oneOverDaycount_ = 1.0/daycount;
                // tenor basis calculation
                D_ = (1.0 + daycount*liborForward) * discYTS->discount(endDate) / discYTS->discount(startDate);
            }
            virtual ~LiborRate() = default;
            inline virtual ActiveType at(const ext::shared_ptr<PathType>& p) {
                return ( p->zeroBond(fixingTime_,startTime_) / p->zeroBond(fixingTime_,endTime_) * D_ - 1.0 ) * oneOverDaycount_;
            }
            inline virtual ext::shared_ptr<PayoffType> at(const DateType t) {
                if ((iborIndex_!=0)&&(!discYTS_.empty())) return ext::shared_ptr<PayoffType>(new LiborRate(t, iborIndex_, discYTS_));
                QL_FAIL("Can not clone Libor rate");
            }
        };

        // Libor rate based for hybrid models
        class LiborRateCcy : public LiborRate {
        private:
            std::string alias_;
        public:
            LiborRateCcy(const DateType                        fixingTime,
                         const ext::shared_ptr<IborIndex>&     iborIndex,
                         const Handle<YieldTermStructure>&     discYTS,
                         const std::string                     alias)
            : LiborRate(fixingTime, iborIndex, discYTS), alias_(alias) {}
            virtual ~LiborRateCcy() = default;
            inline virtual ActiveType at(const ext::shared_ptr<PathType>& p) {
                return (p->zeroBond(LiborRate::fixingTime_, LiborRate::startTime_, alias_) / 
                        p->zeroBond(LiborRate::fixingTime_, LiborRate::endTime_, alias_) * LiborRate::D_ - 1.0) * 
                       LiborRate::oneOverDaycount_;
            }
            inline virtual ext::shared_ptr<PayoffType> at(const DateType t) {
                if ((LiborRate::iborIndex_ != 0) && (!LiborRate::discYTS_.empty()))
                    return ext::shared_ptr<PayoffType>(new LiborRateCcy(t, LiborRate::iborIndex_, LiborRate::discYTS_, alias_));
                QL_FAIL("Can not clone Libor rate");
            }

        };


        //class SwapRate : public MCPayoffT<DateType, PassiveType, ActiveType> {
        //protected:
        //	ext::shared_ptr<MCPayoffT<DateType, PassiveType, ActiveType>::SwapRate > swaprate_;
        //public:
        //	SwapRate( const DateType                        fixingTime,
        //			  const ext::shared_ptr<SwapIndex>&   swapIndex,
        //			  const Handle<YieldTermStructure>&     discYTS    ) : TemplateMCPayoff(fixingTime) {
        //		Date today      = discYTS->referenceDate(); // check if this is the correct date...
        //		Date fixingDate = today + ((BigInteger)ClosestRounding(0)(fixingTime*365.0)); // assuming act/365 day counting
        //		SwapCashFlows scf(swapIndex->underlyingSwap(fixingDate),discYTS,true);  // assume continuous tenor spreads
        //		swaprate_ = ext::shared_ptr<TemplateMCPayoff::SwapRate>(new TemplateMCPayoff::SwapRate(
        //			fixingTime, scf.floatTimes(), scf.floatWeights(), scf.fixedTimes(), scf.annuityWeights() ) );
        //	}
        //	inline virtual ActiveType at(const ext::shared_ptr<PathType>& p) { return swaprate_->at(p); }
        //};


        // CashFlow decorating a payoff with start and pay date for organisation in legs
        class CashFlow : public MCPayoffT<DateType, PassiveType, ActiveType> {
        protected:
            ext::shared_ptr<PayoffType> x_;
            DateType                     startTime_;  // on exercise only cash flows with startTime >= exerciseTime will be considered
            DateType                     payTime_;
            bool                         applyZCBAdjuster_;
        public:
            CashFlow ( const ext::shared_ptr<PayoffType>&    x,
                       const DateType                        startTime,
                       const DateType                        payTime,
                       const bool                            applyZCBAdjuster = false)
                : PayoffType(payTime), x_(x), startTime_(startTime), payTime_(payTime), applyZCBAdjuster_(applyZCBAdjuster) {}
            CashFlow ( const ext::shared_ptr<PayoffType>&   x,
                       const bool                            applyZCBAdjuster = false)
                : PayoffType(x->observationTime()), x_(x),
                  startTime_(x->observationTime()), payTime_(x->observationTime()), applyZCBAdjuster_(applyZCBAdjuster) {}
            virtual ~CashFlow() = default;
            inline virtual ActiveType at(const ext::shared_ptr<PathType>& p) { 
                return (applyZCBAdjuster_) ? (p->zeroBond(payTime_,payTime_)*x_->at(p)) : (x_->at(p))  ;
            }
            // inspectors
            inline const DateType startTime() const { return startTime_; }
            inline const DateType payTime()   const { return payTime_;   }
            inline virtual std::set<DateType> observationTimes() { return unionTimes(PayoffType::observationTimes(), x_->observationTimes()); }
        };

        // a CashFlow leg as a ordered list of CashFlows
        class Leg : public std::vector< ext::shared_ptr<CashFlow> > {
        public:
            static bool firstIsLess( const ext::shared_ptr<CashFlow>& first, const ext::shared_ptr<CashFlow>& second ) { return first->startTime() < second->startTime(); } 
            Leg ( const std::vector< ext::shared_ptr<CashFlow> >&   cashflows ) : std::vector< ext::shared_ptr<CashFlow> >(cashflows.begin(),cashflows.end())  {
                // check that all CashFlows are consistent
                // sort by start time
                std::sort(this->begin(),this->end(),firstIsLess );
            }
        };

        // a swap as a set of CashFlow legs (e.g. structured, funding, notional exchanges)
        class Swap : public std::vector< ext::shared_ptr<Leg> > {
        public:
            Swap ( const std::vector< ext::shared_ptr<Leg> >&   legs ) : std::vector< ext::shared_ptr<Leg> >(legs.begin(),legs.end())  { }
        };


        // this is the key structure for AMC valuation
        class CancellableNote {
        private:
            // underlying
            std::vector< ext::shared_ptr<Leg> >  underlyings_;       // the underlying CashFlow legs
            // call features
            std::vector< DateType >                callTimes_;            // exercise times
            std::vector< ext::shared_ptr<Leg> >  earlyRedemptions_;     // strikes payed at exercise
            std::vector< ext::shared_ptr<Leg> >  regressionVariables_;  // regression variables at exercise
        public:
            CancellableNote ( const std::vector< ext::shared_ptr<Leg> >&    underlyings,
                              const std::vector< DateType >&                callTimes,
                              const std::vector< ext::shared_ptr<Leg> >&    earlyRedemptions,
                              const std::vector< ext::shared_ptr<Leg> >&    regressionVariables )
                              : underlyings_(underlyings), callTimes_(callTimes), earlyRedemptions_(earlyRedemptions), regressionVariables_(regressionVariables) {
                // sanity checks
            }
            // inspectors
            inline const std::vector< ext::shared_ptr<Leg> >& underlyings()           const { return underlyings_;         }
            inline const std::vector< DateType >&               callTimes()           const { return callTimes_;           }
            inline const std::vector< ext::shared_ptr<Leg> >& earlyRedemptions()      const { return earlyRedemptions_;    }
            inline const std::vector< ext::shared_ptr<Leg> >& regressionVariables()   const { return regressionVariables_; }
        };


        // further payoffs

        // annuity
        class Annuity : public MCPayoffT<DateType, PassiveType, ActiveType> {
        protected:
            std::vector<DateType>    payTimes_;
            std::vector<PassiveType> payWeights_;  // these are typically year fractions
        public:
            Annuity( DateType                        obsTime, 
                     const std::vector<DateType>&    payTimes,
                     const std::vector<PassiveType>& payWeights)
                : PayoffType(obsTime), payTimes_(payTimes), payWeights_(payWeights) { }
            virtual ~Annuity() = default;
            inline virtual ActiveType at(const ext::shared_ptr<PathType>& p) {
                ActiveType res = 0.0;
                size_t N = (payTimes_.size()<payWeights_.size()) ? (payTimes_.size()) : (payWeights_.size());
                for (size_t k=0; k<N; ++k) {
                    if (payTimes_[k]>PayoffType::observationTime()) {
                        res += payWeights_[k] * p->zeroBond(PayoffType::observationTime(),payTimes_[k]);
                    }
                }
                return res;
            }
        };

        // prototypical physically settled European swaption
        class ModelSwaption : public MCPayoffT<DateType, PassiveType, ActiveType> {
        protected:
            std::vector<DateType>    times_;        // T_0, .., T_N
            std::vector<PassiveType> payWeights_;   // tau_0, .., tau_N-1
            PassiveType              strikeRate_;   // option strike
            PassiveType              payOrRec_;     // call (+1) or put (-1) option on swap rate
            bool                     isConsistent_; // check consistency of input
        public:
            ModelSwaption( DateType                        obsTime, 
                           const std::vector<DateType>&    times,
                           const std::vector<PassiveType>& payWeights,
                           PassiveType                     strikeRate,
                           PassiveType                     payOrRec      )
                : PayoffType(obsTime), times_(times), payWeights_(payWeights), strikeRate_(strikeRate), payOrRec_(payOrRec) { 
                isConsistent_ = true;
                if (times_.size()<2) isConsistent_ = false;
                for (size_t k=0; k<times_.size(); ++k) if (times_[k]<PayoffType::observationTime()) isConsistent_ = false;
                // default weights
                if (payWeights_.size()!=times_.size()-1) {
                    payWeights_.resize(times_.size()-1);
                    for (size_t k=0; k<payWeights_.size(); ++k) payWeights_[k] = times_[k+1] - times_[k];
                }
                // finished
            }
            virtual ~ModelSwaption() = default;
            inline virtual ActiveType at(const ext::shared_ptr<PathType>& p) {
                ActiveType res = 0.0;
                if (!isConsistent_) return res;
                // annuity...
                for (size_t k=0; k<payWeights_.size(); ++k) res += payWeights_[k] * p->zeroBond(PayoffType::observationTime(),times_[k+1]);
                // floatleg - fixedleg...
                res = p->zeroBond(PayoffType::observationTime(),times_[0]) - p->zeroBond(PayoffType::observationTime(),times_[times_.size()-1]) - strikeRate_*res;
                // payer or receiver swap...
                res *= payOrRec_;
                // exercise option...
                res = (res>0) ? res : 0.0;
                return res;
            }
        };



        // undiscounted correlation between prototypical physically settled European swaption
        class ModelCorrelation : public MCPayoffT<DateType, PassiveType, ActiveType> {
        protected:
            std::vector<DateType> times_;
            DateType T1_, T2_;
            ActiveType swapRate(const ext::shared_ptr<PathType>& p, const DateType t, const DateType TN ) {
                ActiveType num = p->zeroBond(t,t) - p->zeroBond(t,TN);
                ActiveType den = 0.0;
                for (ActiveType Ti = t; Ti<TN; Ti+=1.0) {
                    ActiveType T = (Ti+1.0>TN) ? (TN) : (Ti+1.0);
                    den += (T-Ti) * p->zeroBond(t,T);
                }
                return num / den;
            }
        public:
            ModelCorrelation( const std::vector<DateType>&    times,   // observation times
                              const DateType                  T1,      // swap term one
                              const DateType                  T2 )     // swap term two
                              : PayoffType(0.0), times_(times), T1_(T1), T2_(T2) {
                QL_REQUIRE(times_.size()>1,"ModelCorrelation: At least two observation times required.");
            }
            virtual ~ModelCorrelation() = default;
            // payoff should NOT be discounted
            inline virtual ActiveType discountedAt(const ext::shared_ptr<PathType>& p) { return at(p); }
            inline virtual ActiveType at(const ext::shared_ptr<PathType>& p) {
                std::vector<ActiveType> dS1(times_.size()-1), dS2(times_.size()-1);
                ActiveType EdS1 = 0.0, EdS2 = 0.0; 
                for (size_t i=1; i<times_.size(); ++i) {
                    dS1[i-1] =  swapRate(p,times_[i],times_[i]+T1_) - swapRate(p,times_[i-1],times_[i-1]+T1_);
                    dS2[i-1] =  swapRate(p,times_[i],times_[i]+T2_) - swapRate(p,times_[i-1],times_[i-1]+T2_);
                    EdS1 += dS1[i-1];
                    EdS2 += dS2[i-1];
                }
                EdS1 /= dS1.size();
                EdS2 /= dS2.size();
                ActiveType Var1=0.0, Var2=0.0, Cov=0.0;
                for (size_t i=0; i<times_.size()-1; ++i) {
                    Var1 += (dS1[i] - EdS1)*(dS1[i] - EdS1);
                    Var2 += (dS2[i] - EdS2)*(dS2[i] - EdS2);
                    Cov  += (dS1[i] - EdS1)*(dS2[i] - EdS2);
                }
                return Cov / sqrt(Var1*Var2);
            }
            inline virtual std::set<DateType> observationTimes() {
                std::set<DateType> s;
                for (size_t k = 0; k < times_.size(); ++k) s.insert(times_[k]);
                return s;
            }
        };

        // undiscounted correlation between forward rates
        class ForwardRateCorrelation : public  MCPayoffT<DateType, PassiveType, ActiveType> {
        protected:
            std::vector<DateType> times_;
            DateType T1_, Term1_; 
            DateType T2_, Term2_;
            ActiveType fwSwapRate(const ext::shared_ptr<PathType>& p, const DateType t, const DateType TSettle, const DateType Term ) {
                ActiveType num = p->zeroBond(t,TSettle) - p->zeroBond(t,TSettle+Term);
                ActiveType den = 0.0;
                for (ActiveType Ti = TSettle; Ti<TSettle+Term; Ti+=1.0) {
                    ActiveType T = (Ti+1.0>TSettle+Term) ? (TSettle+Term) : (Ti+1.0);
                    den += (T-Ti) * p->zeroBond(t,T);
                }
                return num / den;
            }
            ActiveType fraRate(const ext::shared_ptr<PathType>& p, const DateType t, const DateType TSettle, const DateType Term ) {
                ActiveType rate = (p->zeroBond(t,TSettle) / p->zeroBond(t,TSettle+Term) - 1.0)/Term;
                return rate;
            }
        public:
            ForwardRateCorrelation( const std::vector<DateType>&    times,   // observation times
                                    const DateType                  T1,      // fixing date one
                                    const DateType                  Term1,   // tenor one
                                    const DateType                  T2,      // fixing date two
                                    const DateType                  Term2)   // tenor two
                                    : PayoffType(0.0), times_(times), T1_(T1), Term1_(Term1), T2_(T2), Term2_(Term2) {
                QL_REQUIRE(times_.size()>1,"ModelCorrelation: At least two observation times required.");
            }
            virtual ~ForwardRateCorrelation() = default;
            // payoff should NOT be discounted
            inline virtual ActiveType discountedAt(const ext::shared_ptr<PathType>& p) { return at(p); }
            inline virtual ActiveType at(const ext::shared_ptr<PathType>& p) {
                std::vector<ActiveType> dS1(times_.size()-1), dS2(times_.size()-1);
                ActiveType EdS1 = 0.0, EdS2 = 0.0; 
                for (size_t i=1; i<times_.size(); ++i) {
                    dS1[i-1] =  fraRate(p,times_[i],T1_,Term1_) - fraRate(p,times_[i-1],T1_,Term1_);
                    dS2[i-1] =  fraRate(p,times_[i],T2_,Term2_) - fraRate(p,times_[i-1],T2_,Term2_);
                    EdS1 += dS1[i-1];
                    EdS2 += dS2[i-1];
                }
                EdS1 /= dS1.size();
                EdS2 /= dS2.size();
                ActiveType Var1=0.0, Var2=0.0, Cov=0.0;
                for (size_t i=0; i<times_.size()-1; ++i) {
                    Var1 += (dS1[i] - EdS1)*(dS1[i] - EdS1);
                    Var2 += (dS2[i] - EdS2)*(dS2[i] - EdS2);
                    Cov  += (dS1[i] - EdS1)*(dS2[i] - EdS2);
                }
                return Cov / sqrt(Var1*Var2);
            }
            inline virtual std::set<DateType> observationTimes() {
                std::set<DateType> s;
                for (size_t k = 0; k < times_.size(); ++k) s.insert(times_[k]);
                return s;
            }
        };

    }; // class RatesPayoffT

}

#endif  /* ifndef quantlib_templatemcswap_hpp */ 
