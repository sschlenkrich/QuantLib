/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2014 Sebastian Schlenkrich

 */

/*! \file basisswap.hpp
    \brief (Cross currency) Interest rate swap consisting of several legs
*/


#include <ql/experimental/basismodels/liborbasisswapratehelpers.hpp>

namespace QuantLib {

    // dummy
    namespace TenorSwapRateHelperDummyNamespace {
        void no_deletion(YieldTermStructure*) {}
    }


    TenorSwapRateHelper::TenorSwapRateHelper( 
                             Real                          rate,
                             const Period&                 tenor,
                             // swap conventions
                             const Period&                 fwdStart,
                             const Calendar&               paymentCalendar,
                             BusinessDayConvention         paymentBDC,
                             bool                          spreadOnRecLeg,
                             // pay leg details are taken from IborIndex
                             const ext::shared_ptr<IborIndex>&  payIndex,
                             // rec leg details are taken from IborIndex
                             const ext::shared_ptr<IborIndex>&  recIndex,
                             // discount curve, forward curves are taken from Ibor indices
                             const Handle<YieldTermStructure>& discountCurve )
            : BasisSwapRateHelper(rate),
              tenor_(tenor), fwdStart_(fwdStart),
              paymentCalendar_(paymentCalendar), paymentBDC_(paymentBDC), spreadOnRecLeg_(spreadOnRecLeg) {

            // define bootstrapping curves
            useForDiscount_   = discountCurve.empty();
            useForPayForward_ = payIndex->forwardingTermStructure().empty();
            useForRecForward_ = recIndex->forwardingTermStructure().empty();
            // we need to check that we set at least one curve for bootstrapping
            QL_REQUIRE(useForDiscount_ || useForPayForward_ || useForRecForward_,
                "No yield curve to be set for bootstrapping.");
            // we need to set up indices and discount curves with relinkable handles to set them appropriately later
            if (!useForDiscount_)   discountRelinkableHandle_ = RelinkableHandle<YieldTermStructure>(discountCurve.currentLink());
            if (!useForPayForward_) payRelinkableHandle_      = RelinkableHandle<YieldTermStructure>(payIndex->forwardingTermStructure().currentLink());
            if (!useForRecForward_) recRelinkableHandle_      = RelinkableHandle<YieldTermStructure>(recIndex->forwardingTermStructure().currentLink());
                
            payIndex_                 = payIndex->clone(payRelinkableHandle_);
            recIndex_                 = recIndex->clone(recRelinkableHandle_);

            // make sure ratehelper reacts on changes
            //registerWith(payIndex_);
            //registerWith(recIndex_);
            //registerWith(discountRelinkableHandle_);
            initializeDates();
        }

    void TenorSwapRateHelper::initializeDates(){
        // set up reference swap
        // first we need a fixed and floating schedude
        Date referenceDate = Settings::instance().evaluationDate();
        Natural fixingDays = std::max(payIndex_->fixingDays(),recIndex_->fixingDays());
        Date spotDate      = paymentCalendar_.advance(referenceDate, fixingDays*Days);
        Date startDate     = spotDate+fwdStart_;
        Date endDate       = startDate+tenor_;
        Schedule   paySchedule = MakeSchedule()
                                 .from(startDate)
                                 .to(endDate)
                                 .withTenor(      payIndex_->tenor()                 )
                                 .withCalendar(   payIndex_->fixingCalendar()        )
                                 .withConvention( payIndex_->businessDayConvention() );
        Schedule   recSchedule = MakeSchedule()
                                 .from(startDate)
                                 .to(endDate)
                                 .withTenor(      recIndex_->tenor()                 )
                                 .withCalendar(   recIndex_->fixingCalendar()        )
                                 .withConvention( recIndex_->businessDayConvention() );
       Leg payLeg = IborLeg(paySchedule,payIndex_)
                    .withNotionals(1.0)
                    .withPaymentDayCounter(payIndex_->dayCounter())
                    .withPaymentAdjustment(paymentBDC_)
                    .withFixingDays(payIndex_->fixingDays());
       Leg recLeg = IborLeg(recSchedule,recIndex_)
                    .withNotionals(1.0)
                    .withPaymentDayCounter(recIndex_->dayCounter())
                    .withPaymentAdjustment(paymentBDC_)
                    .withFixingDays(recIndex_->fixingDays());
       std::vector<Leg> legs(2);
       legs[0]  = payLeg;
       legs[1]  = recLeg;
       std::vector<bool> payer(2);
       payer[0] = true;
       payer[1] = false;
       basisSwap_ = ext::shared_ptr<BasisSwap>(new BasisSwap(legs, payer, (spreadOnRecLeg_ ? 1 : 0), true));
       earliestDate_ = basisSwap_->startDate();
       latestDate_   = basisSwap_->maturityDate();
       // we set up pricing engine here coz it's associated to instrument
       basisSwap_->setPricingEngine(ext::shared_ptr<PricingEngine>(new DiscountingSwapEngine(discountRelinkableHandle_, false)));
    }

    void TenorSwapRateHelper::setTermStructure(YieldTermStructure* ts) {
        // this initialisation is intended to be done only once
        ext::shared_ptr<YieldTermStructure> temp(ts, TenorSwapRateHelperDummyNamespace::no_deletion);
        if (useForDiscount_)    discountRelinkableHandle_.linkTo(temp,false);
        if (useForPayForward_)  payRelinkableHandle_.linkTo(temp,false);
        if (useForRecForward_)  recRelinkableHandle_.linkTo(temp,false);
        // parent class method...
        RelativeDateRateHelper::setTermStructure(ts);
    }

    Real TenorSwapRateHelper::impliedQuote() const {
        QL_REQUIRE(termStructure_ != 0, "term structure not set");
        // we didn't register as observers - force calculation
        basisSwap_->recalculate();
        return basisSwap_->fairRate();
    }


    XCCYSwapRateHelper::XCCYSwapRateHelper(
                             Real                          rate,
                             const Period&                 tenor,
                             // swap conventions
                             const Period&                 fwdStart,
                             const Calendar&               spotStartCalendar,
                             BusinessDayConvention         payBDC,
                             BusinessDayConvention         recBDC,
                             bool                          spreadOnRecLeg,
                             // pay leg details are taken from IborIndex
                             const ext::shared_ptr<IborIndex>&  payIndex,
                             // rec leg details are taken from IborIndex
                             const ext::shared_ptr<IborIndex>&  recIndex,
                             // discount curve, forward curves are taken from Ibor indices
                             const Handle<YieldTermStructure>&    payDiscCurve,
                             const Handle<YieldTermStructure>&    recDiscCurve,
                             // today's fx rates FOR/DOM
                             Real                          payFxForDom,
                             Real                          recFxForDom,
                             // should the swap be fx resetable on non-spreaded leg
                             bool                          fxResetable
                             )
        : BasisSwapRateHelper(rate), tenor_(tenor), fwdStart_(fwdStart),
          spotStartCalendar_(spotStartCalendar), payBDC_(payBDC), recBDC_(recBDC), spreadOnRecLeg_(spreadOnRecLeg),
          payFxForDom_(payFxForDom), recFxForDom_(recFxForDom), fxResetable_(fxResetable) {
        // define bootstrapping curves
        useForPayDiscount_   = payDiscCurve.empty();
        useForRecDiscount_   = recDiscCurve.empty();
        useForPayForward_    = payIndex->forwardingTermStructure().empty();
        useForRecForward_    = recIndex->forwardingTermStructure().empty();
        // we need to check that we set at least one curve for bootstrapping
        QL_REQUIRE(useForPayDiscount_ || useForRecDiscount_ || useForPayForward_ || useForRecForward_,
            "No yield curve to be set for bootstrapping.");
        // we need to set up indices and discount curves with relinkable handles to set them appropriately later
        if (!useForPayDiscount_) payDisRelinkableHandle_ = RelinkableHandle<YieldTermStructure>(payDiscCurve.currentLink());
        if (!useForRecDiscount_) recDisRelinkableHandle_ = RelinkableHandle<YieldTermStructure>(recDiscCurve.currentLink());
        if (!useForPayForward_)  payForRelinkableHandle_ = RelinkableHandle<YieldTermStructure>(payIndex->forwardingTermStructure().currentLink());
        if (!useForRecForward_)  recForRelinkableHandle_ = RelinkableHandle<YieldTermStructure>(recIndex->forwardingTermStructure().currentLink());		
        payIndex_ = payIndex->clone(payForRelinkableHandle_);
        recIndex_ = recIndex->clone(recForRelinkableHandle_);

        // make sure ratehelper reacts on changes
        //registerWith(payIndex_);
        //registerWith(recIndex_);
        //registerWith(payDisRelinkableHandle_);
        //registerWith(recDisRelinkableHandle_);
        initializeDates();
    }

    void XCCYSwapRateHelper::initializeDates() {
        // we may set up a non-resetable swap to get start/end date
        // for resetable swaps we need the disc curve for non-spread leg available
        // hence the swap may need to be updated within each bootstrapping iteration
        setUpBasisSwap(false);
        earliestDate_ = basisSwap_->startDate();
        latestDate_   = basisSwap_->maturityDate();
    }

    void XCCYSwapRateHelper::setTermStructure(YieldTermStructure* ts) {
        // this initialisation is intended to be done only once
        ext::shared_ptr<YieldTermStructure> temp(ts, TenorSwapRateHelperDummyNamespace::no_deletion);
        if (useForPayDiscount_) payDisRelinkableHandle_.linkTo(temp,false);
        if (useForRecDiscount_) recDisRelinkableHandle_.linkTo(temp,false);
        if (useForPayForward_)  payForRelinkableHandle_.linkTo(temp,false);
        if (useForRecForward_)  recForRelinkableHandle_.linkTo(temp,false);
        // parent class method...
        RelativeDateRateHelper::setTermStructure(ts);
    }

    Real XCCYSwapRateHelper::impliedQuote() const {
        QL_REQUIRE(termStructure_ != 0, "term structure not set");
        // in case of resetables we need to recalculate notionals after curve change
        if (fxResetable_) setUpBasisSwap(true);
        // we didn't register as observers - force calculation
        basisSwap_->recalculate();
        return basisSwap_->fairRate();
    }


    void XCCYSwapRateHelper::setUpBasisSwap(bool fxResetable) const {
        // set up reference swap
        // first we need a fixed and floating schedude
        Date referenceDate = Settings::instance().evaluationDate();
        Natural fixingDays = std::max(payIndex_->fixingDays(),recIndex_->fixingDays());
        Date spotDate      = spotStartCalendar_.advance(referenceDate, fixingDays*Days);
        Date startDate     = spotDate+fwdStart_;
        Date endDate       = startDate+tenor_;
        Schedule   paySchedule = MakeSchedule()
                                 .from(startDate)
                                 .to(endDate)
                                 .withTenor(      payIndex_->tenor()                 )
                                 .withCalendar(   payIndex_->fixingCalendar()        )
                                 .withConvention( payIndex_->businessDayConvention() );
        Schedule   recSchedule = MakeSchedule()
                                 .from(startDate)
                                 .to(endDate)
                                 .withTenor(      recIndex_->tenor()                 )
                                 .withCalendar(   recIndex_->fixingCalendar()        )
                                 .withConvention( recIndex_->businessDayConvention() );
        // we consider 1 unit domestic currency as notional and exchange it to rec/pay foreign notional
        // this aims at keeping potential future flexibility for further developments
        // note, (current) par spread valuation is independent of fx
        // initially we assume non-resetable swap with fx(T_0)=fx(t)!!
        std::vector<Real> payNotionals(paySchedule.size()-1,1.0/payFxForDom_);
        std::vector<Real> recNotionals(recSchedule.size()-1,1.0/recFxForDom_);
        Leg payLeg = IborLeg(paySchedule,payIndex_)
                     .withNotionals(payNotionals)
                     .withPaymentDayCounter(payIndex_->dayCounter())
                     .withPaymentAdjustment(payBDC_)
                     .withFixingDays(payIndex_->fixingDays());
        Leg recLeg = IborLeg(recSchedule,recIndex_)
                     .withNotionals(recNotionals)
                     .withPaymentDayCounter(recIndex_->dayCounter())
                     .withPaymentAdjustment(recBDC_)
                     .withFixingDays(recIndex_->fixingDays());
        // notional payment legs...
        Leg payNtlExchange, recNtlExchange;
        ext::shared_ptr<Coupon> coupon; // temp variable for easy access
        coupon = ext::dynamic_pointer_cast<Coupon>(payLeg.front());
        payNtlExchange.push_back(ext::shared_ptr<CashFlow>(new SimpleCashFlow(-payNotionals[0],coupon->accrualStartDate())));
        coupon = ext::dynamic_pointer_cast<Coupon>(payLeg.back());
        payNtlExchange.push_back(ext::shared_ptr<CashFlow>(new SimpleCashFlow(+payNotionals[0],coupon->date())));
        coupon = ext::dynamic_pointer_cast<Coupon>(recLeg.front());
        recNtlExchange.push_back(ext::shared_ptr<CashFlow>(new SimpleCashFlow(-recNotionals[0],coupon->accrualStartDate())));
        coupon = ext::dynamic_pointer_cast<Coupon>(recLeg.back());
        recNtlExchange.push_back(ext::shared_ptr<CashFlow>(new SimpleCashFlow(+recNotionals[0],coupon->date())));
        if (fxResetable) {
            // we need to calculate variable notionals and rebuild the respective legs
            // note, this approach neglegts convexity between FX and OIS-Libor-spread and FX timing adjustment
            // for details see Fujii/Shimada/Takahashi, "Collateral Posting and Choice of Collateral Currency", 2010
            // http://ssrn.com/abstract=1601866
            // we distinguish expected notional and interest payments for detailed cash flow analysis
            if (spreadOnRecLeg_) { // adjust the pay leg
                coupon = ext::dynamic_pointer_cast<Coupon>(payLeg.front());
                payNotionals[0] = 1.0/payFxForDom_ * recDisRelinkableHandle_->discount(coupon->accrualStartDate())/payDisRelinkableHandle_->discount(coupon->accrualStartDate());
                for (Size i=1; i<payNotionals.size(); ++i) {
                    payNotionals[i] = 1.0/payFxForDom_ * recDisRelinkableHandle_->discount(payLeg[i-1]->date())/payDisRelinkableHandle_->discount(payLeg[i-1]->date());
                }
                payLeg = IborLeg(paySchedule,payIndex_)
                     .withNotionals(payNotionals)
                     .withPaymentDayCounter(payIndex_->dayCounter())
                     .withPaymentAdjustment(payBDC_)
                     .withFixingDays(payIndex_->fixingDays());
                payNtlExchange.clear();
                payNtlExchange.push_back(ext::shared_ptr<CashFlow>(new SimpleCashFlow(-payNotionals[0],coupon->accrualStartDate())));
                for (Size i=1; i<payNotionals.size(); ++i) {
                    payNtlExchange.push_back(ext::shared_ptr<CashFlow>(new SimpleCashFlow(payNotionals[i-1]-payNotionals[i],payLeg[i-1]->date())));
                }
                payNtlExchange.push_back(ext::shared_ptr<CashFlow>(new SimpleCashFlow(payNotionals[payNotionals.size()-1],payLeg[payLeg.size()-1]->date())));
            } else { // adjust the receive leg
                coupon = ext::dynamic_pointer_cast<Coupon>(recLeg.front());
                recNotionals[0] = 1.0/recFxForDom_ * payDisRelinkableHandle_->discount(coupon->accrualStartDate())/recDisRelinkableHandle_->discount(coupon->accrualStartDate());
                for (Size i=1; i<recNotionals.size(); ++i) {
                    recNotionals[i] = 1.0/recFxForDom_ * payDisRelinkableHandle_->discount(recLeg[i-1]->date())/recDisRelinkableHandle_->discount(recLeg[i-1]->date());
                }
                recLeg = IborLeg(recSchedule,recIndex_)
                     .withNotionals(recNotionals)
                     .withPaymentDayCounter(recIndex_->dayCounter())
                     .withPaymentAdjustment(recBDC_)
                     .withFixingDays(recIndex_->fixingDays());
                recNtlExchange.clear();
                recNtlExchange.push_back(ext::shared_ptr<CashFlow>(new SimpleCashFlow(-recNotionals[0],coupon->accrualStartDate())));
                for (Size i=1; i<recNotionals.size(); ++i) {
                    recNtlExchange.push_back(ext::shared_ptr<CashFlow>(new SimpleCashFlow(recNotionals[i-1]-recNotionals[i],recLeg[i-1]->date())));
                }
                recNtlExchange.push_back(ext::shared_ptr<CashFlow>(new SimpleCashFlow(recNotionals[recNotionals.size()-1],recLeg[recLeg.size()-1]->date())));
            }
        }
        // now we can assemble the swap
        std::vector<Leg> legs(4);
        legs[0]  = payLeg;
        legs[1]  = recLeg;
        legs[2]  = payNtlExchange;
        legs[3]  = recNtlExchange;
        std::vector<bool> payer(4);
        payer[0] = true;
        payer[1] = false;
        payer[2] = true;
        payer[3] = false;
        basisSwap_ = ext::shared_ptr<BasisSwap>(new BasisSwap(legs, payer, (spreadOnRecLeg_ ? 1 : 0), true));
        // we set up pricing engine here coz it's associated to instrument
        std::vector<Handle<YieldTermStructure>> discCurves(4);
        discCurves[0] = payDisRelinkableHandle_;
        discCurves[1] = recDisRelinkableHandle_;
        discCurves[2] = payDisRelinkableHandle_;
        discCurves[3] = recDisRelinkableHandle_;
        std::vector<Real> fxForDom(4);
        fxForDom[0] = payFxForDom_;
        fxForDom[1] = recFxForDom_;
        fxForDom[2] = payFxForDom_;
        fxForDom[3] = recFxForDom_;
        basisSwap_->setPricingEngine(ext::shared_ptr<PricingEngine>(new 
            BasisSwapEngine(discCurves, fxForDom, false)));
    }


}


