/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2014 Sebastian Schlenkrich

 */

/*! \file basisswap.hpp
    \brief (Cross currency) Interest rate swap consisting of several legs
*/

#ifndef quantlib_liborbasisswapratehelpers_hpp
#define quantlib_liborbasisswapratehelpers_hpp

#include <ql/termstructures/bootstraphelper.hpp>
#include <ql/termstructures/yield/ratehelpers.hpp>
#include <ql/instruments/vanillaswap.hpp>
#include <ql/time/calendar.hpp>
#include <ql/time/daycounter.hpp>
#include <ql/cashflows/coupon.hpp>
#include <ql/cashflows/simplecashflow.hpp>
#include <ql/cashflows/iborcoupon.hpp>
#include <ql/pricingengines/swap/discountingswapengine.hpp>

#include <ql/experimental/basismodels/basisswap.hpp>
#include <ql/experimental/basismodels/basisswapengine.hpp>

namespace QuantLib {

    // interface for RateHelpers basd on BasisSwap
    class BasisSwapRateHelper : public RelativeDateRateHelper {
    protected:
        mutable ext::shared_ptr<BasisSwap> basisSwap_; // reference instrument
    public:
        BasisSwapRateHelper( Real rate ) : RelativeDateRateHelper(rate) { }
        const ext::shared_ptr<BasisSwap> basisSwap() { return basisSwap_; }
    };

    // assume a tenor swap pay xM IborLeg vs. receive yM IborLeg
    // spread is quoted on pay or receive leg
    // optionally provide discount curve and pay/receive forward curve
    class TenorSwapRateHelper : public BasisSwapRateHelper {
    protected:
        // swap details
        Period                 tenor_;
        Period                 fwdStart_;
        Calendar               paymentCalendar_;
        BusinessDayConvention  paymentBDC_;
        bool                   spreadOnRecLeg_;
        ext::shared_ptr<IborIndex>  payIndex_;
        ext::shared_ptr<IborIndex>  recIndex_;
        // discount curve, forward curves are taken from Ibor indices
        RelinkableHandle<YieldTermStructure> discountRelinkableHandle_;
        RelinkableHandle<YieldTermStructure> payRelinkableHandle_;
        RelinkableHandle<YieldTermStructure> recRelinkableHandle_;
        // define where to use the bootstrapping curve
        bool useForDiscount_, useForPayForward_, useForRecForward_;
        // RelativeDateBootstrapHelper interface
        void initializeDates();
    public:
        // BootstrapHelper interface
        virtual Real impliedQuote() const;
        virtual void setTermStructure(YieldTermStructure* ts);
        TenorSwapRateHelper( Real                          rate,
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
                             const Handle<YieldTermStructure>& discountCurve );
    };

    class XCCYSwapRateHelper : public BasisSwapRateHelper {
    protected:
        // swap details
        Period                 tenor_;
        Period                 fwdStart_;
        Calendar               spotStartCalendar_;
        BusinessDayConvention  payBDC_;
        BusinessDayConvention  recBDC_;
        bool                   spreadOnRecLeg_;
        ext::shared_ptr<IborIndex>  payIndex_;
        ext::shared_ptr<IborIndex>  recIndex_;
        // discount curve, forward curves are taken from Ibor indices
        RelinkableHandle<YieldTermStructure> payDisRelinkableHandle_;
        RelinkableHandle<YieldTermStructure> recDisRelinkableHandle_;
        RelinkableHandle<YieldTermStructure> payForRelinkableHandle_;
        RelinkableHandle<YieldTermStructure> recForRelinkableHandle_;
        // today's FX rates (either pay or receive should be 1 unless valued in third currency)
        Real payFxForDom_, recFxForDom_;
        // cross currency swaps may be resetable on the non-spreaded leg
        bool fxResetable_;
        // define where to use the bootstrapping curve
        bool useForPayDiscount_, useForRecDiscount_, useForPayForward_, useForRecForward_;
        // encapsulate swap set up
        void setUpBasisSwap(bool fxResetable=false) const;
        // RelativeDateBootstrapHelper interface
        void initializeDates();
    public:
        // BootstrapHelper interface
        virtual Real impliedQuote() const;
        virtual void setTermStructure(YieldTermStructure* ts);
        XCCYSwapRateHelper(  Real                          rate,
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
                             bool                          fxResetable = false
                             );
    };


}

#endif
