/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2016 Andre Miemiec

 This file is part of QuantLib, a free-software/open-source library
 for financial quantitative analysts and developers - http://quantlib.org/

 QuantLib is free software: you can redistribute it and/or modify it
 under the terms of the QuantLib license.  You should have received a
 copy of the license along with this program; if not, please email
 <quantlib-dev@lists.sf.net>. The license is also available online at
 <http://quantlib.org/license.shtml>.

 This program is distributed in the hope that it will be useful, but WITHOUT
 ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 FOR A PARTICULAR PURPOSE.  See the license for more details.
*/

#include <ql/experimental/basismodels/fxfwdratehelper.hpp>
#include <ql/quote.hpp>
#include <ql/currency.hpp>
#include <ql/cashflows/cashflows.hpp>


namespace QuantLib {

    namespace FxFwdRateHelperDummyNamespace {
        void no_deletion(YieldTermStructure*) {}
    }


    FxFwdRateHelper::FxFwdRateHelper( Currency baseCurrency,
                                      Currency counterCurrency,
                                      Rate   fxSpot,
                                      Natural spotLag,
                                      Calendar spotLagCal,
                                      BusinessDayConvention spotLagConv,
                                      Period swapTerm,
                                      Handle<Quote> points,
                                      Real   unit,
                                      Handle<YieldTermStructure> baseCcyDiscTermStructureHandle,
                                      Handle<YieldTermStructure> cntrCcyDiscTermStructureHandle,
                                      BootstrapType bootstrapBaseOrCounter)
      : RelativeDateRateHelper(points),       
      baseCurrency_(baseCurrency), counterCurrency_(counterCurrency), fxSpot_(fxSpot),  
      spotLag_(spotLag), spotLagCalendar_(spotLagCal), spotLagConvention_(spotLagConv), 
      swapTerm_(swapTerm), unit_(unit),
      baseCcyDiscTermStructureHandle_(baseCcyDiscTermStructureHandle), cntrCcyDiscTermStructureHandle_(cntrCcyDiscTermStructureHandle), 
      bootstrapBaseOrCounter_(bootstrapBaseOrCounter)
    {
        QL_REQUIRE( (baseCurrency_ != counterCurrency_),"Currencies should differ!");

        registerWith(baseCcyDiscTermStructureHandle_);
        registerWith(cntrCcyDiscTermStructureHandle_);

        initializeDates();
    }

    void FxFwdRateHelper::initializeDates() {

        Date today = Settings::instance().evaluationDate();

        Date effectiveDate = spotLagCalendar_.advance(today,spotLag_,Days,spotLagConvention_,false);

        Date terminationDate = effectiveDate + swapTerm_;

        earliestDate_ = effectiveDate;
        latestDate_   = spotLagCalendar_.adjust(terminationDate,spotLagConvention_);

    }

    void FxFwdRateHelper::setTermStructure(YieldTermStructure* t) {
        // do not set the relinkable handle as an observer -
        // force recalculation when needed
        bool observer = false;

        ext::shared_ptr<YieldTermStructure> temp(t, FxFwdRateHelperDummyNamespace::no_deletion);

        if(bootstrapBaseOrCounter_ == FxFwdRateHelper::Base){
           cntrCcyDiscRelinkableHandle_.linkTo(*cntrCcyDiscTermStructureHandle_, observer);
           baseCcyDiscRelinkableHandle_.linkTo(temp, observer);
        } else {
           cntrCcyDiscRelinkableHandle_.linkTo(temp, observer);
           baseCcyDiscRelinkableHandle_.linkTo(*baseCcyDiscTermStructureHandle_, observer);
        }

        RelativeDateRateHelper::setTermStructure(t);
    }


    Real FxFwdRateHelper::impliedQuote() const {

        //discount factors at maturity 
        Real baseCcyDiscountT_ = baseCcyDiscRelinkableHandle_->discount(latestDate_);
        Real cntrCcyDiscountT_ = cntrCcyDiscRelinkableHandle_->discount(latestDate_);

        //discount factors at spot 
        Real baseCcyDiscountS_ = baseCcyDiscRelinkableHandle_->discount(earliestDate_);
        Real cntrCcyDiscountS_ = cntrCcyDiscRelinkableHandle_->discount(earliestDate_);

        Real impicitQuote  = 0.0;

        impicitQuote  = unit_*fxSpot_*(baseCcyDiscountT_/cntrCcyDiscountT_*cntrCcyDiscountS_/baseCcyDiscountS_ - 1.0);

        return impicitQuote;
         
    }

    void FxFwdRateHelper::accept(AcyclicVisitor& v) {
        Visitor<FxFwdRateHelper>* v1 =
            dynamic_cast<Visitor<FxFwdRateHelper>*>(&v);
        if (v1 != 0)
            v1->visit(*this);
        else
            RateHelper::accept(v);
    }

}
