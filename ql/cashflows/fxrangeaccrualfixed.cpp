/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2006, 2007 Giorgio Facchinetti
 Copyright (C) 2006, 2007 Mario Pucci

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


#include <ql/cashflows/cashflowvectors.hpp>
#include <ql/cashflows/fxrangeaccrualfixed.hpp>
#include <ql/indexes/fxindex.hpp>
#include <ql/math/distributions/normaldistribution.hpp>
#include <ql/pricingengines/blackformula.hpp>
#include <ql/termstructures/yieldtermstructure.hpp>
#include <ql/time/schedule.hpp>
#include <cmath>
#include <utility>

namespace QuantLib {

    FxRangeAccrualFixedCoupon::FxRangeAccrualFixedCoupon(
        // FixedRateCoupon
        const Date& paymentDate,
        Real nominal,
        Real rate,
        const DayCounter& dayCounter,
        const Date& accrualStartDate,
        const Date& accrualEndDate,
        // RA feature
        ext::shared_ptr<Schedule> observationsSchedule,
        ext::shared_ptr<FxIndex> fxIndex,
        Real lowerTrigger,
        Real upperTrigger,
        // optional FixedRateCoupon
        const Date& refPeriodStart,
        const Date& refPeriodEnd,
        const Date& exCouponDate)
    : FixedRateCoupon(paymentDate,
                      nominal,
                      rate,
                      dayCounter,
                      accrualStartDate,
                      accrualEndDate,
                      refPeriodStart,
                      refPeriodEnd,
                      exCouponDate),
      observationsSchedule_(std::move(observationsSchedule)), fxIndex_(std::move(fxIndex)),
      lowerTrigger_(lowerTrigger), upperTrigger_(upperTrigger), rangeAccrual_(0.0) {} 


    FxRangeAccrualFixedCoupon::FxRangeAccrualFixedCoupon(
        // FixedRateCoupon
        const Date& paymentDate,
        Real nominal,
        Real rate,
        const DayCounter& dayCounter,
        const Date& accrualStartDate,
        const Date& accrualEndDate,
        // RA feature
        // calculate observation schedule from coupon
        ext::shared_ptr<FxIndex> fxIndex,
        Real lowerTrigger,
        Real upperTrigger,
        // optional FixedRateCoupon
        const Date& refPeriodStart,
        const Date& refPeriodEnd,
        const Date& exCouponDate)
    : FixedRateCoupon(paymentDate,
                      nominal,
                      rate,
                      dayCounter,
                      accrualStartDate,
                      accrualEndDate,
                      refPeriodStart,
                      refPeriodEnd,
                      exCouponDate),
    observationsSchedule_(ext::make_shared<Schedule>(MakeSchedule()
                                    .from(accrualStartDate)
                                    .to(accrualEndDate)
                                    .withFrequency(Daily)
                                    .withCalendar(fxIndex->fixingCalendar())
                                    .withConvention(Following))),
    fxIndex_(std::move(fxIndex)),
    lowerTrigger_(lowerTrigger), upperTrigger_(upperTrigger), rangeAccrual_(0.0) {}


    void FxRangeAccrualFixedCoupon::performCalculations() const {
        FixedRateCoupon::performCalculations();
        // more stuff
        Real inRange = 0.0;
        for (auto d : observationsSchedule()->dates()) {
            auto indexObservation = fxIndex()->fixing(d);
            if (indexObservation >= lowerTrigger() && indexObservation <= upperTrigger())
                inRange += 1.0;
        }
        rangeAccrual_ = inRange / observationsSchedule()->dates().size();
    }

    Real FxRangeAccrualFixedCoupon::rangeAccrual() const {
        calculate();  // make sure member is calculated
        return rangeAccrual_;
    }


    Real FxRangeAccrualFixedCoupon::amount() const {
        calculate();
        return rangeAccrual_ * FixedRateCoupon::amount();
    }


    void FxRangeAccrualFixedCoupon::accept(AcyclicVisitor& v) {
        auto* v1 = dynamic_cast<Visitor<FxRangeAccrualFixedCoupon>*>(&v);
        if (v1 != nullptr)
            v1->visit(*this);
        else
            FixedRateCoupon::accept(v);
    }


}
