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

#include <ql/patterns/visitor.hpp>
#include <ql/patterns/lazyobject.hpp>

#include <ql/cashflows/cashflowvectors.hpp>
#include <ql/cashflows/fxrangeaccrualfixed.hpp>
#include <ql/indexes/fxindex.hpp>
#include <ql/math/distributions/normaldistribution.hpp>
#include <ql/pricingengines/blackformula.hpp>
#include <ql/termstructures/yieldtermstructure.hpp>
#include <ql/termstructures/volatility/equityfx/blackvoltermstructure.hpp>
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
      lowerTrigger_(lowerTrigger), upperTrigger_(upperTrigger), pricer_(0), rangeAccrual_(0.0) {
        QL_REQUIRE(observationsSchedule_, "observationsSchedule_ required.");
        QL_REQUIRE(fxIndex_, "fxIndex_ required.");
        QL_REQUIRE(lowerTrigger_ > 0.0, "lowerTrigger_ > 0.0 required.");
        QL_REQUIRE(lowerTrigger_ < upperTrigger_, "lowerTrigger_ < upperTrigger_ required.");
    }


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
    lowerTrigger_(lowerTrigger), upperTrigger_(upperTrigger), pricer_(0), rangeAccrual_(0.0) {
        QL_REQUIRE(observationsSchedule_, "observationsSchedule_ required.");
        QL_REQUIRE(fxIndex_, "fxIndex_ required.");
        QL_REQUIRE(lowerTrigger_ > 0.0, "lowerTrigger_ > 0.0 required.");
        QL_REQUIRE(lowerTrigger_ < upperTrigger_, "lowerTrigger_ < upperTrigger_ required.");
    }


    void FxRangeAccrualFixedCoupon::performCalculations() const {
        FixedRateCoupon::performCalculations();
        if (pricer_) {
            pricer_->initialize(*this);
            rangeAccrual_ = pricer_->rangeAccrual();
            additionalResults_ = pricer_->additionalResults();
        } else {
            // calculate fall-back via intrinsic value
            Real inRange = 0.0;
            for (auto d : observationsSchedule()->dates()) {
                auto indexObservation = fxIndex()->fixing(d);
                if (indexObservation >= lowerTrigger() && indexObservation <= upperTrigger())
                    inRange += 1.0;
            }
            rangeAccrual_ = inRange / observationsSchedule()->dates().size();
        }
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

    void FxRangeAccrualFixedCoupon::setPricer(const ext::shared_ptr<FxRangeAccrualFixedCouponPricer>& pricer){
        if (pricer_ != nullptr)
            unregisterWith(pricer_);
        pricer_ = pricer;
        if (pricer_ != nullptr)
            registerWith(pricer_);
        update();
    }

    FxRangeAccrualFixedCouponPricer::FxRangeAccrualFixedCouponPricer(
            Handle<BlackVolTermStructure> fxVolatility
    ) : fxVolatility_(fxVolatility), rangeAccrual_(Null<Real>()) {}

    void FxRangeAccrualFixedCouponPricer::initialize(const FxRangeAccrualFixedCoupon& coupon) {
        additionalResults_.clear();
        // implement rangeAccrual_ calculation...
        Real minStd = 0.0005;  // 1% * sqrt(1d)
        CumulativeNormalDistribution Phi;
        Real daysInRange = 0.0;
        for (auto d : coupon.observationsSchedule()->dates()) {
            std::ostringstream s;
            s << io::iso_date(d);
            std::string date_s = s.str();
            Real indexObservation = coupon.fxIndex()->fixing(d);
            Real standardDevLow = 0.0;
            Real standardDevUpp = 0.0;
            if (d > fxVolatility_->referenceDate()) {
                standardDevLow = std::sqrt(std::max(fxVolatility_->blackVariance(d, coupon.lowerTrigger(), true), 0.0));
                standardDevUpp = std::sqrt(std::max(fxVolatility_->blackVariance(d, coupon.upperTrigger(), true), 0.0));
            }
            Real inRangeProbability = 0.0;
            if (standardDevLow < minStd) {  // calculate intrinsic value
                if (indexObservation >= coupon.lowerTrigger() &&
                    indexObservation <= coupon.upperTrigger())
                    inRangeProbability = 1.0;
            } else {  // put option spread valuation
                Real d2Lower = std::log(indexObservation / coupon.lowerTrigger()) / standardDevLow -
                               0.5 * standardDevLow;
                Real d2Upper = std::log(indexObservation / coupon.upperTrigger()) / standardDevUpp -
                               0.5 * standardDevUpp;
                Real putLower = Phi(-d2Lower);
                Real putUpper = Phi(-d2Upper);
                inRangeProbability = putUpper - putLower;
            }
            daysInRange += inRangeProbability;
            //
            additionalResults_["indexObservation_" + date_s] = indexObservation;
            additionalResults_["standardDevLow_" + date_s] = standardDevLow;
            additionalResults_["standardDevUpp_" + date_s] = standardDevUpp;
            additionalResults_["inRangeProbability_" + date_s] = inRangeProbability;
        }
        rangeAccrual_ = daysInRange / coupon.observationsSchedule()->dates().size();
        additionalResults_["daysInRange"] = daysInRange;
        additionalResults_["observationDays"] = (Real) coupon.observationsSchedule()->dates().size();
    }

    Real FxRangeAccrualFixedCouponPricer::rangeAccrual() const {
        return rangeAccrual_;
    }

}
