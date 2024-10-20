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
#include <ql/cashflows/cmscoupon.hpp>
#include <ql/cashflows/cmsrangeaccrualfixed.hpp>
#include <ql/cashflows/conundrumpricer.hpp>
#include <ql/indexes/swapindex.hpp>
#include <ql/math/distributions/normaldistribution.hpp>
#include <ql/termstructures/yieldtermstructure.hpp>
#include <ql/termstructures/volatility/swaption/swaptionvolstructure.hpp>
#include <ql/time/schedule.hpp>
#include <ql/time/daycounters/actual360.hpp>
#include <cmath>
#include <utility>

namespace QuantLib {

    CmsRangeAccrualFixedCoupon::CmsRangeAccrualFixedCoupon(
        // FixedRateCoupon
        const Date& paymentDate,
        Real nominal,
        Real rate,
        const DayCounter& dayCounter,
        const Date& accrualStartDate,
        const Date& accrualEndDate,
        // RA feature
        ext::shared_ptr<Schedule> observationsSchedule,
        ext::shared_ptr<SwapIndex> cmsIndex,
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
      observationsSchedule_(std::move(observationsSchedule)), cmsIndex_(std::move(cmsIndex)),
      lowerTrigger_(lowerTrigger), upperTrigger_(upperTrigger), pricer_(0), rangeAccrual_(0.0) {
        QL_REQUIRE(observationsSchedule_, "observationsSchedule_ required.");
        QL_REQUIRE(cmsIndex_, "cmsIndex_ required.");
        QL_REQUIRE(lowerTrigger_ > 0.0, "lowerTrigger_ > 0.0 required.");
        QL_REQUIRE(lowerTrigger_ < upperTrigger_, "lowerTrigger_ < upperTrigger_ required.");
    }


    CmsRangeAccrualFixedCoupon::CmsRangeAccrualFixedCoupon(
        // FixedRateCoupon
        const Date& paymentDate,
        Real nominal,
        Real rate,
        const DayCounter& dayCounter,
        const Date& accrualStartDate,
        const Date& accrualEndDate,
        // RA feature
        // calculate observation schedule from coupon
        ext::shared_ptr<SwapIndex> cmsIndex,
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
                                    .withCalendar(cmsIndex->fixingCalendar())
                                    .withConvention(Following))),
    cmsIndex_(std::move(cmsIndex)),
    lowerTrigger_(lowerTrigger), upperTrigger_(upperTrigger), pricer_(0), rangeAccrual_(0.0) {
        QL_REQUIRE(observationsSchedule_, "observationsSchedule_ required.");
        QL_REQUIRE(cmsIndex_, "cmsIndex_ required.");
        QL_REQUIRE(lowerTrigger_ < upperTrigger_, "lowerTrigger_ < upperTrigger_ required.");
    }


    void CmsRangeAccrualFixedCoupon::performCalculations() const {
        FixedRateCoupon::performCalculations();
        if (pricer_) {
            pricer_->initialize(*this);
            rangeAccrual_ = pricer_->rangeAccrual();
            additionalResults_ = pricer_->additionalResults();
        } else {
            // calculate fall-back via intrinsic value
            Real inRange = 0.0;
            for (auto d : observationsSchedule()->dates()) {
                auto indexObservation = cmsIndex()->fixing(d);
                if (indexObservation >= lowerTrigger() && indexObservation <= upperTrigger())
                    inRange += 1.0;
            }
            rangeAccrual_ = inRange / observationsSchedule()->dates().size();
        }
    }

    Real CmsRangeAccrualFixedCoupon::rangeAccrual() const {
        calculate();  // make sure member is calculated
        return rangeAccrual_;
    }


    Real CmsRangeAccrualFixedCoupon::amount() const {
        calculate();
        return rangeAccrual_ * FixedRateCoupon::amount();
    }


    void CmsRangeAccrualFixedCoupon::accept(AcyclicVisitor& v) {
        auto* v1 = dynamic_cast<Visitor<CmsRangeAccrualFixedCoupon>*>(&v);
        if (v1 != nullptr)
            v1->visit(*this);
        else
            FixedRateCoupon::accept(v);
    }

    void CmsRangeAccrualFixedCoupon::setPricer(
        const ext::shared_ptr<CmsRangeAccrualFixedCouponPricer>& pricer) {
        if (pricer_ != nullptr)
            unregisterWith(pricer_);
        pricer_ = pricer;
        if (pricer_ != nullptr)
            registerWith(pricer_);
        update();
    }

    CmsRangeAccrualFixedCouponPricer::CmsRangeAccrualFixedCouponPricer(
        Handle<SwaptionVolatilityStructure> swaptionVolatility)
    : swaptionVolatility_(swaptionVolatility), haganPricer_(),
      rangeAccrual_(Null<Real>()) {}

    CmsRangeAccrualFixedCouponPricer::CmsRangeAccrualFixedCouponPricer(
        const ext::shared_ptr<CmsCouponPricer> cmsCouponPricer)
    : swaptionVolatility_(cmsCouponPricer->swaptionVolatility()), haganPricer_(),
      rangeAccrual_(Null<Real>()) {
        // We allow a CmsCouponPricer here to enable a general interface.
        // However, for our implementation, we require a HaganPricer.
        // Consequently, we need to down-cast the pricer.
        QL_REQUIRE(cmsCouponPricer, "cmsCouponPricer is required.");
        haganPricer_ = ext::dynamic_pointer_cast<HaganPricer>(cmsCouponPricer);
        QL_REQUIRE(haganPricer_, "Cannot down-cast cmsCouponPricer to HaganPricer.");
    }

    Real CmsRangeAccrualFixedCouponPricer::cmsPutOption(
        const ext::shared_ptr<SwapIndex>& cmsIndex,
        const Date& exerciseDate,
        const Date& paymentDate,
        const Real optionStrike,
        const Real spreadWidth // 1bp
    ) {
        QL_REQUIRE(haganPricer_, "haganPricer_ required.");
        QL_REQUIRE(spreadWidth > 0.0, "spreadWidth > 0.0 required.");
        CmsCoupon cmsCoupon(
            paymentDate,
            1.0,               // nominal
            exerciseDate,      // startDate
            exerciseDate + 1,  // endDate
            0,                 // fixingDays
            cmsIndex,
            1.0,               // gearing
            0.0,               // spread
            Date(),            // refPeriodStart
            Date(),            // refPeriodEnd
            Actual360()
        );
        cmsCoupon.setPricer(haganPricer_);
        cmsCoupon.performCalculations();
        Real swapRate = cmsIndex->fixing(exerciseDate);
        Real putSprd = 0.0;
        if (optionStrike > swapRate) {  // calculate call spread to improve numerical stability
            Real calPlus = haganPricer_->capletRate(optionStrike + 0.5 * spreadWidth);
            Real calMins = haganPricer_->capletRate(optionStrike - 0.5 * spreadWidth);
            putSprd = 1.0 - (calMins - calPlus) / spreadWidth;
        } else {
            Real putPlus = haganPricer_->floorletRate(optionStrike + 0.5 * spreadWidth);
            Real putMins = haganPricer_->floorletRate(optionStrike - 0.5 * spreadWidth);
            putSprd = (putPlus - putMins) / spreadWidth;
        }
        return putSprd;
    }


    void CmsRangeAccrualFixedCouponPricer::initialize(const CmsRangeAccrualFixedCoupon& coupon) {
        additionalResults_.clear();
        Period swapTerm = coupon.cmsIndex()->tenor();
        // implement rangeAccrual_ calculation...
        Real minStd = 0.000005;  // 1bp * sqrt(1d)
        Real strikeLow = coupon.lowerTrigger();
        Real strikeUpp = coupon.upperTrigger();
        CumulativeNormalDistribution Phi;
        Real daysInRange = 0.0;
        for (auto d : coupon.observationsSchedule()->dates()) {
            std::ostringstream s;
            s << io::iso_date(d);
            std::string date_s = s.str();
            Real indexObservation = coupon.cmsIndex()->fixing(d);
            Real standardDevLow = 0.0;
            Real standardDevUpp = 0.0;
            //
            Real putLow = 0.0;
            Real putUpp = 0.0;
            if (d > swaptionVolatility_->referenceDate()) {
                standardDevLow = std::sqrt(std::max(swaptionVolatility_->blackVariance(d, swapTerm, coupon.lowerTrigger(), true), 0.0));
                standardDevUpp = std::sqrt(std::max(swaptionVolatility_->blackVariance(d, swapTerm, coupon.upperTrigger(), true), 0.0));
            }
            Real inRangeProbability = 0.0;
            //
            if (standardDevLow < minStd) { // calculate intrinsic value
                putLow = (indexObservation < strikeLow) ? (1.0) : (0.0);
            } else {
                if (haganPricer_) { // replication
                    putLow = cmsPutOption(coupon.cmsIndex(), d, coupon.date(), strikeLow);
                } else { // fall-back to Bachelier w/o CMS adjustment
                    putLow = Phi((strikeLow - indexObservation) / standardDevLow);
                }
            }
            //
            if (standardDevUpp < minStd) { // calculate intrinsic value
                putUpp = (indexObservation < strikeUpp) ? (1.0) : (0.0);
            } else {
                if (haganPricer_) { // replication
                    putUpp = cmsPutOption(coupon.cmsIndex(), d, coupon.date(), strikeUpp);
                } else { // fall-back to Bachelier w/o CMS adjustment
                    putUpp = Phi((strikeUpp - indexObservation) / standardDevUpp);
                }
            }
            inRangeProbability = putUpp - putLow;
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

    Real CmsRangeAccrualFixedCouponPricer::rangeAccrual() const {
        return rangeAccrual_;
    }

}
