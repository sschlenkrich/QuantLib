/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2024 Sebastian Schlenkrich

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

#include "toplevelfixture.hpp"
#include "utilities.hpp"

#include <ql/types.hpp>
#include <ql/compounding.hpp>
#include <ql/time/all.hpp>
#include <ql/termstructures/yield/all.hpp>
#include <ql/termstructures/volatility/all.hpp>

#include <ql/indexes/fxindex.hpp>
#include <ql/cashflows/fxrangeaccrualfixed.hpp>


using namespace QuantLib;
using namespace boost::unit_test_framework;

BOOST_FIXTURE_TEST_SUITE(QuantLibTests, TopLevelFixture)

BOOST_AUTO_TEST_SUITE(FxRangeAccrualTests)

namespace {
    
    Period termsData[] = {
        Period( 0,Days),
        Period( 1,Years),
        Period( 2,Years),
        Period( 3,Years),
        Period( 5,Years),
        Period( 7,Years),
        Period(10,Years),
        Period(15,Years),
        Period(20,Years),
        Period(61,Years)   // avoid extrapolation issues with 30y caplets
    };
    std::vector<Period> terms(termsData, termsData+10);

    Real domDiscRatesData[] = {
        0.0300,
        0.0300,
        0.0300,
        0.0300,
        0.0300,
        0.0300,
        0.0300,
        0.0300,
        0.0300,
        0.0300
    };
    std::vector<Real> domDiscRates(domDiscRatesData, domDiscRatesData + 10);

    Real forDiscRatesData[] = {
        0.0400,
        0.0400,
        0.0400,
        0.0400,
        0.0400,
        0.0400,
        0.0400,
        0.0400,
        0.0400,
        0.0400
    };
    std::vector<Real> forDiscRates(forDiscRatesData, forDiscRatesData + 10);

    Handle<YieldTermStructure> getYtsh(const std::vector<Period>& terms,
                                      const std::vector<Real>& rates,
                                      const Real spread = 0.0) {
        Date today = Settings::instance().evaluationDate();
        std::vector<Date> dates;
        dates.reserve(terms.size());
        for (auto term : terms)
            dates.push_back(NullCalendar().advance(today, term, Unadjusted));
        std::vector<Real> ratesPlusSpread(rates);
        for (double& k : ratesPlusSpread)
            k += spread;
        ext::shared_ptr<YieldTermStructure> ts =
            ext::shared_ptr<YieldTermStructure>(new InterpolatedZeroCurve<Cubic>(
                dates, ratesPlusSpread, Actual365Fixed(), NullCalendar()));
        return RelinkableHandle<YieldTermStructure>(ts);
    }


} // namespace


//******************************************************************************************//
//******************************************************************************************//

BOOST_AUTO_TEST_CASE(testCouponSetup)  {

    BOOST_TEST_MESSAGE("Testing FX range accrual coupon without pricer...");

    Date today = Settings::instance().evaluationDate();

    BOOST_TEST_MESSAGE("Today: " << today);

    Handle<YieldTermStructure> domYtsh = getYtsh(terms, domDiscRates);
    Handle<YieldTermStructure> forYtsh = getYtsh(terms, forDiscRates);

    BOOST_TEST_MESSAGE("domYtsh at today: " << domYtsh->discount(today));
    BOOST_TEST_MESSAGE("forYtsh at today: " << forYtsh->discount(today));

    RelinkableHandle<Quote> spotHandle; // = RelinkableHandle<Quote>(ext::make_shared<SimpleQuote>(0.0));

    ext::shared_ptr<FxIndex> fxIndex =
        ext::make_shared<FxIndex>(
            "EUR-USD", TARGET(), domYtsh, forYtsh, spotHandle);

    Schedule fixingSchedule =
        Schedule(Date(1, January, 2015), Date(31, December, 2015), Period(1, Days),
                              TARGET(), Following, Following, DateGeneration::Forward, false);

    Real fx = 1.0;
    for (auto date : fixingSchedule.dates()) {
        fxIndex->addFixing(date, fx);
        fx += 0.001;
    }

    auto startDate = Date(31, August, 2015);
    auto endDate = Date(30, September, 2015);
    auto payDate = Date(30, September, 2015);
    Real notional = 100.0;
    Real fixedRate = 0.01;
    auto dc = Actual360();

    auto make_coupon = [payDate, notional, fixedRate, dc, startDate, endDate,
                        fxIndex](Real lowerTrigger, Real upperTrigger) {
        return FxRangeAccrualFixedCoupon(payDate, notional, fixedRate, dc, startDate, endDate,
                                         fxIndex, lowerTrigger, upperTrigger);
    };

    auto coupon = make_coupon(1.10, 1.175);

    for (auto date : coupon.observationsSchedule()->dates()) {
        BOOST_TEST_MESSAGE("ObsDate: " << date << ", Index: " << fxIndex->fixing(date));
    }

    std::vector<FxRangeAccrualFixedCoupon> cps;
    cps.push_back(make_coupon(1.10, 1.15));
    cps.push_back(make_coupon(1.15, 1.20));
    cps.push_back(make_coupon(1.20, 1.22));
    cps.push_back(make_coupon(1.15, 1.175));  // RA 8/23 = 0.347826087

    for (auto cp : cps) {
        BOOST_TEST_MESSAGE("Rate: " << cp.rate() << ", RA: " << cp.rangeAccrual() << ", Amount: " << cp.amount());
    }
}


BOOST_AUTO_TEST_CASE(testCouponPricing)  {

    BOOST_TEST_MESSAGE("Testing FX range accrual coupon with pricer...");

    Date today = Settings::instance().evaluationDate();

    BOOST_TEST_MESSAGE("Today: " << today);

    Handle<YieldTermStructure> domYtsh = getYtsh(terms, domDiscRates);
    Handle<YieldTermStructure> forYtsh = getYtsh(terms, forDiscRates);

    BOOST_TEST_MESSAGE("domYtsh at today: " << domYtsh->discount(today));
    BOOST_TEST_MESSAGE("forYtsh at today: " << forYtsh->discount(today));

    RelinkableHandle<Quote> spotHandle; // = RelinkableHandle<Quote>(ext::make_shared<SimpleQuote>(0.0));

    ext::shared_ptr<FxIndex> fxIndex =
        ext::make_shared<FxIndex>(
            "EUR-USD", TARGET(), domYtsh, forYtsh, spotHandle);

    Schedule fixingSchedule =
        Schedule(Date(1, January, 2015), Date(31, December, 2015), Period(1, Days),
                              TARGET(), Following, Following, DateGeneration::Forward, false);

    Real fx = 1.0;
    for (auto date : fixingSchedule.dates()) {
        fxIndex->addFixing(date, fx);
        fx += 0.001;
    }


    RelinkableHandle<Quote> volQuote = RelinkableHandle<Quote>(ext::make_shared<SimpleQuote>(0.25));
    ext::shared_ptr<BlackVolTermStructure> volTs = ext::shared_ptr<BlackVolTermStructure>(
        new BlackConstantVol(today, TARGET(), volQuote, Actual365Fixed())
    );
    RelinkableHandle<BlackVolTermStructure> volTsh(volTs);
    ext::shared_ptr<FxRangeAccrualFixedCouponPricer> pricer =
        ext::make_shared<FxRangeAccrualFixedCouponPricer>(volTsh);


    auto startDate = Date(31, August, 2015);
    auto endDate = Date(30, September, 2015);
    auto payDate = Date(30, September, 2015);
    Real notional = 100.0;
    Real fixedRate = 0.01;
    auto dc = Actual360();

    auto make_coupon = [payDate, notional, fixedRate, dc, startDate, endDate,
                        fxIndex](Real lowerTrigger, Real upperTrigger) {
        return FxRangeAccrualFixedCoupon(payDate, notional, fixedRate, dc, startDate, endDate,
                                         fxIndex, lowerTrigger, upperTrigger);
    };

    auto coupon = make_coupon(1.10, 1.175);

    for (auto date : coupon.observationsSchedule()->dates()) {
        BOOST_TEST_MESSAGE("ObsDate: " << date << ", Index: " << fxIndex->fixing(date));
    }

    std::vector<FxRangeAccrualFixedCoupon> cps;
    cps.push_back(make_coupon(1.10, 1.15));
    cps.push_back(make_coupon(1.15, 1.20));
    cps.push_back(make_coupon(1.20, 1.22));
    cps.push_back(make_coupon(1.15, 1.175));  // RA 8/23 = 0.347826087

    for (FxRangeAccrualFixedCoupon& cp : cps) { // make sure we get a reference and not a copy
        cp.setPricer(pricer);
    }

    for (auto cp : cps) {
        BOOST_TEST_MESSAGE("Rate: " << cp.rate() << ", RA: " << cp.rangeAccrual() << ", Amount: " << cp.amount());
    }
}



BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE_END()
