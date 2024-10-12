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

#include <ql/indexes/all.hpp>
#include <ql/cashflows/cmsrangeaccrualfixed.hpp>

#include <iostream>

using namespace QuantLib;
using namespace boost::unit_test_framework;

BOOST_FIXTURE_TEST_SUITE(QuantLibTests, TopLevelFixture)

BOOST_AUTO_TEST_SUITE(CmsRangeAccrualTests)

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

    Real discRatesData[] = {
        0.0250,
        0.0250,
        0.0250,
        0.0250,
        0.0250,
        0.0250,
        0.0250,
        0.0250,
        0.0250,
        0.0250
    };
    std::vector<Real> discRates(discRatesData, discRatesData + 10);

    Real projRatesData[] = {
        0.0280,
        0.0280,
        0.0280,
        0.0280,
        0.0280,
        0.0280,
        0.0280,
        0.0300,
        0.0400,
        0.0400
    };
    std::vector<Real> projRates(projRatesData, projRatesData + 10);

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

    BOOST_TEST_MESSAGE("Testing CMS range accrual coupon without pricer...");

    Date today = Settings::instance().evaluationDate();

    BOOST_TEST_MESSAGE("Today: " << today);

    Handle<YieldTermStructure> discYtsh = getYtsh(terms, discRates);
    Handle<YieldTermStructure> projYtsh = getYtsh(terms, projRates);

    BOOST_TEST_MESSAGE("discYtsh at today: " << discYtsh->discount(today));
    BOOST_TEST_MESSAGE("projYtsh at today: " << projYtsh->discount(today));

    ext::shared_ptr<SwapIndex> swapIndex =
        ext::make_shared<EuriborSwapIsdaFixA>(Period(10, Years), projYtsh, discYtsh);

    Schedule fixingSchedule =
        Schedule(Date(1, January, 2015), Date(31, December, 2015), Period(1, Days),
                              TARGET(), Following, Following, DateGeneration::Forward, false);

    Real swapRate = 0.0100;
    for (auto date : fixingSchedule.dates()) {
        swapIndex->addFixing(date, swapRate);
        swapRate += 0.0001;
    }

    auto startDate = Date(31, August, 2015);
    auto endDate = Date(30, September, 2015);
    auto payDate = Date(30, September, 2015);
    Real notional = 100.0;
    Real fixedRate = 0.01;
    auto dc = Actual360();

    auto make_coupon = [payDate, notional, fixedRate, dc, startDate, endDate,
                        swapIndex](Real lowerTrigger, Real upperTrigger) {
        return CmsRangeAccrualFixedCoupon(payDate, notional, fixedRate, dc, startDate, endDate,
                                         swapIndex, lowerTrigger, upperTrigger);
    };

    auto coupon = make_coupon(0.0100, 0.0300);

    for (auto date : coupon.observationsSchedule()->dates()) {
        BOOST_TEST_MESSAGE("ObsDate: " << date << ", Index: " << swapIndex->fixing(date));
    }

    std::vector<CmsRangeAccrualFixedCoupon> cps;
    cps.push_back(make_coupon(0.0250, 0.0260));  // RA 0
    cps.push_back(make_coupon(0.0260, 0.0275));  // RA 8/23
    cps.push_back(make_coupon(0.0275, 0.0280));  // RA 5/23
    cps.push_back(make_coupon(0.0280, 0.0290));  // RA 10/23
    cps.push_back(make_coupon(0.0290, 0.0300));  // RA 0
    
    for (auto cp : cps) {
        BOOST_TEST_MESSAGE("Rate: " << cp.rate() << ", RA: " << cp.rangeAccrual() << ", Amount: " << cp.amount());
    }
}


BOOST_AUTO_TEST_CASE(testCouponPricing)  {

    BOOST_TEST_MESSAGE("Testing CMS range accrual coupon with pricer...");

    // std::cout << "Press Enter to Continue";  // use interrupt to attach debugger
    // std::cin.ignore();

    Date today = Settings::instance().evaluationDate();

    BOOST_TEST_MESSAGE("Today: " << today);

    Handle<YieldTermStructure> discYtsh = getYtsh(terms, discRates);
    Handle<YieldTermStructure> projYtsh = getYtsh(terms, projRates);

    BOOST_TEST_MESSAGE("discYtsh at today: " << discYtsh->discount(today));
    BOOST_TEST_MESSAGE("projYtsh at today: " << projYtsh->discount(today));


    ext::shared_ptr<SwapIndex> swapIndex =
        ext::make_shared<EuriborSwapIsdaFixA>(Period(10, Years), projYtsh, discYtsh);

    Schedule fixingSchedule =
        Schedule(Date(1, January, 2015), Date(31, December, 2015), Period(1, Days),
                              TARGET(), Following, Following, DateGeneration::Forward, false);

    Real swapRate = 0.0100;
    for (auto date : fixingSchedule.dates()) {
        swapIndex->addFixing(date, swapRate);
        swapRate += 0.0001;
    }


    RelinkableHandle<Quote> volQuote = RelinkableHandle<Quote>(ext::make_shared<SimpleQuote>(0.0050));
    VolatilityType vType = Normal;
    //
    ext::shared_ptr<SwaptionVolatilityStructure> volTs =
        ext::shared_ptr<ConstantSwaptionVolatility>(
            new ConstantSwaptionVolatility(
                today, TARGET(), Following, volQuote, Actual365Fixed(), vType)
    );
    RelinkableHandle<SwaptionVolatilityStructure> volTsh(volTs);
    //
    ext::shared_ptr<CmsRangeAccrualFixedCouponPricer> pricer =
        ext::make_shared<CmsRangeAccrualFixedCouponPricer>(volTsh);

    auto startDate = Date(31, August, 2015);
    auto endDate = Date(30, September, 2015);
    auto payDate = Date(30, September, 2015);
    Real notional = 100.0;
    Real fixedRate = 0.01;
    auto dc = Actual360();
    
    auto make_coupon = [payDate, notional, fixedRate, dc, startDate, endDate,
                        swapIndex](Real lowerTrigger, Real upperTrigger) {
        return CmsRangeAccrualFixedCoupon(payDate, notional, fixedRate, dc, startDate, endDate,
                                         swapIndex, lowerTrigger, upperTrigger);
    };
    
    auto coupon = make_coupon(0.0100, 0.0300);
    
    for (Date date : coupon.observationsSchedule()->dates()) {
        BOOST_TEST_MESSAGE("ObsDate: " << io::iso_date(date) << ", Index: " << swapIndex->fixing(date));
    }
    
    std::vector<CmsRangeAccrualFixedCoupon> cps;
    cps.push_back(make_coupon(0.0000, 0.0250)); // RA 0
    cps.push_back(make_coupon(0.0250, 0.0260)); // RA 0
    cps.push_back(make_coupon(0.0260, 0.0275)); // RA 8/23
    cps.push_back(make_coupon(0.0275, 0.0280)); // RA 5/23
    cps.push_back(make_coupon(0.0280, 0.0290)); // RA 10/23
    cps.push_back(make_coupon(0.0290, 0.0300)); // RA 0
    cps.push_back(make_coupon(0.0300, 0.0600)); // RA 0
    
    for (CmsRangeAccrualFixedCoupon& cp : cps) { // make sure we get a reference and not a copy
        cp.setPricer(pricer);
    }
    
    BOOST_TEST_MESSAGE("Coupon Results:");
    for (CmsRangeAccrualFixedCoupon& cp : cps) {
        BOOST_TEST_MESSAGE("Rate: " << cp.rate() << ", RA: " << cp.rangeAccrual() << ", Amount: " << cp.amount());
    }
    
    BOOST_TEST_MESSAGE("Additional Results 5th coupon:");
    auto additionalResults = cps[4].additionalResults();
    for (auto const& x : additionalResults) {
        std::cout << x.first << " : " << x.second << std::endl;
    }
}



BOOST_AUTO_TEST_CASE(testCouponLeg) {

    BOOST_TEST_MESSAGE("Testing CMS range accrual coupon leg with pricer...");

    // std::cout << "Press Enter to Continue";  // use interrupt to attach debugger
    // std::cin.ignore();

    Date today = Settings::instance().evaluationDate();

    BOOST_TEST_MESSAGE("Today: " << today);

    Handle<YieldTermStructure> discYtsh = getYtsh(terms, discRates);
    Handle<YieldTermStructure> projYtsh = getYtsh(terms, projRates);

    BOOST_TEST_MESSAGE("discYtsh at today: " << discYtsh->discount(today));
    BOOST_TEST_MESSAGE("projYtsh at today: " << projYtsh->discount(today));

    ext::shared_ptr<SwapIndex> swapIndex =
        ext::make_shared<EuriborSwapIsdaFixA>(Period(10, Years), projYtsh, discYtsh);

    Schedule fixingSchedule =
        Schedule(Date(1, January, 2015), Date(31, December, 2015), Period(1, Days), TARGET(),
                 Following, Following, DateGeneration::Forward, false);

    Real swapRate = 0.0100;
    for (auto date : fixingSchedule.dates()) {
        swapIndex->addFixing(date, swapRate);
        swapRate += 0.0001;
    }

    RelinkableHandle<Quote> volQuote =
        RelinkableHandle<Quote>(ext::make_shared<SimpleQuote>(0.0050));
    VolatilityType vType = Normal;
    //
    ext::shared_ptr<SwaptionVolatilityStructure> volTs =
        ext::shared_ptr<ConstantSwaptionVolatility>(new ConstantSwaptionVolatility(
            today, TARGET(), Following, volQuote, Actual365Fixed(), vType));
    RelinkableHandle<SwaptionVolatilityStructure> volTsh(volTs);
    //
    ext::shared_ptr<CmsRangeAccrualFixedCouponPricer> pricer =
        ext::make_shared<CmsRangeAccrualFixedCouponPricer>(volTsh);

    auto startDate = Date(15, January, 2015);
    auto endDate = Date(15, January, 2045);
    auto period = Period(3, Months);
    auto cal = TARGET();
    auto endOfMonth = false;
    Schedule couponSchedule = Schedule(startDate, endDate, period, cal, Following, Following,
                                       DateGeneration::Backward, endOfMonth);

    Real notional = 100.0;
    Real fixedRate = 0.01;
    auto dc = Actual360();

    auto lowerTrigger = 0.0250;
    auto upperTrigger = 0.0350;

    std::vector<CmsRangeAccrualFixedCoupon> leg;
    for (Size k = 0; k < couponSchedule.dates().size()-1; ++k) {
        auto startDate = couponSchedule.dates()[k];
        auto endDate = couponSchedule.dates()[k+1];
        auto payDate = couponSchedule.dates()[k + 1];
        leg.push_back(CmsRangeAccrualFixedCoupon(payDate, notional, fixedRate, dc, startDate,
                                                endDate, swapIndex, lowerTrigger, upperTrigger)
        );
    }

    for (CmsRangeAccrualFixedCoupon& cp : leg) { // make sure we get a reference and not a copy
        cp.setPricer(pricer);
    }

    BOOST_TEST_MESSAGE("Coupon Results:");
    for (CmsRangeAccrualFixedCoupon& cp : leg) {
        BOOST_TEST_MESSAGE("Start: " << io::iso_date(cp.accrualStartDate())
                                     << ", End: " << io::iso_date(cp.accrualEndDate())
                                     << ", RA: " << cp.rangeAccrual()
                                     << ", Rate: " << cp.rate()
                                    << ", Amount: " << cp.amount());
    }
}



BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE_END()
