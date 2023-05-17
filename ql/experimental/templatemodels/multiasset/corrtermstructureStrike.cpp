/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2017 Cord Harms

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

#include <ql/experimental/templatemodels/multiasset/corrtermstructureStrike.hpp>

namespace QuantLib {

	CorrelationTermStructureStrike::CorrelationTermStructureStrike(BusinessDayConvention bdc,
                                                     const DayCounter& dc)
    : TermStructure(dc), bdc_(bdc) {}

	CorrelationTermStructureStrike::CorrelationTermStructureStrike(const Date& referenceDate,
                                                     const Calendar& cal,
                                                     BusinessDayConvention bdc,
                                                     const DayCounter& dc)
	: TermStructure(referenceDate, cal, dc), bdc_(bdc) {}

	CorrelationTermStructureStrike::CorrelationTermStructureStrike(Natural settlementDays,
                                                     const Calendar& cal,
                                                     BusinessDayConvention bdc,
                                                     const DayCounter& dc)
	: TermStructure(settlementDays, cal, dc), bdc_(bdc) {}

    void CorrelationTermStructureStrike::checkStrike(Rate k,
											  Natural ulId,
                                              bool extrapolate) const {
        QL_REQUIRE(extrapolate || allowsExtrapolation() ||
                   (k >= minStrike(ulId) && k <= maxStrike(ulId)),
                   "strike (" << k << ") is outside the curve domain ["
                   << minStrike(ulId) << "," << maxStrike(ulId)<< "] for underlying " << ulId);
    }

	void CorrelationTermStructureStrike::accept(AcyclicVisitor& v) {
		Visitor<CorrelationTermStructureStrike>* v1 =
			dynamic_cast<Visitor<CorrelationTermStructureStrike>*>(&v);
		if (v1 != 0)
			v1->visit(*this);
		else
			QL_FAIL("not a local-correlation term structure visitor");
	}
}
