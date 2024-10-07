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

/*! \file fxindex.hpp
    \brief base class for fx indexes based on equity index
*/

#ifndef quantlib_fxindex_hpp
#define quantlib_fxindex_hpp

#include <ql/indexes/equityindex.hpp>


namespace QuantLib {

    //! Base class for FX indexes
    /*! The FX index object allows to retrieve past fixings,
        as well as project future fixings.
        
        Forward is calculated as:
        \f[
        I(t, T) = I(t, t) \frac{P_{F}(t, T)}{P_{D}(t, T)},
        \f]
        where \f$ I(t, t) \f$ is today's value of the index,
        \f$ P_{F}(t, T) \f$ is a discount factor of the foreign
		currency curve at future time \f$ T \f$, and
		\f$ P_{D}(t, T) \f$ is a discount factor of the domestic
		curve at future time \f$ T \f$.

        To forecast future fixings, the user can either provide a
        handle to the current index spot. If spot handle is empty,
        today's fixing will be used, instead.
    */
    class FxIndex : public EquityIndex {
      public:
        FxIndex(std::string name,
                Calendar fixingCalendar,
                Handle<YieldTermStructure> domInterest,
                Handle<YieldTermStructure> forInterest,
                Handle<Quote> spot = {});

        //! the domestic rate curve used to forecast fixings
        Handle<YieldTermStructure> domesticInterestRateCurve() const {
            return equityInterestRateCurve();
        }
        //! the foreign rate curve used to forecast fixings
        Handle<YieldTermStructure> foreignInterestRateCurve() const {
            return equityDividendCurve();
        }

    };

}

#endif