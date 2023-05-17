/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*

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

/*! \file indextermstructure.hpp
\brief general purpose term structure of interpolated values
*/

#ifndef quantlib_index_term_structure_hpp
#define quantlib_index_term_structure_hpp

#include <ql/termstructure.hpp>
#include <ql/quote.hpp>
#include <ql/math/interpolation.hpp>
#include <ql/math/interpolations/cubicinterpolation.hpp>


namespace QuantLib {

    class IndexTermStructure : public TermStructure {
    protected:
        std::vector<Date>                dates_;
        std::vector<Real>                values_;
        ext::shared_ptr<Interpolation>   timeInterpol_;
        inline void setupInterpolation() {
            std::vector<Time> times(dates_.size());
            for (Size i=0; i<dates_.size(); ++i) times[i] = dayCounter().yearFraction(referenceDate(),dates_[i]);
            // for the moment we leave it with hard-coded cubic interpolation
            timeInterpol_ = ext::shared_ptr<Interpolation>( new Interpolation(Cubic().interpolate(times.begin(),times.end(),values_.begin())) );
        }
    public:
        /*! \name Constructors
        See the TermStructure documentation for issues regarding
        constructors.
        */
        //@{
        IndexTermStructure(
            const DayCounter&                  dc     = DayCounter(),
            const std::vector<Date>&           dates  = std::vector<Date>(),
            const std::vector<Real>&           values = std::vector<Real>(),
            const ext::shared_ptr<Interpolation>&   timeInterpol = ext::shared_ptr<Interpolation>()
            ) : TermStructure(dc), dates_(dates), values_(values), timeInterpol_(timeInterpol) {
                setupInterpolation();
        }
        IndexTermStructure(
            const Date&                        referenceDate,
            const Calendar&                    cal    = Calendar(),
            const DayCounter&                  dc     = DayCounter(),
            const std::vector<Date>&           dates  = std::vector<Date>(),
            const std::vector<Real>&           values = std::vector<Real>()
            ) : TermStructure(referenceDate,cal,dc), dates_(dates), values_(values) {
                setupInterpolation();
        }
        IndexTermStructure(
            Natural                            settlementDays,
            const Calendar&                    cal,
            const DayCounter&                  dc     = DayCounter(),
            const std::vector<Date>&           dates  = std::vector<Date>(),
            const std::vector<Real>&           values = std::vector<Real>()
            ) : TermStructure(settlementDays,cal,dc), dates_(dates), values_(values) {
                setupInterpolation();
        }
        // let the user specify the time interpolation; only act as a wrapper
        IndexTermStructure(
            const Date&                        referenceDate,
            const Calendar&                    cal    = Calendar(),
            const DayCounter&                  dc     = DayCounter(),
            const ext::shared_ptr<Interpolation>&   timeInterpol = ext::shared_ptr<Interpolation>()
            ) : TermStructure(referenceDate,cal,dc), timeInterpol_(timeInterpol) { }

        inline Real value ( Time t, bool extrapolate = true) const {
            return timeInterpol_->operator()(t,extrapolate);
        }

        inline Real value (const Date& d, bool extrapolate = true) const {
            Time t = dayCounter().yearFraction(referenceDate(),d);
            return value(t,extrapolate);
        }

        // TermStructure interface

        //! the latest date for which the curve can return values
        virtual Date maxDate() const { return dates_.back(); }

    };


}

#endif //quantlib_index_term_structure_hpp
