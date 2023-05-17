/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2010 Sebastian Schlenkrich
*/

/*! \file spreadyieldtermstructure.hpp
    \brief compounded yield term structure evaluating exp{ - ( z1 + alpha z2 ) t }
*/

#ifndef quantlib_spreadyieldtermstructure_hpp
#define quantlib_spreadyieldtermstructure_hpp

#include <ql/handle.hpp>
#include <ql/termstructures/yieldtermstructure.hpp>

namespace QuantLib {

    class SpreadYTS : public YieldTermStructure {
    private:
        Handle<YieldTermStructure> baseCurve_;
        Handle<YieldTermStructure> sprdCurve_;
        Real alpha_;
    public :
        SpreadYTS( const Handle<YieldTermStructure>& baseCurve = Handle<YieldTermStructure>(),
                   const Handle<YieldTermStructure>& sprdCurve = Handle<YieldTermStructure>(),
                   const Real alpha = 1.0 )
                   : baseCurve_(baseCurve), sprdCurve_(sprdCurve), alpha_(alpha),
                       YieldTermStructure(baseCurve->referenceDate(),baseCurve->calendar(),baseCurve->dayCounter()) 
        {
            registerWith(baseCurve_);
            registerWith(sprdCurve_);
        }
        Date maxDate() const { return baseCurve_->maxDate(); };
    protected :
        inline DiscountFactor discountImpl(Time t) const {
            if (alpha_ == 1.0)  return baseCurve_->discount(t) * sprdCurve_->discount(t);
            if (alpha_ ==-1.0)  return baseCurve_->discount(t) / sprdCurve_->discount(t);
            if (alpha_ == 0.0)  return baseCurve_->discount(t);
            return                     baseCurve_->discount(t) * pow(sprdCurve_->discount(t),alpha_);
        }

    };

}

#endif
