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

/*! \file fxfwdratehelper.hpp
    \brief FX-Points rate helpers
*/

#ifndef quantlib_fxfwdratehelper_hpp
#define quantlib_fxfwdratehelper_hpp

#include <ql/termstructures/bootstraphelper.hpp>
#include <ql/time/calendar.hpp>
#include <ql/time/daycounter.hpp>
#include <ql/cashflows/iborcoupon.hpp>

namespace QuantLib {

    class Quote;

    typedef BootstrapHelper<YieldTermStructure> RateHelper;
    typedef RelativeDateBootstrapHelper<YieldTermStructure>
                                                        RelativeDateRateHelper;

 

    //! Rate helper for bootstrapping fx points relying on interest rate parity
    //!
    //! For a currency pair CCY1/CCY2 the first currency is the base currency, while the second is the counter currency. 
    //! The fx rate quotes the second (counter) currency in terms of one unit of the first (base) currency. 
    //!
    //! Neglecting a spot lag from formal interest rate parity one obtains the number of fx points as: 
    //!
    //! FX-Points/units := FX_t - FX_0 = FX_0*[DiscountFactorBase(T)/DiscountFactorCounter(T) - 1]
    //!
    //! FX-Points are quoted in terms of a currency pair specific unit, e.g. for EUR/USD this unit is 10000.
    //!
    //! The yield curves entered into the rate helper usually are the appropriate OIS-curves of the two currencies. 
    //!
    //! The enum bootstrapBaseOrCounter controls, which cross currency curve isgoing to be bootstrapped. The choice 'Base' 
    //! determines a cross currency curve for the base currency while applying the OIS curve to the counter currency and for
    //! the choice 'Counter' it is the other way around. 


    class FxFwdRateHelper : public RelativeDateRateHelper {
      public:
        enum BootstrapType { Base = -1, Counter = 1};
        FxFwdRateHelper( Currency baseCurrency,
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
                         BootstrapType bootstrapBaseOrCounter
                         );

        //! \name RateHelper interface
        //@{
        Real impliedQuote() const;
        void setTermStructure(YieldTermStructure*);
        //@}
        //! \name Visitability
        //@{
        void accept(AcyclicVisitor&);
        //@}
    protected:
        //currency pair
        Currency baseCurrency_;
        Currency counterCurrency_;
        Rate   fxSpot_;
        //general parameters
        Natural spotLag_;
        Calendar spotLagCalendar_;
        BusinessDayConvention spotLagConvention_;
        Period swapTerm_;
        Real   unit_;

        Handle<YieldTermStructure> baseCcyDiscTermStructureHandle_;
        Handle<YieldTermStructure> cntrCcyDiscTermStructureHandle_;

        //Bootstrap direction, i.e. which curve is going to be bootstrapped
        BootstrapType bootstrapBaseOrCounter_;

        //curve going to be optimized
        RelinkableHandle<YieldTermStructure> baseCcyDiscRelinkableHandle_;
        RelinkableHandle<YieldTermStructure> cntrCcyDiscRelinkableHandle_;

        void initializeDates();
     };

}

#endif
