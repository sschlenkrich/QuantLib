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

/*! \file rangeaccrual.hpp
    \brief range-accrual coupon
*/

#ifndef quantlib_cms_range_accrual_fixed_h
#define quantlib_cms_range_accrual_fixed_h

#include <ql/cashflows/couponpricer.hpp>
#include <ql/cashflows/fixedratecoupon.hpp>
#include <ql/time/schedule.hpp>
#include <vector>

namespace QuantLib {

    class SwapIndex;
    class SwaptionVolatilityStructure;
    class CmsRangeAccrualFixedCouponPricer;

    class CmsRangeAccrualFixedCoupon: public FixedRateCoupon {

      public:
        CmsRangeAccrualFixedCoupon(
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
            const Date& refPeriodStart = Date(),
            const Date& refPeriodEnd = Date(),
            const Date& exCouponDate = Date());

        CmsRangeAccrualFixedCoupon(
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
            const Date& refPeriodStart = Date(),
            const Date& refPeriodEnd = Date(),
            const Date& exCouponDate = Date());

        //@}
        //! \name LazyObject interface
        //@{
        void performCalculations() const override;
        //@}
        //! \name CashFlow interface
        //@{
        Real amount() const override;
        //@}


        ext::shared_ptr<Schedule> observationsSchedule() const { return observationsSchedule_; }
        ext::shared_ptr<SwapIndex> cmsIndex() const { return cmsIndex_; }
        Real lowerTrigger() const { return lowerTrigger_; }
        Real upperTrigger() const { return upperTrigger_; }
        Real rangeAccrual() const;

        //! \name Visitability
        //@{
        void accept(AcyclicVisitor&) override;
        //@}

        void setPricer(const ext::shared_ptr<CmsRangeAccrualFixedCouponPricer>&);

        std::map<std::string, Real>& additionalResults() const { return additionalResults_; }

      private:

        const ext::shared_ptr<Schedule> observationsSchedule_;
        ext::shared_ptr<SwapIndex> cmsIndex_;
        std::vector<Date> observationDates_;
        Real lowerTrigger_;
        Real upperTrigger_;

        ext::shared_ptr<CmsRangeAccrualFixedCouponPricer> pricer_;
        mutable Real rangeAccrual_;
        mutable std::map<std::string, Real> additionalResults_;
     };


    class CmsRangeAccrualFixedCouponPricer
    : public virtual Observer,
      public virtual Observable {
      public:

        CmsRangeAccrualFixedCouponPricer(
            Handle<SwaptionVolatilityStructure> swaptionVolatility
        );

        void initialize(const CmsRangeAccrualFixedCoupon& coupon);

        Real rangeAccrual() const;

        std::map<std::string, Real>& additionalResults() const { return additionalResults_; }

      //! \name Observer interface
      //@{
      void update() override { notifyObservers(); }
      //@}

      protected:
        Handle<SwaptionVolatilityStructure> swaptionVolatility_;
        Real rangeAccrual_;
        mutable std::map<std::string, Real> additionalResults_;
    };

}


#endif
