/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2010 Sebastian Schlenkrch
*/

/*! \file fixedratebondoption.cpp
    \brief Bermudan) fixed-rate bond option
*/


#include <ql/settings.hpp>
#include <ql/exercise.hpp>
#include <ql/cashflows/coupon.hpp>
#include <ql/cashflows/simplecashflow.hpp>
#include <ql/cashflows/fixedratecoupon.hpp>
#include <ql/experimental/templatemodels/hullwhite/fixedratebondoption.hpp>


namespace QuantLib {

        // constructor to map a swaption to bond option according to spread model
    FixedRateBondOption::FixedRateBondOption (
                        const ext::shared_ptr<Swaption>& swaption,
                        const Handle<YieldTermStructure>& discountCurve,
                        bool                              contTenorSpread ) {
        // evaluate callOrPut_, exerciseDates_, dirtyStrikeValues_, and cashFlows_
        callOrPut_ = (swaption->underlyingSwap()->type()==VanillaSwap::Type::Receiver) ? Option::Type::Call : Option::Type::Put;
        // consider only future exercise dates
        for (Size k=0; k<swaption->exercise()->dates().size(); ++k) {
            if (swaption->exercise()->dates()[k]>Settings::instance().evaluationDate()) {
                exerciseDates_.push_back(swaption->exercise()->dates()[k]);
            }
        }
        Leg floatLeg = swaption->underlyingSwap()->floatingLeg();
        // evaluate strike paid at exercise, assume deterministic strike paid at next start date (settlement date)
        for (Size k=0; k<exerciseDates_.size(); ++k) {
            Size floatIdx=0;
            while ( (exerciseDates_[k]>(ext::dynamic_pointer_cast<Coupon>(floatLeg[floatIdx]))->accrualStartDate()) && (floatIdx<floatLeg.size()-1)) ++floatIdx;
            if (exerciseDates_[k]>(ext::dynamic_pointer_cast<Coupon>(floatLeg[floatIdx])->accrualStartDate())) {
                dirtyStrikeValues_.push_back(0.0);  // if there is no coupon left the strike is trivially equal to zero
            } else {
                Real nominal      = (ext::dynamic_pointer_cast<Coupon>(floatLeg[floatIdx]))->nominal();
                Real dfExercise   = discountCurve->discount(exerciseDates_[k]);
                Real dfSettlement = discountCurve->discount((ext::dynamic_pointer_cast<Coupon>(floatLeg[floatIdx]))->accrualStartDate());
                dirtyStrikeValues_.push_back(nominal*dfSettlement/dfExercise);				
            }
        }
        // evaluate floating leg deterministic spreads
        Leg spreadLeg;
        for (Size k=0; k<floatLeg.size(); ++k) {
            ext::shared_ptr<Coupon> coupon = ext::dynamic_pointer_cast<Coupon>(floatLeg[k]);
            if (!coupon) QL_FAIL("FloatingLeg CashFlow is no Coupon.");
            Date startDate = coupon->accrualStartDate();
            if (startDate>Settings::instance().evaluationDate()) { // consider only future cash flows
                Date endDate = coupon->accrualEndDate();
                Rate liborForwardRate = coupon->rate();
                Rate discForwardRate = (discountCurve->discount(startDate)/discountCurve->discount(endDate)-1.0)/coupon->accrualPeriod();
                Rate spread;
                Date payDate;
                if (contTenorSpread) {
                    // Db = (1 + Delta L^libor) / (1 + Delta L^ois)
                    // spread (Db - 1) paid at startDate
                    spread = ((1.0 + coupon->accrualPeriod()*liborForwardRate) / (1.0 + coupon->accrualPeriod()*discForwardRate) - 1.0) / coupon->accrualPeriod();
                    payDate = startDate;
                } else {
                    // spread L^libor - L^ois
                    spread = liborForwardRate - discForwardRate;
                    payDate = coupon->date();
                }
                spreadLeg.push_back(ext::shared_ptr<CashFlow>(new FixedRateCoupon(payDate,-1.0*coupon->nominal(),spread,coupon->dayCounter(),startDate,endDate)));
            }  // if ...
        }  // for ...
        // merge fixed leg and spreads according to start date
        Leg fixedLeg = swaption->underlyingSwap()->fixedLeg();
        Size i=0, j=0;
        while (i<fixedLeg.size() || j<spreadLeg.size()) {
            if (i>=fixedLeg.size())  { cashflows_.push_back(spreadLeg[j]); ++j; continue; }
            if (j>=spreadLeg.size()) { cashflows_.push_back(fixedLeg[i]);  ++i; continue; }
            // here we have i<fixedLeg.size() && j<spreadLeg.size()
            ext::shared_ptr<Coupon> fixedCoupon = ext::dynamic_pointer_cast<Coupon>(fixedLeg[i]);
            if (!fixedCoupon) QL_FAIL("FixedLeg CashFlow is no Coupon.");
            ext::shared_ptr<Coupon> spreadCoupon = ext::dynamic_pointer_cast<Coupon>(spreadLeg[j]);
            if (!spreadCoupon) QL_FAIL("SpreadLeg CashFlow is no Coupon.");
            if (fixedCoupon->accrualStartDate()<=spreadCoupon->accrualStartDate()) { cashflows_.push_back(fixedLeg[i]);  ++i; }
            else                                                                   { cashflows_.push_back(spreadLeg[j]); ++j; }
        }  // while ...
        // finally, add the notional at the last date
        ext::shared_ptr<Coupon> lastFloatCoupon = ext::dynamic_pointer_cast<Coupon>(floatLeg.back());
        cashflows_.push_back(ext::shared_ptr<CashFlow>(new SimpleCashFlow(lastFloatCoupon->nominal(),lastFloatCoupon->accrualEndDate())));
    }

     const std::vector< QuantLib::Date > FixedRateBondOption::startDates() {
         std::vector< QuantLib::Date > dates;
         for (Size k=0; k<cashflows_.size(); ++k) {
             ext::shared_ptr<Coupon> coupon = ext::dynamic_pointer_cast<Coupon>(cashflows_[k]);
             Date startDate = cashflows_[k]->date(); // default, e.g. for redemptions
             if (coupon) startDate = coupon->accrualStartDate();
             dates.push_back(startDate);
         }
         return dates;
     }

     const std::vector< QuantLib::Date > FixedRateBondOption::payDates() {
         std::vector< QuantLib::Date > dates;
         for (Size k=0; k<cashflows_.size(); ++k) dates.push_back(cashflows_[k]->date());
         return dates;
     }

     const std::vector< QuantLib::Real > FixedRateBondOption::cashflowValues() {
         std::vector< QuantLib::Real > values;
         for (Size k=0; k<cashflows_.size(); ++k) values.push_back(cashflows_[k]->amount());
         return values;
     }

     void FixedRateBondOption::setupArguments(PricingEngine::arguments* args) const {
        FixedRateBondOption::arguments* arguments = dynamic_cast<FixedRateBondOption::arguments*>(args);
        QL_REQUIRE(arguments != 0, "wrong argument type");

        //arguments->settlementDate = underlyingBond_->settlementDate();
        //arguments->cashflows = underlyingBond_->cashflows();
        //arguments->calendar = underlyingBond_->calendar();
        arguments->cashflows = cashflows_;
        arguments->exerciseDates = exerciseDates_;
        arguments->dirtyStrikeValues = dirtyStrikeValues_;
        arguments->callOrPut = callOrPut_;
    }

}

