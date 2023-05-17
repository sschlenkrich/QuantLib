/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2010 Sebastian Schlenkrich
*/

/*! \file fixedratebondoption.hpp
    \brief (Bermudan) fixed-rate bond option
*/

#ifndef quantlib_fixedratebondoption_hpp
#define quantlib_fixedratebondoption_hpp

#include <ql/instruments/bonds/fixedratebond.hpp>
#include <ql/instruments/swaption.hpp>
#include <ql/termstructures/yieldtermstructure.hpp>
#include <ql/time/date.hpp>
#include <ql/option.hpp>

namespace QuantLib {

    class FixedRateBondOption : public Instrument {
    protected:
        Leg                                  cashflows_;
        std::vector<Date>                    exerciseDates_;
        std::vector<Real>                    dirtyStrikeValues_;
        Option::Type                         callOrPut_;
        void setupArguments(PricingEngine::arguments*) const;
    public:
        FixedRateBondOption ( const ext::shared_ptr<FixedRateBond>& underlyingBond,
                              const std::vector<Date>& exerciseDates,
                              const std::vector<Real>& dirtyStrikeValues,
                              const Option::Type callOrPut )
                              : cashflows_(underlyingBond->cashflows()), exerciseDates_(exerciseDates),
                                dirtyStrikeValues_(dirtyStrikeValues),  callOrPut_(callOrPut) {}

        // constructor to map a swaption to bond option according to spread model
        FixedRateBondOption ( const ext::shared_ptr<Swaption>& swaption,
                              const Handle<YieldTermStructure>& discountCurve,
                              bool                              contTenorSpread = true );

        bool isExpired() const { // Instrument interface
            return exerciseDates_.back()<= Settings::instance().evaluationDate();
        }

        // inspectors
        const std::vector<Date>&  exerciseDates()      const { return exerciseDates_; }
        const std::vector<Real>&  dirtyStrikeValues()  const { return dirtyStrikeValues_; }
        const Option::Type        callOrPut()          const { return callOrPut_; }
        const std::vector< QuantLib::Date > startDates();
        const std::vector< QuantLib::Date > payDates();
        const std::vector< QuantLib::Real > cashflowValues();

        class arguments : public PricingEngine::arguments {
        public:
            Leg cashflows;
            std::vector<Date> exerciseDates;
            std::vector<Real> dirtyStrikeValues;
            Option::Type callOrPut;
            void validate() const {} // add some meaningfull checks here
        };

        class results : public Instrument::results {
        public:
            // value (NPV), errorEstimate, additionalResults and valuationDate 
            // are declared in Instrument::results
            // here we may add the sensitivities
            // fetchResults is inherited from Instrument and currently NOT overloaded
        };

        class engine : public GenericEngine<arguments,results> {};
        // derived engines have to evaluate
        // Real value, i.e., NPV
        // Real errorEstimate, i.e., numerical vs. analytical results
        // Date valuationDate, date until discounted
        // std::map<std::string,boost::any> additionalResults, european reference prices

    };

}

#endif /* #ifndef quantlib_fixedratebondoption_hpp */
