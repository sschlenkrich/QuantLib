/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2010 Sebastian Schlenkrich
*/

/*! \file bondoptionengine.cpp
    \brief engine for bond options with Hull White model
*/


#include <time.h>
#include <boost/lexical_cast.hpp>

#include <ql/settings.hpp>
#include <ql/experimental/templatemodels/hullwhite/fixedratebondoption.hpp>
#include <ql/experimental/templatemodels/hullwhite/bondoptionengine.hpp>
#include <ql/cashflows/coupon.hpp>
#include <ql/cashflows/fixedratecoupon.hpp>
#include <ql/cashflows/simplecashflow.hpp>

namespace QuantLib {

    // calibrate model based on given swaptions
    void BondOptionEngine::calibrateModel( std::vector< ext::shared_ptr<Swaption> >          swaptions,
                                           const bool                                        contTenorSpread,
                                           const Real                                        tolVola) {
        std::sort(swaptions.begin(), swaptions.end(), lessByExerciseFirstDate);
        // set up inputs for model calibration
        std::vector<Real>                  exercTimes, strikeVals, b76Prices;
        std::vector< std::vector<Real> >   startTimes, payTimes, cashFlows;
        std::vector<Option::Type>          callOrPut;
        // get discount curve and conventions from model
        Handle<YieldTermStructure> discCurve = model_->termStructure();
        DayCounter                 dc        = model_->termStructure()->dayCounter();
        Date                       today     = model_->termStructure()->referenceDate();
        // iterate swaptions...
        referenceSwaptions_.clear();
        for (Size k=0; k<swaptions.size(); ++k) {
            // skip swaptions with equal exercise date since we calibrate by bootstrapping
            if ( (k>0) && (swaptions[k]->exercise()->date(0)==swaptions[k-1]->exercise()->date(0)) ) continue;
            // save reference swaption
            referenceSwaptions_.push_back(swaptions[k]);
            // build an equivalent bond option
            FixedRateBondOption bondOption(swaptions[k], discCurve, contTenorSpread);
            // extract option details
            exercTimes.push_back(dc.yearFraction(today,bondOption.exerciseDates()[0]));
            strikeVals.push_back(bondOption.dirtyStrikeValues()[0]);
            b76Prices.push_back(swaptions[k]->NPV());  // assume we have an engine and market data attached
            callOrPut.push_back(bondOption.callOrPut());
            // extract underlying details
            std::vector<Real> bondStartTimes, bondPayTimes;
            for (Size i=0; i<bondOption.cashflowValues().size(); ++i) {
                bondStartTimes.push_back( dc.yearFraction(today,bondOption.startDates()[i]) );
                bondPayTimes.push_back( dc.yearFraction(today,bondOption.payDates()[i]) );
            }
            startTimes.push_back(bondStartTimes);
            payTimes.push_back(bondPayTimes);
            cashFlows.push_back(bondOption.cashflowValues());
        }
        // calibrate Hull White model
        model_->BermudanCalibration(exercTimes, strikeVals, b76Prices, startTimes, payTimes, cashFlows, callOrPut, tolVola);
    }

    void BondOptionEngine::calculate() const {
        std::vector<Time> startTimes;
        std::vector<Time> payTimes;
        std::vector<Real> cashFlowValues;
        std::vector<Time> exerciseTimes;
        std::vector<Real> strikeValues;
        DayCounter dayCounter = model_->termStructure()->dayCounter();
        Date today = model_->termStructure()->referenceDate();
        Date startDate;
        ext::shared_ptr<Coupon> coupon;
        // set up cash flows
        for (Size i=0; i<arguments_.cashflows.size(); ++i) {
            coupon = ext::dynamic_pointer_cast<Coupon>(arguments_.cashflows[i]);
            if (coupon) { // cast is ok
                startDate = coupon->accrualStartDate();
            }
            else { // cash flow is no coupon, assume redemption payment, startDate = payDate
                startDate = arguments_.cashflows[i]->date();
            }
            if (startDate>today) { // consider only coupons with startDate later than today
                startTimes.push_back(dayCounter.yearFraction(today,startDate));
                payTimes.push_back(dayCounter.yearFraction(today,arguments_.cashflows[i]->date()));
                cashFlowValues.push_back(arguments_.cashflows[i]->amount());
            }
        }
        // set up exercises
        Size Nexc = std::min(arguments_.exerciseDates.size(), arguments_.dirtyStrikeValues.size());
        for (Size i=0; i<Nexc; ++i) {
            if (arguments_.exerciseDates[i]>today) { // consider only exercises later than today
                exerciseTimes.push_back(dayCounter.yearFraction(today,arguments_.exerciseDates[i]));
                strikeValues.push_back(arguments_.dirtyStrikeValues[i]);
            }
        }
        // do not calibrate in pricing engine unless there is a calibrator

        // evaluate Bermudan bond option
        clock_t start_clock = clock();
        results_.value = model_->BermudanBondOption(exerciseTimes,strikeValues,
            startTimes,payTimes,cashFlowValues,arguments_.callOrPut,dimension_,gridRadius_,bermudanTolerance_);
        clock_t end_clock = clock();
        // report additional results here ...
        results_.additionalResults["runtime"] = (Real)(end_clock - start_clock)/CLOCKS_PER_SEC;
        std::vector<Real> europeansAnalytical = model_->europeansAnalytical();
        std::vector<Real> europeansNumerical  = model_->europeansNumerical();
        // (absolute) error is estimated by corresponding European prices
        Real errorEstimate = 0.0;
        for (Size i=0; i<std::min(europeansAnalytical.size(),europeansNumerical.size()); ++i) {
            errorEstimate = std::max(errorEstimate,fabs(europeansNumerical[i]-europeansAnalytical[i]));
        }
        results_.errorEstimate = errorEstimate;
        // report adittional results for
        for (Size i=0; i<referenceSwaptions_.size(); ++i) {
            std::string name = "black76price_";
            name += boost::lexical_cast<std::string>( i+1 );
            results_.additionalResults[name] = referenceSwaptions_[i]->NPV();
        }
        for (Size i=0; i<referenceSwaptions_.size(); ++i) {
            std::string name = "black76vola_";
            name += boost::lexical_cast<std::string>( i+1 );
            std::map<std::string,boost::any>::const_iterator stDev = referenceSwaptions_[i]->additionalResults().find("stdDev");
            if (stDev!=referenceSwaptions_[i]->additionalResults().end()) {
                results_.additionalResults[name] = 
                    boost::any_cast<Real>(stDev->second) / std::sqrt(model_->termStructure()->dayCounter().yearFraction(today,referenceSwaptions_[i]->exercise()->dates()[0]));
            }
        }
        for (Size i=0; i<referenceSwaptions_.size(); ++i) {
            std::string name = "black76vega_";
            name += boost::lexical_cast<std::string>( i+1 );
            std::map<std::string,boost::any>::const_iterator vega = referenceSwaptions_[i]->additionalResults().find("vega");
            if (vega!=referenceSwaptions_[i]->additionalResults().end()) {
                results_.additionalResults[name] = 	boost::any_cast<Real>(vega->second);
            }
        }
        // if we have an AD-enabeled model report vega(s) here...
        ext::shared_ptr<MinimADHullWhiteModel> amodel = ext::dynamic_pointer_cast<MinimADHullWhiteModel>(model_);
        if (amodel) {
            // derivative of Bermudan price w.r.t. short rate vola
            std::vector<QuantLib::Real> vegas = amodel->bermudanVega();
            // differentiate calibration; short rate vola w.r.t. B76 prices
            for (Size j=vegas.size(); j>0; --j) {
                for (Size i=j; i<vegas.size(); ++i) vegas[j-1] -= vegas[i]*amodel->calibrationJacobian()[i][j-1];
                vegas[j-1] /= amodel->calibrationJacobian()[j-1][j-1];
            }
            // finally differentiate reference prices w.r.t. Black'76 volas
            for (Size i=0; i<std::min(vegas.size(),referenceSwaptions_.size()); ++i) {
                std::map<std::string,boost::any>::const_iterator vega = referenceSwaptions_[i]->additionalResults().find("vega");
                if (vega!=referenceSwaptions_[i]->additionalResults().end()) {
                    vegas[i] *= boost::any_cast<Real>(vega->second);
                }
            }
            // the sum of vegas represents the sensitivity w.r.t. to a parallel shift of the B76 vola surface
            QuantLib::Real vega=0.0;
            for (Size i=0; i<vegas.size(); ++i) vega += vegas[i];
            // store Bermudan and reference European vega
            results_.additionalResults["vega"] = vega;
            results_.additionalResults["vegas_size"] = (1.0*vegas.size());
            for (Size i=0; i<vegas.size(); ++i) {
                std::string name = "vegas_";
                name += boost::lexical_cast<std::string>( i+1 );
                results_.additionalResults[name] = vegas[i]; 
            }
        }
    }

}

