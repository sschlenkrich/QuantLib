/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2010 Sebastian Schlenkrich

*/

/*! \file hullwhitemodels.cpp
    \brief evaluate vegas in the Hull White model using MinimAD
*/


#include <ql/mathconstants.hpp>
#include <ql/experimental/templatemodels/hullwhite/hullwhitemodels.hpp>



namespace QuantLib {
    void MinimADHullWhiteModel::cloneModel() {
        // initialise independent variables
        avolaValues_.resize(this->volaValues().size());
        for (Size k=0; k<this->volaValues().size(); ++k) avolaValues_[k] = this->volaValues()[k];
        // clone the passive model
        amodel_ = ext::shared_ptr< HullWhiteModelT<QuantLib::Time,QuantLib::Real,ActiveType> >(
            new HullWhiteModelT<QuantLib::Time,QuantLib::Real,ActiveType>(this->termStructure(), this->mean(), this->volaTimes(), avolaValues_));
        amodel_->setEstimateAccuracy(this->estimateAccuracy());
    }

    QuantLib::Real MinimADHullWhiteModel::BermudanBondOption( 
                                const std::vector<QuantLib::Time>&  exercDates,    // option's exercise dates (equal settlment)
                                const std::vector<QuantLib::Real>&  strikeValues,  // strike payed at excercise dates
                                // option's underlying
                                const std::vector<QuantLib::Time>&  startDates,   // start dates of coupon period
                                const std::vector<QuantLib::Time>&  payDates,     // pay dates of coupon perid
                                const std::vector<QuantLib::Real>&  cashFlows,    // fixed coupon payments (absolut value)
                                const Option::Type                  cop,          // call (1) or put (-1) option
                                // discretisation properties
                                const size_t                        dim,
                                const QuantLib::Real                gridRadius,
                                const QuantLib::Real                tol ) {
        this->cloneModel();
        // evaluate the function
        ActiveType aresult = amodel_->BermudanBondOption( exercDates, strikeValues, startDates, payDates, cashFlows, cop, dim, gridRadius, tol);
        // fetch reference European prices
        this->europeansAnalytical_.resize(amodel_->europeansAnalytical().size());
        for (Size k=0; k<this->europeansAnalytical_.size(); ++k) this->europeansAnalytical_[k] = amodel_->europeansAnalytical()[k].val();
        this->europeansNumerical_.resize(amodel_->europeansNumerical().size());
        for (Size k=0; k<this->europeansNumerical_.size(); ++k) this->europeansNumerical_[k] = amodel_->europeansNumerical()[k].val();
        // evaluate derivatives
        europeansAnalyticalVega_.resize(amodel_->europeansAnalytical().size());
        for (Size k=0; k<europeansAnalyticalVega_.size(); ++k) {
            europeansAnalyticalVega_[k].resize(amodel_->volaValues().size());
            for (Size i=0; i<europeansAnalyticalVega_[k].size(); ++i) {
                europeansAnalyticalVega_[k][i] = 
                    amodel_->europeansAnalytical()[k] % amodel_->volaValues()[i];
            }
        }
        europeansNumericalVega_.resize(amodel_->europeansNumerical().size());
        for (Size k=0; k<europeansNumericalVega_.size(); ++k) {
            europeansNumericalVega_[k].resize(amodel_->volaValues().size());
            for (Size i=0; i<europeansNumericalVega_[k].size(); ++i) {
                europeansNumericalVega_[k][i] = 
                    amodel_->europeansNumerical()[k] % amodel_->volaValues()[i];
            }
        }
        bermudanVega_.resize(amodel_->volaValues().size());
        for (Size i=0; i<bermudanVega_.size(); ++i) bermudanVega_[i] = aresult % amodel_->volaValues()[i];

        // finally return passive result
        return aresult.val();
    }

    const std::vector<QuantLib::Real>& MinimADHullWhiteModel::BermudanCalibration (
                                       const std::vector<QuantLib::Time>&  exercDates,   // option's exercise dates (equal settlment)
                                       const std::vector<QuantLib::Real>&  strikeValues, // strike payed at excercise dates
                                       const std::vector<QuantLib::Real>&  b76Prices,    // reference European prices
                                       // option's underlying
                                       const std::vector< std::vector<QuantLib::Time> >&  startDates,   // start dates of coupon period
                                       const std::vector< std::vector<QuantLib::Time> >&  payDates,     // pay dates of coupon perid
                                       const std::vector< std::vector<QuantLib::Real> >&  cashFlows,    // fixed coupon payments (absolut value)
                                       const std::vector< Option::Type >&                 cop,          // call (1) or put (-1) option
                                       // calibration parameters
                                       const QuantLib::Real                tol_vola ) {  // absolut tolerance in short rate volatility
        // solve passive calibration problem to get fixed point
        RealHullWhiteModel::BermudanCalibration(exercDates,strikeValues,b76Prices,startDates,payDates,cashFlows,cop,tol_vola);
        // initialise active HW model
        cloneModel();
        // evaluate active coupon bond options equivalent to calibration
        Size Nexc;
        Nexc = std::min(exercDates.size(),strikeValues.size());
        Nexc = std::min(Nexc,b76Prices.size());
        std::vector<ActiveType> bondOptions(Nexc);
        for (Size k=0; k<Nexc; ++k) {
            bondOptions[k] = amodel_->CouponBondOption( exercDates[k], strikeValues[k], startDates[std::min(k,startDates.size()-1)], payDates[std::min(k,payDates.size()-1)], cashFlows[std::min(k,cashFlows.size()-1)], cop[std::min(k,cop.size()-1)]);
        }
        // evaluate derivatives
        calibrationJacobian_.resize(bondOptions.size());
        for (Size k=0; k<Nexc; ++k) {
            calibrationJacobian_[k].resize(amodel_->volaValues().size());
            for (Size i=0; i<calibrationJacobian_[k].size(); ++i) {
                calibrationJacobian_[k][i] =
                    bondOptions[k] % amodel_->volaValues()[i];
            }
        }
        return volaValues_;
    }


}


