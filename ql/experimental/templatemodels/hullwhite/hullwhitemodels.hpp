/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2010 Sebastian Schlenkrich

*/

/*! \file hullwhitemodels.hpp
    \brief evaluate vegas in the Hull White model using MinimAD2
*/


#ifndef quantlib_templatehullwhitemodels_hpp
#define quantlib_templatehullwhitemodels_hpp

//#include <ql/experimental/templatehullwhite/adtageo/adtageo.hpp>
#include <ql/experimental/templatemodels/auxilliaries/minimADVariable2T.hpp>
#include <ql/experimental/templatemodels/hullwhite/hullwhitemodelT.hpp>
#include <ql/experimental/templatemodels/hullwhite/g2ppT.hpp>


namespace QuantLib {

    //typedef ADTAGEO::daglad ActiveType;
    typedef MinimAD::Variable<QuantLib::Real> ActiveType;

    // basic binding of template parameters
    typedef HullWhiteModelT<QuantLib::Time,QuantLib::Real,QuantLib::Real> RealHullWhiteModel;

    class MinimADHullWhiteModel : public RealHullWhiteModel {
    private:
        // independant variables
        std::vector<ActiveType> avolaValues_;
        // a clone of the model with active data types
        ext::shared_ptr< HullWhiteModelT<QuantLib::Time,QuantLib::Real,ActiveType> > amodel_;
        // derivatives evaluated
        std::vector< std::vector<QuantLib::Real> > calibrationJacobian_;
        std::vector< std::vector<QuantLib::Real> > europeansAnalyticalVega_;
        std::vector< std::vector<QuantLib::Real> > europeansNumericalVega_;
        std::vector<QuantLib::Real>                bermudanVega_;
        void cloneModel();
    public:
        MinimADHullWhiteModel(  const Handle<YieldTermStructure>&   termStructure,
                                const QuantLib::Real                mean,
                                const std::vector<QuantLib::Time>&  volaDates,
                                const std::vector<QuantLib::Real>&  volaValues )
            : RealHullWhiteModel(termStructure, mean, volaDates, volaValues) { }
        // Bermudan bond option plus vegas
        QuantLib::Real BermudanBondOption( 
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
                                const QuantLib::Real                tol );
        const std::vector<QuantLib::Real>& BermudanCalibration (
                                       const std::vector<QuantLib::Time>&  exercDates,   // option's exercise dates (equal settlment)
                                       const std::vector<QuantLib::Real>&  strikeValues, // strike payed at excercise dates
                                       const std::vector<QuantLib::Real>&  b76Prices,    // reference European prices
                                       // option's underlying
                                       const std::vector< std::vector<QuantLib::Time> >&  startDates,   // start dates of coupon period
                                       const std::vector< std::vector<QuantLib::Time> >&  payDates,     // pay dates of coupon perid
                                       const std::vector< std::vector<QuantLib::Real> >&  cashFlows,    // fixed coupon payments (absolut value)
                                       const std::vector< Option::Type >&                 cop,          // call (1) or put (-1) option
                                       // calibration parameters
                                       const QuantLib::Real                tol_vola );   // absolut tolerance in short rate volatility
        // Inspectors
        const std::vector< std::vector<QuantLib::Real> >& calibrationJacobian()     const { return calibrationJacobian_; }
        const std::vector< std::vector<QuantLib::Real> >& europeansAnalyticalVega() const { return europeansAnalyticalVega_; }
        const std::vector< std::vector<QuantLib::Real> >& europeansNumericalVega()  const { return europeansNumericalVega_; }
        const std::vector<QuantLib::Real>&                bermudanVega()            const { return bermudanVega_; }
    };

    typedef G2ppT<QuantLib::Time, QuantLib::Real, QuantLib::Real> RealG2ppModel;
}

#endif  /* ifndef quantlib_templatehullwhitemodel_hpp */
