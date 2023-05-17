/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2015 Sebastian Schlenkrich

*/

/*! \file montecarlomodells.hpp
    \brief bindings of templatised montecarlo simulation and pricing
*/


#ifndef quantlib_montecarlomodells_hpp
#define quantlib_montecarlomodells_hpp

//#include <ql/experimental/templatehullwhite/adtageo/adtageo.hpp>
//#include <ql/experimental/template/auxilliaries/MinimADVariable2.hpp>

#include <ql/experimental/templatemodels/stochasticprocessT.hpp>
#include <ql/experimental/templatemodels/montecarlo/mcsimulationT.hpp>
#include <ql/experimental/templatemodels/montecarlo/mcpayoffT.hpp>
#include <ql/experimental/templatemodels/montecarlo/ratespayoffT.hpp>
#include <ql/experimental/templatemodels/montecarlo/commoditypayoffT.hpp>
#include <ql/experimental/templatemodels/montecarlo/amcpricerT.hpp>
#include <ql/experimental/templatemodels/auxilliaries/regressionT.hpp>


namespace QuantLib {

    // basic binding of template parameters

    typedef MCSimulationT<QuantLib::Time,QuantLib::Real,QuantLib::Real> RealMCSimulation;

    typedef MCPayoffT< QuantLib::Time,QuantLib::Real,QuantLib::Real> RealMCPayoff;

    typedef BasePayoffT< QuantLib::Time,QuantLib::Real,QuantLib::Real> RealMCBase;

    typedef RatesPayoffT< QuantLib::Time,QuantLib::Real,QuantLib::Real> RealMCRates;

    typedef CommodityPayoffT< QuantLib::Time,QuantLib::Real,QuantLib::Real> RealMCCommodity;

    typedef AMCPricerT< QuantLib::Time,QuantLib::Real,QuantLib::Real> RealAMCPricer;

    typedef RealMCPayoff::Pricer          RealMCPayoffPricer;
    typedef RealMCRates::CashFlow         RealMCCashFlow;
    typedef RealMCRates::Leg              RealMCLeg;
    typedef RealMCRates::CancellableNote  RealMCCancellableNote;

    typedef TemplateAuxilliaries::Regression<QuantLib::Real> RealRegression;

}

#endif  /* ifndef quantlib_templatequasigaussianmodels_hpp */
