/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2014 Sebastian Schlenkrich

 */

/*! \file basisswapengine.hpp
    \brief Discounting swap engine for basis swaps
*/

#ifndef quantlib_basisswapengine_hpp
#define quantlib_basisswapengine_hpp

#include <ql/handle.hpp>
#include <ql/pricingengine.hpp>
#include <ql/termstructures/yieldtermstructure.hpp>
#include <ql/experimental/basismodels/basisswap.hpp>


namespace QuantLib {

    class BasisSwapEngine : public BasisSwap::engine {
    protected:
        std::vector<Handle<YieldTermStructure>> discCurves_;
        std::vector<Real>                       fxForDom_;
        boost::optional<bool>                   includeSettlementDateFlows_;
        Date                                    settlementDate_, npvDate_;
    public:
        BasisSwapEngine( 
            const std::vector<Handle<YieldTermStructure>>& discCurves,
            const std::vector<Real>&                       fxForDom,
            boost::optional<bool>                          includeSettlementDateFlows = boost::none,
            Date                                           settlementDate             = Date(),
            Date                                           npvDate                    = Date()
            );
        
        void calculate() const;
    };


}

#endif
