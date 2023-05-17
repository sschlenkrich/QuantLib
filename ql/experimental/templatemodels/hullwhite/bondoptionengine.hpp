/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2010 Sebastian Schlenkrich
*/

/*! \file bondoptionengine.hpp
    \brief engine for bond options with Hull White model
*/

#ifndef quantlib_bond_option_engine_hpp
#define quantlib_bond_option_engine_hpp


#include <ql/experimental/templatemodels/hullwhite/fixedratebondoption.hpp>
#include <ql/experimental/templatemodels/hullwhite/hullwhitemodels.hpp>

#include <ql/indexes/iborindex.hpp>
#include <ql/termstructures/yieldtermstructure.hpp>
#include <ql/handle.hpp>
#include <ql/exercise.hpp>
#include <ql/instruments/swaption.hpp>
#include <ql/termstructures/volatility/swaption/swaptionvolstructure.hpp>
#include <ql/pricingengines/swaption/blackswaptionengine.hpp>

namespace QuantLib {

    class BondOptionEngine : public FixedRateBondOption::engine {
    private:
        ext::shared_ptr<RealHullWhiteModel>        model_;                // Hull White model used
        // discretisation parameters for the numerical solution of Bermudan bond options
        Size                                       dimension_;            // discretisation of numerical solution
        Real                                       gridRadius_;           // radius of short rate grid
        Real                                       bermudanTolerance_;    // tolerance for numerical integration
        // calibration instruments
        std::vector< ext::shared_ptr<Swaption> >   referenceSwaptions_;

        // utility function to compare swaptions
        static bool lessByExerciseFirstDate (  ext::shared_ptr<Swaption> a,  ext::shared_ptr<Swaption> b) { return a->exercise()->date(0) < b->exercise()->date(0); }

    public:

        // constructor with given model and no calibration
        BondOptionEngine( const ext::shared_ptr<RealHullWhiteModel>&    model,
                          const Size                                    dimension,
                          const Real                                    gridRadius,
                          const Real                                    bermudanTolerance )
                          : model_(model), dimension_(dimension), gridRadius_(gridRadius), bermudanTolerance_(bermudanTolerance) { }

        void calculate() const;

        const ext::shared_ptr<RealHullWhiteModel>& model() const { return model_; }

        // calibrate model based on given swaptions
        void calibrateModel( std::vector< ext::shared_ptr<Swaption> >          swaptions,
                             const bool                                        contTenorSpread,
                             const Real                                        tolVola);

    };

}

#endif
