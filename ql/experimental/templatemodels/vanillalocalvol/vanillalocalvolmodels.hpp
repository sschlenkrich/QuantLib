/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2010 Sebastian Schlenkrich

*/



#ifndef quantlib_vanillalocalvolmodels_hpp
#define quantlib_vanillalocalvolmodels_hpp

#include <ql/types.hpp>

//#include <ql/experimental/templatemodels/stochasticprocessT.hpp>
#include <ql/experimental/templatemodels/vanillalocalvol/vanillalocalvolmodelT.hpp>

namespace QuantLib {

    typedef VanillaLocalVolModelT<QuantLib::Time, QuantLib::Real, QuantLib::Real> VanillaLocalVolModelDev;

}

#endif  /* ifndef quantlib_vanillalocalvolmodels_hpp */
