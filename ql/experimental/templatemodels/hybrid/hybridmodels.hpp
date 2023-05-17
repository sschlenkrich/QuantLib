/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2019 Sebastian Schlenkrich

*/

/*! \file hybridmodels.hpp
    \brief bind template types to double
*/


#ifndef quantlib_hybridmodels_hpp
#define quantlib_hybridmodels_hpp

#include <ql/experimental/templatemodels/hybrid/assetmodelT.hpp>
#include <ql/experimental/templatemodels/hybrid/hybridmodelT.hpp>
#include <ql/experimental/templatemodels/hybrid/spreadmodelT.hpp>


namespace QuantLib {

    typedef AssetModelT<QuantLib::Time, QuantLib::Real, QuantLib::Real> AssetModel;

    typedef HybridModelT<QuantLib::Time, QuantLib::Real, QuantLib::Real> HybridModel;

    typedef SpreadModelT<QuantLib::Time, QuantLib::Real, QuantLib::Real> SpreadModel;

}

#endif  /* ifndef quantlib_hybridmodelshpp */
