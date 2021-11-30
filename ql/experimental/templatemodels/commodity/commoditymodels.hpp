/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2015 Sebastian Schlenkrich

*/

/*! \file commoditymodels.hpp
    \brief bindings of templatised commodity models
*/


#ifndef quantlib_templatecommoditymodels_hpp
#define quantlib_templatecommoditymodels_hpp

#include <ql/experimental/templatemodels/commodity/twofactornormalmodelT.hpp>
#include <ql/experimental/templatemodels/commodity/twofactorlognormalmodelT.hpp>

namespace QuantLib {

    // basic binding of template parameters
    typedef TwoFactorMeanReversionModelT<QuantLib::Time,QuantLib::Real,QuantLib::Real> Real2FMeanReversionModel;

    // basic binding of template parameters
    typedef TwoFactorNormalModelT<QuantLib::Time,QuantLib::Real,QuantLib::Real> Real2FNormalModel;

    // basic binding of template parameters
    typedef TwoFactorLognormalModelT<QuantLib::Time,QuantLib::Real,QuantLib::Real> Real2FLognormalModel;

}

#endif  /* ifndef quantlib_templatecommoditymodels_hpp */
