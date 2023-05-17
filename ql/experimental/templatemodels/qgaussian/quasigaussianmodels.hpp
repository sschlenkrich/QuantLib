/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2015 Sebastian Schlenkrich

*/

/*! \file quasigaussianmodels.hpp
    \brief bindings of templatised quasi Gaussian models
*/


#ifndef quantlib_templatequasigaussianmodels_hpp
#define quantlib_templatequasigaussianmodels_hpp

//#include <ql/experimental/templatehullwhite/adtageo/adtageo.hpp>
//#include <ql/experimental/template/auxilliaries/MinimADVariable2.hpp>

#include <ql/types.hpp>

#include <ql/experimental/templatemodels/qgaussian/quasigaussianmodelT.hpp>
#include <ql/experimental/templatemodels/qgaussian/quasigaussianmodelabcdT.hpp>
#include <ql/experimental/templatemodels/stochasticprocessT.hpp>
#include <ql/experimental/templatemodels/qgaussian/qgswaptionmodelT.hpp>

namespace QuantLib {


    // basic binding of template parameters
    typedef QuasiGaussianModelT<QuantLib::Time,QuantLib::Real,QuantLib::Real> RealQuasiGaussianModel;
    typedef QuasiGaussianModelAbcdT<QuantLib::Time,QuantLib::Real,QuantLib::Real> RealQuasiGaussianModelAbcd;

    typedef QGSwaptionModelT<QuantLib::Time,QuantLib::Real,QuantLib::Real> RealQGSwaptionModel;

}

#endif  /* ifndef quantlib_templatequasigaussianmodels_hpp */
