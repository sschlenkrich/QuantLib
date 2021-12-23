/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2018 Sebastian Schlenkrich

*/



#ifndef quantlib_qglsvmodel_hpp
#define quantlib_qglsvmodel_hpp

#include <ql/experimental/templatemodels/qgaussian2/qglocalvolmodel.hpp>


namespace QuantLib {


    // alternative backward-looking calibration methodology for QGLocalvolModel
    class QGLSVModel : public QGLocalvolModel {
    protected:
        size_t nStrikes_;
        Real   svKernelScaling_;
    public:
        QGLSVModel(
            const Handle<YieldTermStructure>&                      termStructure,
            const Handle<SwaptionVolatilityStructure>&             volTS,
            const Real                                             chi,
            const Real                                             theta,
            const Real                                             eta,
            const ext::shared_ptr<SwapIndex>&                      swapIndex,
            const std::vector<Real>&                               times,
            const size_t                                           nStrikes,
            const bool                                             calcStochVolAdjustment,
            const Real                                             kernelWidth,       // 1.06 N^0.2, see Silverman's rule of thumb
            const Real                                             svKernelScaling,   // ~3.0, smooth conditional expectation typically requires larger kernel width
            const size_t                                           nPaths,
            const BigNatural                                       seed = 1234,
            const size_t                                           debugLevel = 1)
            : QGLocalvolModel(termStructure, volTS, chi, theta, eta, swapIndex, times, std::vector<Real>(), calcStochVolAdjustment, kernelWidth, nPaths, seed, debugLevel), nStrikes_(nStrikes), svKernelScaling_(svKernelScaling) {}
        // do the actual calculation

        virtual ~QGLSVModel() = default;
        
        virtual void simulateAndCalibrate();
    };




}

#endif  /* ifndef quantlib_qglsvmodel_hpp */
