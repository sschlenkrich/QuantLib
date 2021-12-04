/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2017 Cord Harms

 This file is part of QuantLib, a free-software/open-source library
 for financial quantitative analysts and developers - http://quantlib.org/

 QuantLib is free software: you can redistribute it and/or modify it
 under the terms of the QuantLib license.  You should have received a
 copy of the license along with this program; if not, please email
 <quantlib-dev@lists.sf.net>. The license is also available online at
 <http://quantlib.org/license.shtml>.

 This program is distributed in the hope that it will be useful, but WITHOUT
 ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 FOR A PARTICULAR PURPOSE.  See the license for more details.
*/

/*! \file ParticleMethodUtils.hpp
    \brief particle method for local correlation derived ....
*/

#ifndef quantlib_particleMethodUtils_hpp
#define quantlib_particleMethodUtils_hpp

#include <ql/experimental/termstructures/Helper/KernelInterface.hpp>
#include <ql/experimental/termstructures/localCorrFX/localcorrsurfaceabfFX.hpp>
#include <ql/handle.hpp>
#include <string.h>

namespace QuantLib {
  
    /*! 
			Calibration of a correlation termstructure LocalCorrSurfaceABFFX using the Particle Method (Monte Carlo calibration)
			cf. J. Guyon, 2013, A new Class of local correlation models
	*/
    class ParticleMethodUtils {
      public:

		  static void calibrateFX(Handle<LocalCorrSurfaceABFFX> surface,const std::string& kernelIn, unsigned int numberOfPaths, Time maxTime,
			  Time deltaT, Time tMin, Real kappa, Real sigmaAVR, Real exponentN, Real gridMinQuantile,
			  Real gridMaxQuantile, unsigned int ns1, unsigned int ns2);
		  static Real getCrossFX(Real asset1, Real asset2);
      protected:
		  
	  private:
		  static Real bandwidth(Time t, Real s0, Real kappa, Real sigmaAVR, Real tMin, unsigned int numberOfPaths, Real exponentN);
		  static Real kernel(Real bandwidth, Real x, ext::shared_ptr<KernelInterface>& kernel);
		  static size_t numberStrikeGrid(Time t, unsigned int ns1, unsigned int ns2);
		  
    };

}

#endif
