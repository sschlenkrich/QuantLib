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

/*! \file localcorrsurfaceabfFX.hpp
    \brief Local Correlation surface derived.
	J. Guyon, 2013, A new Class of local correlation models
*/

#ifndef quantlib_localcorrsurfaceabffx_hpp
#define quantlib_localcorrsurfaceabffx_hpp

#include <ql/experimental/termstructures/localcorrsurfaceabf.hpp>

namespace QuantLib {

    //! Local Correlation surface derived 
    /*! For details about this implementation refer to
        
        \bug this class is untested, probably unreliable.
    */
    class LocalCorrSurfaceABFFX : public LocalCorrSurfaceABF {
      public:
        LocalCorrSurfaceABFFX(const std::vector<ext::shared_ptr<QuantLib::GeneralizedBlackScholesProcess>>& processes,
							  const ext::shared_ptr<QuantLib::GeneralizedBlackScholesProcess>&			  processToCal);
		LocalCorrSurfaceABFFX(const std::vector<ext::shared_ptr<QuantLib::HestonSLVProcess>>& processes,
			const ext::shared_ptr<QuantLib::GeneralizedBlackScholesProcess>&			  processToCal,
			const RealStochasticProcess::MatA											  correlations);

		//@}
		//! \name Visitability
		//@{
		virtual void accept(AcyclicVisitor&);
        //@}
		virtual QuantLib::Real localA(Time t, const RealStochasticProcess::VecA& assets,
			bool extrapolate = false) const = 0;
		virtual QuantLib::Real localB(Time t, const RealStochasticProcess::VecA& assets,
			bool extrapolate = false) const = 0;

		Matrix getLocalCorrelationSurface(Time t, std::vector<Real> assetGrid1, std::vector<Real> assetGrid2);
      protected:
		  virtual Real localCorrImplTeq0(Time t, const RealStochasticProcess::VecA& X0, bool extrapolate = false);
		  virtual QuantLib::Real localFStrike(Time t, const RealStochasticProcess::VecA& X0);
		  virtual QuantLib::Real checkLambdaValue(QuantLib::Real lambda);
	  private:
		  RealStochasticProcess::MatA getPureHestonImpliedCorrelationMatrix();
    };
}

#endif
