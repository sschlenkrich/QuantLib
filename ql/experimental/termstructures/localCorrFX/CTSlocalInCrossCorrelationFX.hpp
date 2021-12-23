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

/*! \file CTSlocalInCrossCorrelationFX.hpp
    \brief Local Correlation surface derived ....
*/

#ifndef quantlib_localInCrossCorrelation_hpp
#define quantlib_localInCrossCorrelation_hpp

#include <ql/experimental/termstructures/localCorrFX/localcorrsurfaceabfFX.hpp>

namespace QuantLib {

    //! Local Correlation surface derived 
    /*! For details about this implementation refer to
        
        \bug this class is untested, probably unreliable.
    */
    class CTSlocalInCrossCorrelationFX : public LocalCorrSurfaceABFFX {
      public:
		  CTSlocalInCrossCorrelationFX(const std::vector<ext::shared_ptr<QuantLib::GeneralizedBlackScholesProcess>>& processes,
								  const ext::shared_ptr<QuantLib::GeneralizedBlackScholesProcess>&					processToCal);
		  CTSlocalInCrossCorrelationFX(const std::vector<ext::shared_ptr<QuantLib::HestonSLVProcess>>& processes,
			  const ext::shared_ptr<QuantLib::GeneralizedBlackScholesProcess>&						processToCal,
			  const RealStochasticProcess::MatA											    correlation);

		  //@}
		  //! \name Visitability
		  //@{
		  virtual void accept(AcyclicVisitor&);
        //@}
		  virtual QuantLib::Real localA(Time t, const RealStochasticProcess::VecA& assets,
			  bool extrapolate = false) const;
		  virtual QuantLib::Real localB(Time t, const RealStochasticProcess::VecA& assets,
			  bool extrapolate = false) const;
      protected:
		  //virtual void initializeF();
	  private:
        
    };

}

#endif
