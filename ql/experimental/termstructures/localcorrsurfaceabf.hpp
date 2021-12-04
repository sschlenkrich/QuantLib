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

/*! \file localcorrsurface.hpp
    \brief Local Correlation surface derived ....
*/

#ifndef quantlib_localcorrsurfaceabf_hpp
#define quantlib_localcorrsurfaceabf_hpp

#include <ql/experimental/termstructures/localcorrtermstructure.hpp>
#include <ql/math/interpolation.hpp>

namespace QuantLib {

    //! Local Correlation surface derived 
    /*! For details about this implementation refer to
        
        \bug this class is untested, probably unreliable.

		J. Guyon, 2013, A new Class of local correlation models
    */
	class LocalCorrSurfaceABF : public LocalCorrTermStructure {
	public:
		LocalCorrSurfaceABF(const std::vector<ext::shared_ptr<QuantLib::GeneralizedBlackScholesProcess>>& processes,
			const ext::shared_ptr<QuantLib::GeneralizedBlackScholesProcess>&		  	    processToCal);
		LocalCorrSurfaceABF(const std::vector<ext::shared_ptr<QuantLib::HestonSLVProcess>>& processes,
			const ext::shared_ptr<QuantLib::GeneralizedBlackScholesProcess>&		  	    processToCal);


		//@}
		//! \name Visitability
		//@{
		virtual void accept(AcyclicVisitor&);

		//@}
		//@}
		//! \name Modifiers
		//@{
		template <class Interpolator>
		void setInterpolationStrike(size_t idx_time, const Interpolator& i = Interpolator()) {
			std::vector<std::vector<Real> >::iterator strike_iterator;
			std::vector<std::vector<Real> >::iterator surface_iterator = surfaceF_.begin();

			size_t j = 0;
			bool stop = false;

			for (strike_iterator = strikes_.begin();strike_iterator != strikes_.end() && ! stop;++strike_iterator, ++j, ++surface_iterator) {
				if (j == idx_time) {
					interpolatorStrikesF_[j] =
						i.interpolate((*strike_iterator).begin(), (*strike_iterator).end(), (*surface_iterator).begin());
					stop = true;
				}
			}
			notifyObservers();
		}
		template <class Interpolator>
		void setInterpolationTime(const Interpolator& i = Interpolator()) {
			valuesSecInt_.resize(times_.size());
			interpolatorStrikesF_.resize(times_.size());
			interpolatorTimesF_ = i.interpolate(times_.begin(), times_.end(), valuesSecInt_.begin());
			notifyObservers();
		}

		Date maxDate() const {
			return maxDate_;
		}

		virtual QuantLib::Real localA(Time t, const RealStochasticProcess::VecA& assets,
			bool extrapolate = false) const = 0;
		virtual QuantLib::Real localB(Time t, const RealStochasticProcess::VecA& assets,
			bool extrapolate = false) const = 0;

		std::vector<Time>& getTimes() { return times_; };
		std::vector<std::vector<Real>>& getStrikes() {return strikes_;};
		std::vector<std::vector<Real>>& getSurfaceF() { return surfaceF_; };

      protected:
		  virtual Real localCorrImplTeq0(Time t, const RealStochasticProcess::VecA& X0, bool extrapolate = false) = 0;
		  void localCorrImpl(RealStochasticProcess::MatA& corrMatrix, Time t, const RealStochasticProcess::VecA& X0,
			  bool extrapolate = false);
		  QuantLib::Real localF(Time t, const RealStochasticProcess::VecA& X0,
			  bool extrapolate = false);
		  virtual QuantLib::Real localFStrike(Time t, const RealStochasticProcess::VecA& X0) = 0;
		  virtual QuantLib::Real checkLambdaValue(QuantLib::Real lambda) = 0;

		  RealStochasticProcess::MatA corr0_;
		  RealStochasticProcess::MatA corr1_;

		  std::vector<Interpolation> interpolatorStrikesF_; //interpolated first
		  Interpolation interpolatorTimesF_; //interpolated second

		  Date maxDate_;
		  std::vector<Time> times_;
		  std::vector<std::vector<Real>> strikes_; //each time step has its own strike grid
		  std::vector<std::vector<Real>> surfaceF_;
		  std::vector<Real> valuesSecInt_;

		  std::vector<Real> assetTemp_;

	  private:
        
    };

}

#endif
