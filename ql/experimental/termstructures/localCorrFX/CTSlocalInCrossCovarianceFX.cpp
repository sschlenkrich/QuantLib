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

#include <ql/experimental/termstructures/localCorrFX/CTSlocalInCrossCovarianceFX.hpp>

#include <ql/experimental/termstructures/Helper/ParticleMethodUtils.hpp>

namespace QuantLib {

	CTSlocalInCrossCovarianceFX::CTSlocalInCrossCovarianceFX(
		const std::vector<ext::shared_ptr<QuantLib::GeneralizedBlackScholesProcess>>& processes,
		const ext::shared_ptr<QuantLib::GeneralizedBlackScholesProcess>&			    processToCal)
    : LocalCorrSurfaceABFFX(processes,processToCal){
		//initializeF();
		//setInterpolation<Linear>();
	}

	CTSlocalInCrossCovarianceFX::CTSlocalInCrossCovarianceFX(
		const std::vector<ext::shared_ptr<QuantLib::HestonSLVProcess>>&				processes,
		const ext::shared_ptr<QuantLib::GeneralizedBlackScholesProcess>&			    processToCal,
		const RealStochasticProcess::MatA											    correlation)
		: LocalCorrSurfaceABFFX(processes, processToCal,correlation) {
		//initializeF();
		//setInterpolation<Linear>();
	}

	QuantLib::Real CTSlocalInCrossCovarianceFX::localA(Time t, const RealStochasticProcess::VecA& assets,
		bool extrapolate) const {
		return 0;
	}

	QuantLib::Real CTSlocalInCrossCovarianceFX::localB(Time t, const RealStochasticProcess::VecA& assets,
		bool extrapolate) const {
		return processes_[0]->leverageFct()->localVol(t, assets[0], true) * processes_[1]->leverageFct()->localVol(t, assets[1], true);
	}
	void CTSlocalInCrossCovarianceFX::accept(AcyclicVisitor& v) {
		Visitor<CTSlocalInCrossCovarianceFX>* v1 =
			dynamic_cast<Visitor<CTSlocalInCrossCovarianceFX>*>(&v);
		if (v1 != 0)
			v1->visit(*this);
		else
			LocalCorrSurfaceABFFX::accept(v);
	}


	//void CTSlocalInCrossCovarianceFX::initializeF() {
	//	calibratorLocalCorr_->calibrateFX(strikes_, times_, surfaceF_, processes_, processToCal_);
	//	interpolatorStrikesF_.resize(times_.size());
	//	valuesSecInt_.resize(times_.size());
	//}
}

