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

#include <ql/experimental/templatemodels/multiasset/localCorrFX/CTSlocalInCrossNegSkewFX.hpp>

#include <ql/experimental/templatemodels/multiasset/Helper/ParticleMethodUtils.hpp>
#include <math.h>

namespace QuantLib {

	CTSlocalInCrossNegSkewFX::CTSlocalInCrossNegSkewFX(
		const std::vector<ext::shared_ptr<QuantLib::GeneralizedBlackScholesProcess>>& processes,
		const ext::shared_ptr<QuantLib::GeneralizedBlackScholesProcess>&			    processToCal, double beta)
    : LocalCorrSurfaceABFFX(processes,processToCal), beta_(beta){
		//initializeF();
		//setInterpolation<Linear>();
	}

	CTSlocalInCrossNegSkewFX::CTSlocalInCrossNegSkewFX(
		const std::vector<ext::shared_ptr<QuantLib::HestonSLVProcess>>&				processes,
		const ext::shared_ptr<QuantLib::GeneralizedBlackScholesProcess>&			    processToCal, double beta,
		const RealStochasticProcess::MatA											    correlation)
		: LocalCorrSurfaceABFFX(processes, processToCal,correlation), beta_(beta) {
		//initializeF();
		//setInterpolation<Linear>();
	}

	QuantLib::Real CTSlocalInCrossNegSkewFX::localA(Time t, const RealStochasticProcess::VecA& assets,
		bool extrapolate) const {
		return pow(assets[0]*assets[1],beta_);
	}

	QuantLib::Real CTSlocalInCrossNegSkewFX::localB(Time t, const RealStochasticProcess::VecA& assets,
		bool extrapolate) const {
		return 1;
	}
	void CTSlocalInCrossNegSkewFX::accept(AcyclicVisitor& v) {
		Visitor<CTSlocalInCrossNegSkewFX>* v1 =
			dynamic_cast<Visitor<CTSlocalInCrossNegSkewFX>*>(&v);
		if (v1 != 0)
			v1->visit(*this);
		else
			LocalCorrSurfaceABFFX::accept(v);
	}


	//void CTSlocalInCrossNegSkewFX::initializeF() {
	//	calibratorLocalCorr_->calibrateFX(strikes_, times_, surfaceF_, processes_, processToCal_);
	//	interpolatorStrikesF_.resize(times_.size());
	//	valuesSecInt_.resize(times_.size());
	//}
}

