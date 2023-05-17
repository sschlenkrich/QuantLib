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

#include <ql/experimental/templatemodels/multiasset/localcorrsurfaceabf.hpp>

namespace QuantLib {

    LocalCorrSurfaceABF::LocalCorrSurfaceABF(
		const std::vector<ext::shared_ptr<QuantLib::GeneralizedBlackScholesProcess>>& processes,
		const ext::shared_ptr<QuantLib::GeneralizedBlackScholesProcess>&			    processToCal)
    : LocalCorrTermStructure(processes, processToCal){
		assetTemp_ = std::vector<QuantLib::Real>(processes_.size());
    }

	LocalCorrSurfaceABF::LocalCorrSurfaceABF(
		const std::vector<ext::shared_ptr<QuantLib::HestonSLVProcess>>& processes,
		const ext::shared_ptr<QuantLib::GeneralizedBlackScholesProcess>&			    processToCal)
		: LocalCorrTermStructure(processes, processToCal) {
		assetTemp_ = std::vector<QuantLib::Real>(processes_.size());
	}

    void LocalCorrSurfaceABF::accept(AcyclicVisitor& v) {
        Visitor<LocalCorrSurfaceABF>* v1 =
            dynamic_cast<Visitor<LocalCorrSurfaceABF>*>(&v);
        if (v1 != 0)
            v1->visit(*this);
        else
            LocalCorrTermStructure::accept(v);
    }

	void LocalCorrSurfaceABF::localCorrImpl(RealStochasticProcess::MatA& corrMatrix, Time t, const RealStochasticProcess::VecA& X0,
		bool extrapolate) {
		Real lambda;
		QL_REQUIRE(times_.size()>0,"Local correlation model not calibrated yet.");
		if (t < times_[1]) { //first time step in grid does not depend on F,a,b.
			lambda = localCorrImplTeq0(t,X0, extrapolate);
		}
		else
		{
			for (size_t i = 0; i < assetTemp_.size(); i++)
			{
				assetTemp_[i] = processes_[i]->s0()->value()*std::exp(X0[i]);
			}
			lambda = (localF(t, X0, extrapolate) - localA(t, assetTemp_, extrapolate)) / localB(t, assetTemp_, extrapolate);
		}

		lambda = checkLambdaValue(lambda);

		for (size_t i = 0; i < corr0_.size(); i++) //corr0 because corrMatrix may have bigger size (i.e. Heston has volvol correlation)
		{
			for (size_t j = i; j < corr0_[i].size(); j++)
			{
				corrMatrix[i][j] = (1 - lambda)*corr0_[i][j] + lambda  *corr1_[i][j];
				//QL_REQUIRE(corrMatrix[i][j] != 1 || i==j, "correlation is not allowed wo be 1 for i!=j");
				corrMatrix[j][i] = corrMatrix[i][j];
			}
		}
	}

	QuantLib::Real LocalCorrSurfaceABF::localF(Time t, const RealStochasticProcess::VecA& X0,
		bool extrapolate) {
		
		double strike = localFStrike(t, X0);
		double strikeIn;
		double res;

		//first interpolate strike dimension, function F merely called after t1 (before local corr independent of a,b,F)
		for (size_t i = 1; i < times_.size(); i++)
		{
			if (strikes_[i].size() != 0) {
				//no extrapolation for now as proposed in paper....
				strikeIn = strike;
				if (strike > strikes_[i][strikes_[i].size()-1]) strikeIn = strikes_[i][strikes_[i].size()-1];
				if (strike < strikes_[i][0]) strikeIn = strikes_[i][0];

				valuesSecInt_[i] = interpolatorStrikesF_[i](strikeIn, extrapolate);
			}
			else {
				//calibration not completed, this is why size of strike array can be zero for a certain time.
				valuesSecInt_[i] = valuesSecInt_[i - 1]; //constant extrapolation during calibration
			}
		}
		interpolatorTimesF_.update();
		//second interpolate time dimension 
		res = interpolatorTimesF_(t, extrapolate);
		if (res != res)
			QL_FAIL("Interpolation on correlation surface failed. Erroneous surface.");
		return res;
	}
}

