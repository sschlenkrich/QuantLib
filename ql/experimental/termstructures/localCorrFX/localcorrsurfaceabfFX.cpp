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

#include <vector>
#include <ql/errors.hpp>
#include <ql/experimental/templatemodels/auxilliaries/choleskyfactorisationT.hpp>
#include <ql/experimental/termstructures/localCorrFX/localcorrsurfaceabfFX.hpp>
#include <ql/experimental/termstructures/Helper/ParticleMethodUtils.hpp>
#include <ql/math/matrixutilities/pseudosqrt.hpp>
#include <ql/experimental/templatemodels/auxilliaries/svdT.hpp>
#include <ql/math/matrixutilities/symmetricschurdecomposition.hpp>

namespace QuantLib {

	namespace
	{

		// Matrix infinity norm. See Golub and van Loan (2.3.10) or
		// <http://en.wikipedia.org/wiki/Matrix_norm>
		Real normInf(const Matrix& M) {
			Size rows = M.rows();
			Size cols = M.columns();
			Real norm = 0.0;
			for (Size i = 0; i<rows; ++i) {
				Real colSum = 0.0;
				for (Size j = 0; j<cols; ++j)
					colSum += std::fabs(M[i][j]);
				norm = std::max(norm, colSum);
			}
			return norm;
		}

		//J. Higham, Computating the nearest correlation matrix - a problem from finance
		void projectSymmetricToCorrelation(RealStochasticProcess::MatA& Y)
		{
			//Set W:=I and compute Algorithm 3.3

			RealStochasticProcess::MatA S = RealStochasticProcess::MatA(Y);
			Matrix R = Matrix(Y.size(), Y.size());
			Matrix X = Matrix(Y.size(), Y.size());
			Matrix Xp = Matrix(Y.size(), Y.size());
			Matrix SV = Matrix(Y.size(), Y.size(), 0.0);

			int cc = 0;
			Real tol = 10e-8;
			Real err = QL_MAX_REAL;

			while (err > tol)
			{
				for (size_t i = 0; i < Y.size(); i++)
				{
					for (size_t j = 0; j < Y.size(); j++)
					{
						if (cc == 0) S[i][j] = 0; //initialize
						R[i][j] = Y[i][j] - S[i][j]; //Dykstra's correction
						Xp[i][j] = X[i][j];
					}
				}

				//projectino onto S
				//SVD tmp = SVD(R);
				
				SymmetricSchurDecomposition dec = SymmetricSchurDecomposition(R);

				for (size_t i = 0; i < Y.size(); i++)
				{
					SV[i][i] = std::max(dec.eigenvalues()[i], 0.0);
				}
				
				X = dec.eigenvectors()*SV*transpose(dec.eigenvectors());

				for (size_t i = 0; i < Y.size(); i++)
				{
					for (size_t j = 0; j < Y.size(); j++)
					{
						QL_ASSERT(abs(X[i][j] - X[j][i]) < 10e-8, "X not symmetric.");
						//X[i][j] = X[j][i]; //we need QL_EPSILON accuracy
						Y[i][j] = i == j ? 1 : X[i][j]; //projektion onto U
						S[i][j] = X[i][j] - R[i][j];
					}
				}
				cc++;
				QL_ASSERT(cc < 100, "Correlation matrix projection does not converge.");
				err = normInf(X - Xp) / normInf(X);
			}
			//check
			TemplateAuxilliaries::performCholesky(Y, Y.size(),true);
		}
		
	}


    LocalCorrSurfaceABFFX::LocalCorrSurfaceABFFX(
		const std::vector<ext::shared_ptr<QuantLib::GeneralizedBlackScholesProcess>>& processes,
		const ext::shared_ptr<QuantLib::GeneralizedBlackScholesProcess>&  		    processToCal)
    : LocalCorrSurfaceABF(processes, processToCal){
		corr0_ = RealStochasticProcess::MatA(2);
		corr1_ = RealStochasticProcess::MatA(2);

		for (size_t k = 0; k<2; ++k) corr0_[k].resize(2);
		for (size_t k = 0; k<2; ++k) corr1_[k].resize(2);

		//correlation is determined by lambda in abf-class, therefore simple correlation matrices

		corr0_[0][0] = 1;
		corr0_[1][0] = 0;
		corr0_[0][1] = 0;
		corr0_[1][1] = 1;

		corr1_[0][0] = 1;
		corr1_[1][0] = 1;
		corr1_[0][1] = 1;
		corr1_[1][1] = 1;
    }

	LocalCorrSurfaceABFFX::LocalCorrSurfaceABFFX(
		const std::vector<ext::shared_ptr<QuantLib::HestonSLVProcess>>&				processes,
		const ext::shared_ptr<QuantLib::GeneralizedBlackScholesProcess>&  		    processToCal,
		const RealStochasticProcess::MatA											    correlations)
		: LocalCorrSurfaceABF(processes, processToCal) {
		
		if (correlations.size() == 0) {
			corr0_ = getPureHestonImpliedCorrelationMatrix();
			corr1_ = getPureHestonImpliedCorrelationMatrix();
		}
		else {
			corr0_ = RealStochasticProcess::MatA(correlations);
			corr1_ = RealStochasticProcess::MatA(correlations);
		}
		corr0_[0][1] = -1;
		corr0_[1][0] = -1;
		corr1_[0][1] = 1;
		corr1_[1][0] = 1;
		projectSymmetricToCorrelation(corr0_);
		projectSymmetricToCorrelation(corr1_);
	}


	QuantLib::Real LocalCorrSurfaceABFFX::localFStrike(Time t, const RealStochasticProcess::VecA& X0) {
		QL_REQUIRE(X0.size() == 4 || X0.size() == 2, "Local Correlation for FX only works for two dimensional FX model.");//=4 due to heston

		return ParticleMethodUtils::getCrossFX(processes_[0]->s0()->value() * std::exp(X0[0]) , (processes_[1]->s0()->value() * std::exp(X0[1])));
		
	}
	
	void LocalCorrSurfaceABFFX::accept(AcyclicVisitor& v) {
		Visitor<LocalCorrSurfaceABFFX>* v1 =
			dynamic_cast<Visitor<LocalCorrSurfaceABFFX>*>(&v);
		if (v1 != 0)
			v1->visit(*this);
		else
			LocalCorrSurfaceABF::accept(v);
	}

	Real LocalCorrSurfaceABFFX::localCorrImplTeq0(Time t, const RealStochasticProcess::VecA& X0, bool extrapolate) {
		
		//smiled surface will through an error, therefore assume one minute ahead
		t = 1.0 / (365 * 24 * 60);
		
		Real s1 = processes_[0]->s0()->value() * std::exp(X0[0]);
		Real s2 = processes_[1]->s0()->value() * std::exp(X0[1]);
		Real vol1 = processes_[0]->leverageFct()->localVol(t, s1, extrapolate);
		Real vol2 = processes_[1]->leverageFct()->localVol(t, s2, extrapolate);
		
		if (!processes_[0]->isLocalVolProcess()) vol1 *= (X0[2] <= 0 ? 0.001 : std::sqrt(X0[2]));
		if (!processes_[1]->isLocalVolProcess()) vol2 *= (X0[3] <= 0 ? 0.001 : std::sqrt(X0[3]));

		if (vol1 != vol1) QL_FAIL("leverage function of asset 1 does have non-real values");
		if (vol2 != vol2) QL_FAIL("leverage function of asset 1 does have non-real values");

		Real vol3 = processToCal_->localVolatility()->localVol(t, ParticleMethodUtils::getCrossFX(s1 , s2), extrapolate);

		return (vol1*vol1 + vol2*vol2 - vol3*vol3) / (2 * vol1*vol2);
	}
	

	Matrix LocalCorrSurfaceABFFX::getLocalCorrelationSurface(Time t,
		std::vector<Real> assetGrid1, std::vector<Real> assetGrid2) {
		Matrix result(assetGrid1.size(), assetGrid2.size());
		std::vector<std::vector<Real>> corrM;
		if(processes_[0]->isLocalVolProcess())
			corrM = std::vector<std::vector<Real>>(2, std::vector<Real>(2));
		else
			corrM = std::vector<std::vector<Real>>(4, std::vector<Real>(4));

		std::vector<Real> x0(2);

		for (size_t i = 0; i < assetGrid1.size(); i++)
		{
			for (size_t j = 0; j < assetGrid2.size(); j++)
			{
				x0[0] = log(assetGrid1[i] / processes_[0]->s0()->value());
				x0[1] = log(assetGrid2[j] / processes_[1]->s0()->value());

				localCorr(corrM, t, x0, true);
				result[i][j] = corrM[0][1];
			}
		}
		return result;
	}

	RealStochasticProcess::MatA LocalCorrSurfaceABFFX::getPureHestonImpliedCorrelationMatrix()
	{
		RealStochasticProcess::MatA corrM = RealStochasticProcess::MatA(2 * processes_.size());

		for (size_t k = 0; k<corrM.size(); ++k) corrM[k].resize(2 * processes_.size());

		for (size_t i = 0; i < 2 * processes_.size(); i++)
		{
			for (size_t j = 0; j < 2 * processes_.size(); j++)
			{
				if (i == j) {
					corrM[i][j] = 1;
				}
				else if (i == j + processes_.size() || i + processes_.size() == j) {
					int assetIndex = std::min(i, j);
					corrM[i][j] = processes_[assetIndex]->rho();
				}
				else {
					corrM[i][j] = 0;
				}
			}
		}

		return corrM;
	}

	QuantLib::Real LocalCorrSurfaceABFFX::checkLambdaValue(QuantLib::Real lambda) {
	
		if (lambda != lambda)
			QL_FAIL("lambda is erroneous.");

		if (corr0_.size() == 2) {
			//BSprocess, range [-1,1]
			if (lambda > 1) return 0.999;
			if (lambda < -1) return -0.999;
			return lambda;
		}
		else if(corr0_.size() == 4) {
			//transform lambda from correlation to weight from 0 to 1
			//rho = (1-lambda) *corr0 + lambda*corr1
			//=>lambda = (rho-corr0)/(corr1-corr0)
			lambda = (lambda - corr0_[0][1]) / (corr1_[0][1] - corr0_[0][1]);
			if (lambda > 1) return 1;
			if (lambda < 0) return 0;
			return lambda;
		}
		else {
			QL_FAIL("Unexpected error. FXmodel should always be 2 dimensional.");
		}
	
	}

}

