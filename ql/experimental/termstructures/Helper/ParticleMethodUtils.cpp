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

#include <ql/experimental/termstructures/Helper/ParticleMethodUtils.hpp>
#include <ql/termstructures/yield/flatforward.hpp>
#include <ql/experimental/templatemodels/multiasset/localcorrelationSLVmodel.hpp>
#include <ql/experimental/templatemodels/montecarlo/montecarlomodells.hpp>
#include <ql/experimental/templatemodels/montecarlo/mcpayoffT.hpp>
#include <ql/math/interpolations/linearinterpolation.hpp>
#include <math.h>
#include <boost/math/distributions.hpp>
#include <boost/shared_ptr.hpp>

namespace QuantLib {

	void ParticleMethodUtils::calibrateFX(Handle<LocalCorrSurfaceABFFX> surface, const std::string& kernelIn, unsigned int numberOfPaths, Time maxTime,
		Time deltaT, Time tMin, Real kappa, Real sigmaAVR, Real exponentN, Real gridMinQuantile,
		Real gridMaxQuantile, unsigned int ns1, unsigned int ns2) {
		
		ext::shared_ptr<KernelInterface> kernel;

		if (kernelIn == "QuarticKernel") {
			kernel = ext::shared_ptr<KernelInterface>(new QuarticKernel());
		}
		else {
			QL_REQUIRE(false, "Kernel not supported. Supported is: QuarticKernel");
		}

		std::vector<ext::shared_ptr<QuantLib::HestonSLVProcess>> processes = surface->getProcesses();
		ext::shared_ptr<QuantLib::GeneralizedBlackScholesProcess>			    processToCal= surface->getProcessToCal();

		//build assetModel for simulation
		ext::shared_ptr<YieldTermStructure> yld(new FlatForward(processToCal->blackVolatility()->referenceDate(),0, processToCal->blackVolatility()->dayCounter()));
		Handle<YieldTermStructure> yldH = Handle<YieldTermStructure>(yld);
		std::vector<std::string> aliases(2);
		aliases[0] = "fx1";
		aliases[1] = "fx2";
		Handle<LocalCorrTermStructure> surfaceGen(surface.currentLink());
		ext::shared_ptr<LocalCorrelationSLVModel> assetModel = ext::shared_ptr<LocalCorrelationSLVModel>((new LocalCorrelationSLVModel(yldH, aliases, processes, surfaceGen)));

		std::vector<Time>& times = surface->getTimes();
		std::vector<std::vector<Real>>& strikes = surface->getStrikes();
		std::vector<std::vector<Real>>& surfaceF = surface->getSurfaceF();

		//time grid from t=0 to maxTime:
		times.resize(1);
		times[0] = 0;
		size_t i = 0;
		while (times[i] < maxTime) {
			times.push_back(times[i]+deltaT);
			i++;
		}
		surface->setInterpolationTime<Linear>();

		RealMCSimulation simulation(assetModel, times, times, numberOfPaths,1,true,true,false);
		simulation.prepareForSlicedSimulation();
		simulation.simulate(0, false);

		//create strike grid. 
		//the strike grid depends on simulation results (min and max quantile)
		//the simulation itself depends on local correlation which itself is calibrated using this function.
		//Consequence: strike grid can merely be calculated iteratively during calibration
		
		size_t numberStrikes;
		Real strikeStep;
		surfaceF.resize(times.size());
		strikes.resize(times.size());
		std::vector<Real> state;
		std::vector<Real> assets(2);
		std::vector<Real> crossFX(numberOfPaths);

		Real vol1;
		Real vol2;
		std::vector<Real> vol3;
		std::vector<Real> eNum;
		std::vector<Real> eDen;
		std::vector<Real> eScale;
		Real a = 0;
		Real b = 0;
		Real kernelV = 0;
		Real bandwidth = 0;
		Real minPosBw;
		Real maxNegBw;
		Real bwIn;
		Real bwRatio;

		//Calculate local correlation successively over time

		for (size_t i = 1; i < surfaceF.size()-1; i++) //iteration over time, for i=0 nothing to do as correlation independent of a,b,f. In last entry, no additional simulation necessary.
		{
			numberStrikes = numberStrikeGrid(times[i],ns1,ns2);
			QL_REQUIRE(numberStrikes>1,"ns1 or ns2 has to be increased, strike grid cannot be calculated.");
			simulation.simulate(i,false);

			//Now strike grid can be calculated

			surfaceF[i].resize(numberStrikes);
			strikes[i].resize(numberStrikes);
			for (size_t k = 0; k < numberOfPaths; k++)
			{
				state = simulation.state(k, times[i]);
				assets[0] = processes[0]->s0()->value() * std::exp(state[0]);
				assets[1] = processes[1]->s0()->value() * std::exp(state[1]);
				crossFX[k] = getCrossFX(assets[0] , assets[1]);
			}

			std::sort(crossFX.begin(),crossFX.end());

			strikes[i][0] = crossFX[(size_t)(numberOfPaths*gridMinQuantile)];
			strikes[i][strikes[i].size()-1] = crossFX[(size_t)(numberOfPaths*gridMaxQuantile)];

			strikeStep = (strikes[i][strikes[i].size() - 1] - strikes[i][0]) / (numberStrikes-1);

			bandwidth = ParticleMethodUtils::bandwidth(times[i], getCrossFX(processes[0]->s0()->value(),processes[1]->s0()->value()), kappa, sigmaAVR, tMin, numberOfPaths, exponentN);

			for (size_t j = 1; j < strikes[i].size()-1; j++)
			{
				strikes[i][j] = strikes[i][j-1] + strikeStep;
			}

			//Particle method for that time step 

			eNum.resize(strikes[i].size());
			eDen.resize(strikes[i].size());
			eScale.resize(strikes[i].size());
			vol3.resize(strikes[i].size());

			for (size_t j = 0; j < strikes[i].size(); j++)
			{
				eNum[j] = 0;
				eDen[j] = 0;
				eScale[j] = 0;
				vol3[j] = processToCal->localVolatility()->localVol(times[i], strikes[i][j], true);
			}

			for (size_t k = 0; k < numberOfPaths; k++) //particle method: over MC paths
			{
				state = simulation.state(k, times[i]);

				assets[0] = processes[0]->s0()->value() * std::exp(state[0]);
				assets[1] = processes[1]->s0()->value() * std::exp(state[1]);

				state[2] = state[2] < 0 ? 0.0001 : state[2];//heston vol might become negative (dep. on feller constant)
				state[3] = state[3] < 0 ? 0.0001 : state[3];

				vol1 = processes[0]->leverageFct()->localVol(times[i], assets[0], true)*std::sqrt(state[2]); //*ai from Heston
				vol2 = processes[1]->leverageFct()->localVol(times[i], assets[1], true)*std::sqrt(state[3]);

				a = surface->localA(times[i], assets, true);
				b = surface->localB(times[i], assets, true);

				if (vol1 != vol1 || vol2 != vol2)
					QL_FAIL("surface of asset1 oder asset2 not well defined");

				for (size_t j = 0; j < surfaceF[i].size(); j++) //iteration over strike dimension
				{

					minPosBw = 10000000;
					maxNegBw = -10000000;

					bwIn = getCrossFX(assets[0], assets[1]) - strikes[i][j];
					bwRatio = bwIn / bandwidth;

					if (bwRatio > maxNegBw && bwRatio <= 0) maxNegBw = bwRatio;
					if (bwRatio < minPosBw && bwRatio >= 0) minPosBw = bwRatio;

					kernelV = assets[1] * ParticleMethodUtils::kernel(bandwidth, bwIn, kernel);

					eNum[j] += (vol1*vol1 + vol2*vol2 + 2 * a*vol1*vol2 / b)*kernelV;
					eDen[j] += vol1*vol2*kernelV / b;
					eScale[j] += kernelV;
				}
			}
			
			for (size_t j = 0; j < surfaceF[i].size(); j++) //iteration over strike dimension
			{
				QL_REQUIRE(eScale[j] != 0, std::string("ParticleMethodUtils::calibrateFX: resulting bandwidth is too small for calibration (support: < ") + std::to_string(maxNegBw)
					+ std::string(" and > ") + std::to_string(minPosBw)
					+ std::string("). Either decrease number of MC paths or increase exponentN or kappa."));
				if (eNum[j] != eNum[j] || vol3 != vol3 || eScale[j] != eScale[j] || eDen[j] != eDen[j] || eDen[j] == 0)
					QL_FAIL("surface not well defined.");
				surfaceF[i][j] = (eNum[j] - vol3[j] * vol3[j]* eScale[j]) / (2 * eDen[j]);
			}
			//set interpolation on new dimension:
			surface->setInterpolationStrike<Linear>(i);
		}
	}	  

	Real ParticleMethodUtils::bandwidth(Time t, Real s0, Real kappa, Real sigmaAVR, Real tMin, unsigned int numberOfPaths, Real exponentN) {
		Real mult = pow(numberOfPaths, exponentN);
		return kappa*sigmaAVR* s0 * sqrt(t>tMin ? t : tMin) * mult;
	}

	Real ParticleMethodUtils::kernel(Real bandwidth, Real x, ext::shared_ptr<KernelInterface>& kernel) {
		QL_REQUIRE(bandwidth != 0, "Error in ParticleMethodUtils: bandwidth is not allowed to be zero.");
		return kernel->value(x / bandwidth) / bandwidth;
	}
	size_t ParticleMethodUtils::numberStrikeGrid(Time t, unsigned int ns1, unsigned int ns2) {
		Real numberOfStrikes = ns1*sqrt(t);
		return (numberOfStrikes > ns2) ? ((size_t)numberOfStrikes) : (ns2);
	}
	
	Real ParticleMethodUtils::getCrossFX(Real asset1, Real asset2) {
		return asset1 / asset2;
	}
}
