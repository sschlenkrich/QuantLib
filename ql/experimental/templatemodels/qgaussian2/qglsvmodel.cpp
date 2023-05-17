/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2018 Sebastian Schlenkrich

*/

#include <chrono>

#include <ql/math/interpolations/linearinterpolation.hpp>
#include <ql/math/interpolations/cubicinterpolation.hpp>
#include <ql/termstructures/volatility/smilesection.hpp>

#include <ql/experimental/templatemodels/qgaussian2/qglsvmodel.hpp>

namespace QuantLib {

    namespace {
        class LessByVector {
            const std::vector<Real>* vector_;
        public:
            LessByVector(const std::vector<Real>* vector) : vector_(vector) {}
            bool operator() (const size_t &a, const size_t &b) {
                return (*vector_)[a] < (*vector_)[b];
            }
        };

        inline size_t minIdxSorted(const std::vector<Real>& X, const std::vector<size_t>& sort, const Real x) {
            if (x <= X[sort[0]]) return 0;
            if (x > X[sort.back()]) return X.size();
            // bisection search
            size_t a = 0, b = X.size() - 1;
            while (b - a > 1) {
                size_t s = (a + b) / 2;
                if (x <= X[sort[s]]) b = s;
                else           a = s;
            }
            return b;
        }

        inline Real triangularKernel(const Real u) { return ((fabs(u) < 1.0) ? (1.0 - fabs(u)) : (0.0)); }
    }

    void QGLSVModel::simulateAndCalibrate() {
        std::chrono::time_point<std::chrono::steady_clock> startTime; // used for profiling
        std::map<std::string, Real> compTime;
        compTime["Samples"] = 0.0;
        compTime["Sort"]    = 0.0;
        compTime["Pricing"] = 0.0;
        compTime["Interp"]  = 0.0;
        compTime["LocVol"]  = 0.0;
        compTime["SvAdj"]   = 0.0;
        compTime["Simul"]   = 0.0;

        Initialiser init(ext::static_pointer_cast<QGLocalvolModelBackwardFlavor>(shared_from_this()));
        simulation_->simulate(1);

        QL_REQUIRE(nStrikes_ > 1, "nStrikes_ > 1 required");
        Real binSize = (1.0 * simulation_->nPaths()) / (nStrikes_ + 1.0);
        QL_REQUIRE(binSize > 1.0, "binSize > 1.0 required");
        std::vector<size_t> strikeIdx(nStrikes_);  // save the index for convenience
        for (size_t k = 0; k < strikeIdx.size(); ++k) strikeIdx[k] = (size_t)((k + 1) * binSize);

        // run actual simulation and calibration
        for (size_t idx = 1; idx < times().size(); ++idx) {
            // specify swap rate fixing at T[idx], but observed at T[idx-1]
            SwapRate swapRate(this, init.today(), times()[idx]);

            // gather swap rate samples etc. and calculate adjusters
            startTime = std::chrono::steady_clock::now();
            std::vector<Real> oneOverBSample(simulation_->nPaths(), 0.0);
            std::vector<Real> annuitySample(simulation_->nPaths(), 0.0);
            std::vector<Real> swapRateSample(simulation_->nPaths(), 0.0);
            std::vector<Real> stochVarianceSample(simulation_->nPaths(), 0.0);
            ext::shared_ptr<QGLocalvolModel::MCPayoff> mcFloatLeg(new MCAnnuity(times()[idx - 1], swapRate.scf().floatTimes(), swapRate.scf().floatWeights()));
            ext::shared_ptr<QGLocalvolModel::MCPayoff> mcFixedLeg(new MCAnnuity(times()[idx - 1], swapRate.scf().fixedTimes(), swapRate.scf().annuityWeights()));
            ext::shared_ptr<QGLocalvolModel::MCPayoff> one(new QGLocalvolModel::MCBase::FixedAmount(1.0));
            ext::shared_ptr<QGLocalvolModel::MCPayoff> oneAtT(new QGLocalvolModel::MCBase::Pay(one, times()[idx - 1]));
            for (size_t k = 0; k < simulation_->nPaths(); ++k) {
                ext::shared_ptr<MCSimulation::Path> p = simulation_->path(k);
                oneOverBSample[k] = oneAtT->discountedAt(p);
                annuitySample[k] = mcFixedLeg->at(p);
                swapRateSample[k] = mcFloatLeg->at(p) / annuitySample[k];
                stochVarianceSample[k] = simulation_->observedPath(k)[idx][2];
            }
            // calculate adjusters suksessively
            Real mcDF = 0.0;
            for (size_t k = 0; k < simulation_->nPaths(); ++k) mcDF += oneOverBSample[k];
            mcDF /= simulation_->nPaths();
            Real adjOneOverB = termStructure_->discount(times()[idx - 1]) / mcDF;
            for (size_t k = 0; k < simulation_->nPaths(); ++k) oneOverBSample[k] *= adjOneOverB;
            Real mcAnnuity = 0.0;
            for (size_t k = 0; k < simulation_->nPaths(); ++k) mcAnnuity += (annuitySample[k] * oneOverBSample[k]);
            mcAnnuity /= simulation_->nPaths();
            Real adjAnnuity = swapRate.annuity() / mcAnnuity;
            for (size_t k = 0; k < simulation_->nPaths(); ++k) annuitySample[k] *= adjAnnuity;
            Real mcFloat = 0.0;
            for (size_t k = 0; k < simulation_->nPaths(); ++k) mcFloat += (annuitySample[k] * swapRateSample[k] * oneOverBSample[k]);
            mcFloat /= simulation_->nPaths();
            Real adjSwapRate = swapRate.swapRate() - mcFloat / swapRate.annuity();
            for (size_t k = 0; k < simulation_->nPaths(); ++k) swapRateSample[k] += adjSwapRate;
            compTime["Samples"] += std::chrono::duration<double>(std::chrono::steady_clock::now() - startTime).count();

            // calculate sort index
            startTime = std::chrono::steady_clock::now();
            std::vector<size_t> sorted(simulation_->nPaths());
            for (size_t k = 0; k < simulation_->nPaths(); ++k) sorted[k] = k;
            std::sort(sorted.begin(), sorted.end(), LessByVector(&swapRateSample));
            compTime["Sort"] += std::chrono::duration<double>(std::chrono::steady_clock::now() - startTime).count();

            // determine strikes
            std::vector<Real> strikes(nStrikes_);
            for (size_t k = 0; k < strikes.size(); ++k) strikes[k] = swapRateSample[sorted[strikeIdx[k]]];

            // calculate call prices and probability function
            startTime = std::chrono::steady_clock::now();
            Real stdDev = 0.0;
            for (size_t k = 0; k < simulation_->nPaths(); k++) stdDev += swapRateSample[k] * swapRateSample[k] * annuitySample[k] * oneOverBSample[k];  // including R.-N.-derivative
            stdDev = sqrt(stdDev / simulation_->nPaths() / swapRate.annuity() - swapRate.swapRate()*swapRate.swapRate());  // Var[X] = E[X^2] - E[x]^2			
            std::vector<Real> mcCall(nStrikes_,0.0);
            std::vector<Real> mcProb(nStrikes_,0.0);   // this is only for debugging
            std::vector<Real> mcDens(nStrikes_, 0.0);
            Real sumOfRNDerivative = 0.0;
            for (size_t k=strikes.size(); k>0; --k) {
                if (k < strikes.size()) mcCall[k - 1] = mcCall[k] + sumOfRNDerivative*(strikes[k] - strikes[k - 1]);
                if (k < strikes.size()) mcProb[k - 1] = mcProb[k];
                else mcProb[k - 1] = (Real)simulation_->nPaths();
                size_t upIdx = (k==strikes.size()) ? (simulation_->nPaths()) : (strikeIdx[k]);
                for (size_t j=upIdx; j>strikeIdx[k-1]; --j) {
                    Real radonNikodymDerivative = annuitySample[sorted[j - 1]] * oneOverBSample[sorted[j - 1]] / swapRate.annuity();
                    sumOfRNDerivative += radonNikodymDerivative;
                    Real call = swapRateSample[sorted[j - 1]] - strikes[k - 1];
                    QL_REQUIRE(call >= 0.0, "call>=0.0 required.");
                    mcCall[k - 1] += radonNikodymDerivative * call;
                    mcProb[k - 1] -= radonNikodymDerivative;
                }
                // calculate density including change of meassure
                Real h = kernelWidth_ * stdDev;  // effective kernel width for density calculation
                size_t loIdxDens = minIdxSorted(swapRateSample, sorted, strikes[k - 1] - h);
                size_t upIdxDens = minIdxSorted(swapRateSample, sorted, strikes[k - 1] + h);
                for (size_t j = loIdxDens; j < upIdxDens; ++j) {
                    Real radonNikodymDerivative = annuitySample[sorted[j]] * oneOverBSample[sorted[j]] / swapRate.annuity();
                    mcDens[k - 1] += radonNikodymDerivative * triangularKernel((swapRateSample[sorted[j]] - strikes[k - 1]) / h) / h;
                }
            }
            for (size_t k = 0; k < strikes.size(); ++k) mcCall[k] /= simulation_->nPaths();
            for (size_t k = 0; k < strikes.size(); ++k) mcProb[k] /= simulation_->nPaths();
            for (size_t k = 0; k < strikes.size(); ++k) mcDens[k] /= simulation_->nPaths();
            compTime["Pricing"] += std::chrono::duration<double>(std::chrono::steady_clock::now() - startTime).count();

            // calculate density via monotonic interpolaion
            startTime = std::chrono::steady_clock::now();
            Interpolation probInterpolation = AkimaCubicInterpolation(strikes.begin(), strikes.end(), mcProb.begin());
            compTime["Interp"] += std::chrono::duration<double>(std::chrono::steady_clock::now() - startTime).count();

            // calculate local volatility
            startTime = std::chrono::steady_clock::now();
            std::vector<Real> smileCall(nStrikes_, 0.0);
            std::vector<Real> dCalldT(nStrikes_, 0.0);
            std::vector<Real> d2CalldK2(nStrikes_, 0.0);
            std::vector<Real> dProbdK(nStrikes_, 0.0);  // debugging and benchmarking
            std::vector<Real> tmpStrikes;
            std::vector<Real> tmpLocalVol;
            ext::shared_ptr<SmileSection> smileSection = volTS_->smileSection(times()[idx], swapIndex_->tenor(), true);
            for (size_t k = 0; k < strikes.size(); ++k) {
                smileCall[k] = smileSection->optionPrice(strikes[k], Option::Call);
                dCalldT[k] = (smileCall[k] - mcCall[k]) / (times()[idx] - times()[idx - 1]);
                dProbdK[k] = probInterpolation.derivative(strikes[k]); 
                d2CalldK2[k] = mcDens[k];
                if (dCalldT[k] < 1.0e-8) {
                    if (debugLevel_ > 1) debugLog_.push_back("Warning! Skip local vol at T: " + std::to_string(times()[idx - 1]) + ", k: " + std::to_string(k) + ", K: " + std::to_string(strikes[k]) + ", mcCall: " + std::to_string(mcCall[k]) + ", smileCall: " + std::to_string(smileCall[k]));
                    continue;
                }
                if (d2CalldK2[k] < 1.0e-8) {
                    if (debugLevel_ > 1) debugLog_.push_back("Warning! Skip local vol at T: " + std::to_string(times()[idx - 1]) + ", k: " + std::to_string(k) + ", K: " + std::to_string(strikes[k]) + ", d2CalldK2: " + std::to_string(d2CalldK2[k]));
                    continue;
                }
                QL_REQUIRE(dCalldT[k] > 0, "dCalldT[k] > 0 required");
                QL_REQUIRE(d2CalldK2[k] > 0, "d2CalldK2[k] > 0 required");
                tmpStrikes.push_back(strikes[k]);
                tmpLocalVol.push_back(sqrt(2 * dCalldT[k] / d2CalldK2[k]));
            }
            Interpolation tmpVolInterp = Linear().interpolate(tmpStrikes.begin(), tmpStrikes.end(), tmpLocalVol.begin());
            std::vector<Real> localVol(nStrikes_, 0.0);
            for (size_t k = 0; k < strikes.size(); ++k) localVol[k] = (tmpVolInterp(strikes[k], true) < 1.0e-4) ? (1.0e-4) : (tmpVolInterp(strikes[k], true));
            compTime["LocVol"] += std::chrono::duration<double>(std::chrono::steady_clock::now() - startTime).count();

            // calculate conditional expectation of variance
            startTime = std::chrono::steady_clock::now();
            std::vector<Real> expectationZCondS(nStrikes_, 1.0);
            if (calcStochVolAdjustment_) {
                // E^A[.] = E^Q[ z(T)*q(T) | S(T)=K ] / E^Q[ q(T) | S(T)=K ]
                // E^Q[.] = sum{ z_i*q_i*Kernel(S_i) } / n
                // R.-N.-Derivative q(T) = N(0)/An(0) * An(T)/N(T)
                std::vector<Real> zTimesQGrid(nStrikes_, 0.0);  // numerator
                std::vector<Real> qGrid(nStrikes_, 0.0);  // denumerator
                for (size_t k = 0; k < strikes.size(); ++k) {
                    Real h = svKernelScaling_ * kernelWidth_ * stdDev;  // effective kernel width for density calculation
                    size_t loIdx = minIdxSorted(swapRateSample, sorted, strikes[k] - h);
                    size_t upIdx = minIdxSorted(swapRateSample, sorted, strikes[k] + h);
                    for (size_t j = loIdx; j < upIdx; ++j) {
                        Real radonNikodymDerivative = annuitySample[sorted[j]] * oneOverBSample[sorted[j]] / swapRate.annuity();
                        Real kernelAtStrike = triangularKernel((swapRateSample[sorted[j]] - strikes[k])/ h) / h;
                        zTimesQGrid[k] += stochVarianceSample[sorted[j]] * radonNikodymDerivative * kernelAtStrike;
                        qGrid[k] += radonNikodymDerivative * kernelAtStrike;
                    }
                    expectationZCondS[k] = zTimesQGrid[k] / qGrid[k];
                    QL_REQUIRE(zTimesQGrid[k] > 0.0, "zTimesQGrid[k] > 0.0 required");
                    QL_REQUIRE(qGrid[k]> 0.0, "qGrid[k] > 0.0 required");
                }
            }
            compTime["SvAdj"] += std::chrono::duration<double>(std::chrono::steady_clock::now() - startTime).count();

            // adjust local volatility
            for (size_t k = 0; k < strikes.size(); ++k) localVol[k] /= sqrt(expectationZCondS[k]);

            // set up interpolation
            strikeGrid_.push_back(strikes);
            locvolGrid_.push_back(localVol);
            sigmaS_.push_back(Linear().interpolate(strikeGrid_.back().begin(), strikeGrid_.back().end(), locvolGrid_.back().begin()));

            // set up swap rate model here such that we don't need to care about it during sigma_x calculation
            swapRateModel_ = qGSwapRateModel(swapRate.scf(), times()[idx - 1]);

            // simulate the next step T[idx-1] to T[idx]
            startTime = std::chrono::steady_clock::now();
            simulation_->simulate(idx + 1);
            compTime["Simul"] += std::chrono::duration<double>(std::chrono::steady_clock::now() - startTime).count();

            if (debugLevel_ > 0) debugLog_.push_back("T = " + std::to_string(times()[idx-1]) + ", fixingDate = " + std::to_string(swapRate.fixingDate().serialNumber()) + ", swapRate = " + std::to_string(swapRate.swapRate()) + ", annuity = " + std::to_string(swapRate.annuity()));
            if (debugLevel_ > 1) for (size_t k = 0; k < strikes.size(); ++k) debugLog_.push_back("T = " + std::to_string(times()[idx - 1]) + ", k: " + std::to_string(k) + ", K: " + std::to_string(strikes[k]) + ", dCdT: " + std::to_string(dCalldT[k]) + ", d2CdK2: " + std::to_string(d2CalldK2[k]) + ", dProbdK: " + std::to_string(dProbdK[k]) + ", lVol: " + std::to_string(localVol[k]) + ", E[z|S]: " + std::to_string(expectationZCondS[k]));
        }
        if (debugLevel_ > 0) debugLog_.push_back("Profiling (in sec.) Samples: " + std::to_string(compTime["Samples"]) + ", Sort: " + std::to_string(compTime["Sort"]) + ", Pricing: " + std::to_string(compTime["Pricing"]) + ", Interp: " + std::to_string(compTime["Interp"]) + ", LocVol: " + std::to_string(compTime["LocVol"]) + ", SvAdj: " + std::to_string(compTime["SvAdj"]) + ", Simul: " + std::to_string(compTime["Simul"]));
    }



}

