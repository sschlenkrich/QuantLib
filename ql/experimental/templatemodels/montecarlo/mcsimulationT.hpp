/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2015 Sebastian Schlenkrich

*/

/*! \file mcsimulationT.hpp
    \brief simulate and store paths of

    Process dX(t) = a(t,X) dt + b(t,X) dW
    
    Methods drift(t,X) and diffusion(t,X) are provided by ProcessType
    
*/


#ifndef quantlib_templatemcsimulation_hpp
#define quantlib_templatemcsimulation_hpp

#include <ql/math/randomnumbers/rngtraits.hpp>
#include <ql/experimental/templatemodels/stochasticprocessT.hpp>

#include <ql/experimental/templatemodels/auxilliaries/auxilliariesT.hpp>



namespace QuantLib {

    class TemplateSimulation : public virtual Observable { };

    // Declaration of MC simulation class
    template <class DateType, class PassiveType, class ActiveType>
    class MCSimulationT : public TemplateSimulation {
    private:
        class LessByVector {
            const std::vector<Real>& vector_;
        public:
            LessByVector(const std::vector<Real>& vector) : vector_(vector) {}
            bool operator() (const size_t &a, const size_t &b) {
                return vector_[a] < vector_[b];
            }
        };

    protected:

        typedef StochasticProcessT<DateType, PassiveType, ActiveType> ProcessType;

        // container class definitions
        typedef std::vector<DateType>                      VecD;
        typedef std::vector<PassiveType>                   VecP; 
        typedef std::vector<ActiveType>                    VecA;
        typedef std::vector< std::vector<DateType> >       MatD;
        typedef std::vector< std::vector<PassiveType> >    MatP;
        typedef std::vector< std::vector<ActiveType> >     MatA;

        // default pseudo random numbers
        ext::shared_ptr<PseudoRandom::rsg_type> rsg_;

        // the link to the process/model
        ext::shared_ptr<ProcessType>            process_;

        // time grid for process simulation
        VecD simTimes_;
        // time grid for process observation, subset of simTimes (store state at obsTimes_)
        VecD obsTimes_;
        // Monte Carlo properties
        BigNatural  seed_;
        bool        richardsonExtrapolation_;

        // allow time interpolation on path
        bool timeInterpolation_;

        // we may precompute and store the Brownian increments
        bool storeBrownians_;
        // Brownian increments dW (nPaths x k*nSimTimes x nFactors), k = 2 for Richardson Extrapolation
        std::vector< MatD > dW_;

        // paths stored X_[paths][obsTimes][size]
        std::vector< MatA > X_;

        // additional states calculated via interpolation; (*Y_[addObsTimes])[paths][size]
        VecD addObsTimes_;
        std::vector< ext::shared_ptr< MatA > > Y_;

        // adjust numeraire by exp{adj(t)*t} to meet E{ 1/N(t) } = P(0,t)
        bool applyNumeraireAdjuster_;
        VecD numeraireObservTimes_;
        VecA numeraireAdjuster_;

        // adjust simulated zcb by exp{-adj(t,t+dt)(t+dt)} to meet initial yield curve
        bool applyZcbAdjuster_;
        VecD zcbObservTimes_;
        VecD zcbOffsetTimes_;
        MatA zcbAdjuster_;

        // adjust simulated asset to S(t) + adj(t) to meet E[S(t)]
        bool applyAssetAdjuster_;
        std::map<std::string, size_t>  assetIndex_;  // a list of aliases for which adjusters are calculated
        VecD assetObservTimes_;                      // we use a common set of time grid points
        MatA assetAdjuster_;                         // a vector of adjusters for each alias

        // dummy
        PassiveType tmp_;

        // initialisation

        inline void reallocateMemory(size_t nPaths) {
            // clear/free memory
            dW_.clear();
            X_.clear();
            // dW_
            if (storeBrownians_) {
                dW_.resize(nPaths);
            //  this does not help because it is overwritten in getNextBrownianIncrements
            //	size_t ndWt = ((richardsonExtrapolation_) ? 2 : 1) * (simTimes_.size()-1);
            //    for (size_t i=0; i<nPaths; ++i) {
            //		dW_[i].resize(ndWt);
            //		for (size_t j=0; j<ndWt; ++j) dW_[i][j].resize(process_->factors());
            //	}
            }
            else {
                dW_.resize(1); // we only save the most recent path
            }
            // X_
            X_.resize(nPaths);
            for (size_t i=0; i<nPaths; ++i) {
                X_[i].resize(obsTimes_.size());
                for (size_t j=0; j<obsTimes_.size(); ++j) X_[i][j].resize(process_->size());
            }
            return;
        }

        inline void initialiseRSG() { // this is applied for path-wise simulation
            size_t ndWt = ((richardsonExtrapolation_) ? 2 : 1) * (simTimes_.size()-1);
            rsg_ = ext::shared_ptr<PseudoRandom::rsg_type>(new PseudoRandom::rsg_type(PseudoRandom::make_sequence_generator(ndWt * process_->factors(), seed_)));
        }

        // return the next Brownian motion path increments from random sequence generator
        inline const std::vector< std::vector<QuantLib::Real> > getNextBrownianIncrements() {
            std::vector<QuantLib::Real> sequence = rsg_->nextSequence().value;
            size_t ndWt = ((richardsonExtrapolation_) ? 2 : 1) * (simTimes_.size()-1);
            size_t nFac = process_->factors();
            QL_REQUIRE(sequence.size()==ndWt*nFac,"TemplateMCSimulation: dimension mismatch");
            std::vector< std::vector<QuantLib::Real> > dW(ndWt);
            for (size_t i=0; i<ndWt; ++i) {
                dW[i].resize(nFac);
                for (size_t j=0; j<nFac; ++j) {
                    dW[i][j] = sequence[i*nFac+j];
                }
            }
            return dW;
        }

        // return the Brownian motion path increments from from cache or random sequence generator
        inline const std::vector< std::vector<QuantLib::Real> >& getBrownianIncrementsByPath(const size_t path) {
            if (storeBrownians_) {  // in this case we can avoid data copy
                QL_REQUIRE(path<dW_.size(), "TemplateMCSimulation: path index out of bounds.");
                return dW_[path];
                }
            else {
                dW_[0] = getNextBrownianIncrements();
                return dW_[0];
            }
        }


        // simulate a single path X_[path]
        inline void simulatePath(const size_t path) {
            QL_REQUIRE(path<X_.size(),"TemplateMCSimulation: path index out of bounds.");			
            const MatD& dWt = getBrownianIncrementsByPath(path);
            // initialisation
            VecA X0 = process_->initialValues();
            VecA X1(X0.size()), X12(X0.size());
            VecD dW( process_->factors() );
            X_[path][0]    = X0;
            size_t obs_idx = 1;
            for (size_t sim_idx=1; sim_idx<simTimes_.size(); ++sim_idx) {
                DateType dt = simTimes_[sim_idx]-simTimes_[sim_idx-1];
                if (richardsonExtrapolation_) {
                    // full Euler step
                    for (size_t k=0; k<dW.size(); ++k) dW[k] = (dWt[2*(sim_idx-1)][k] + dWt[2*(sim_idx-1)+1][k])/sqrt(2.0);
                    //X1 = evolve(simTimes_[sim_idx-1],X0,dt,dW);
                    process_->evolve(simTimes_[sim_idx-1],X0,dt,dW,X1);
                    // two half size Euler steps
                    //VecA X12 = evolve(simTimes_[sim_idx-1],X0,dt/2.0, dWt[2*(sim_idx-1)]);
                    //X0 = evolve(simTimes_[sim_idx-1]+dt/2.0,X12,dt/2.0,dWt[2*(sim_idx-1)+1]);
                    process_->evolve(simTimes_[sim_idx-1],        X0,  dt/2.0, dWt[2*(sim_idx-1)],   X12);
                    process_->evolve(simTimes_[sim_idx-1]+dt/2.0, X12, dt/2.0, dWt[2*(sim_idx-1)+1], X0);
                    // extrapolation
                    for (size_t k=0; k<X1.size(); ++k) X1[k] = 2*X0[k] - X1[k];
                    // extrapolation may lead to ill-defined states
                    process_->truncate( simTimes_[sim_idx], X1 );
                } else { // only full Euler step
                    process_->evolve(simTimes_[sim_idx-1],X0,dt,dWt[sim_idx-1], X1);
                }
                if ((obs_idx<obsTimes_.size()) && (simTimes_[sim_idx]==obsTimes_[obs_idx])) {
                    X_[path][obs_idx]=X1;
                    ++obs_idx;
                }
                X0 = X1;
            }
        }


    public:

        // a path hides the actual state implementation from the payoff
        // thus payoffs only need to know the path interface and does
        // not need to care about the actual simulation
        class Path {
        protected:
            ext::shared_ptr<ProcessType>    process_;
            MCSimulationT*                    sim_;
            size_t                            idx_;
        public:
            Path(const ext::shared_ptr<ProcessType>           process,
                 MCSimulationT*                                 sim,
                 const size_t                                   idx )
                 : process_(process), sim_(sim), idx_(idx) {}

            // maybe better cache the last state to avoid repeated interpolation in state...

            // we delegate (and decorate) the implementation to the stochastic process

            inline ActiveType numeraire(DateType obsTime) {
                return sim_->numeraireAdjuster(obsTime)*
                    process_->numeraire(obsTime,sim_->state(idx_,obsTime));
            }

            inline ActiveType zeroBond(DateType obsTime, DateType payTime) {
                return sim_->zcbAdjuster(obsTime, payTime) *
                       process_->zeroBond(obsTime, payTime, sim_->state(idx_,obsTime) );
            }

            inline ActiveType zeroBond(DateType obsTime, DateType payTime, const std::string& alias) {
                return sim_->zcbAdjuster(obsTime, payTime, alias) *
                    process_->zeroBond(obsTime, payTime, sim_->state(idx_, obsTime), alias);
            }

            inline ActiveType asset(DateType obsTime, const std::string& alias) {
                    return sim_->assetAdjuster(obsTime,alias) + 
                        process_->asset(obsTime, sim_->state(idx_, obsTime), alias);
            }

            // calculate barrier no-hit probability via Brownian Bridge applied to log-asset
            inline ActiveType assetBarrierNoHit(DateType tStart,
                                                DateType tEnd,
                                                PassiveType downBarrier,
                                                PassiveType upBarrier,
                                                PassiveType downOrUpOrBoth,  // -1 (down) or +1 (up) or 0 (both)
                                                const std::string& alias) {
                // we need to find the grid on which we may calculate the
                std::vector<DateType> times;
                times.push_back(tStart);
                for (size_t k = 0; k < sim_->obsTimes().size(); ++k) {
                    if ((sim_->obsTimes()[k] > tStart) && (sim_->obsTimes()[k] < tEnd)) times.push_back(sim_->obsTimes()[k]);
                    if (sim_->obsTimes()[k] >= tEnd) break;
                }
                times.push_back(tEnd);
                // we initialise with the asset at start date
                ActiveType noHitProb = 1, S2(0.0), sigma2(0.0);
                VecA X2;
                for (size_t k = 0; k < times.size(); ++k) {
                    ActiveType S1(S2), sigma1(sigma2);
                    VecA X1(X2);
                    X2 = sim_->state(idx_, times[k]);
                    S2 = sim_->assetAdjuster(times[k], alias) + process_->asset(times[k], X2, alias); // use low-level calculation to avoid repeated state computation
                    // if we have a discrete hit then we don't need probabilties
                    if ((downOrUpOrBoth <= 0) && (S2 <= downBarrier)) return 0.0;
                    if ((downOrUpOrBoth >= 0) && (S2 >= upBarrier)) return 0.0;
                    sigma2 = process_->assetVolatility(times[k], X2, alias);  // this is the local log-volatility 
                    if (k == 0) continue; // we can't calculate prob's in the first iteration
                    ActiveType variance = (sigma1*sigma1 + sigma2*sigma2) * (times[k] - times[k - 1]) / 2.0;
                    ActiveType downHit = std::exp(-2 * std::log(downBarrier / S1) * std::log(downBarrier / S2) / variance);
                    ActiveType upHit = std::exp(-2 * std::log(upBarrier / S1) * std::log(upBarrier / S2) / variance);
                    // we assme disjoint events for up and down hits
                    ActiveType hit = 0.0;
                    if (downOrUpOrBoth <= 0) hit += downHit;
                    if (downOrUpOrBoth >= 0) hit += upHit;
                    hit = (hit > 1.0) ? (1.0) : (hit);
                    noHitProb *= (1.0 - hit);
                }
                return noHitProb;
            }

            // for commodity payoffs
            inline ActiveType futureAsset(DateType obsTime, DateType settlementTime, const std::string& alias) {
                return	process_->futureAsset(obsTime, settlementTime, sim_->state(idx_, obsTime), alias);
            }
        };


        MCSimulationT( const ext::shared_ptr<ProcessType>   process,
                       const VecD&                          simTimes,
                       const VecD&                          obsTimes,
                       size_t                               nPaths,
                       BigNatural                           seed = 1234,
                       bool                                 richardsonExtrapolation = true,
                       bool                                 timeInterpolation = false,
                       bool                                 storeBrownians = false )
                    : process_(process), seed_(seed), richardsonExtrapolation_(richardsonExtrapolation),
                          timeInterpolation_(timeInterpolation), storeBrownians_(storeBrownians) {
            // check inputs
            QL_REQUIRE(simTimes.size()>0,"TemplateMCSimulation: non-empty simulation times required");
            for (size_t k=1; k<simTimes.size(); ++k) QL_REQUIRE(simTimes[k-1]<simTimes[k],"TemplateMCSimulation: simulation times in ascending order required");
            QL_REQUIRE(obsTimes.size()>0,"TemplateMCSimulation: non-empty observation times required");
            for (size_t k=1; k<obsTimes.size(); ++k) QL_REQUIRE(obsTimes[k-1]<obsTimes[k],"TemplateMCSimulation: observation times in ascending order required");
            // use only obsTimes>0
            obsTimes_.push_back(0.0);
            for (size_t k=0; k<obsTimes.size(); ++k) if (obsTimes[k]>0) obsTimes_.push_back(obsTimes[k]);
            // merge simTimes and obsTimes_ into simTimes_
            size_t sim_idx=0;
            while ((sim_idx<simTimes.size()) && (simTimes[sim_idx]<=0)) ++sim_idx;
            size_t obs_idx=0;
            while ((obs_idx<obsTimes_.size()) && (obsTimes_[obs_idx]<=0)) ++obs_idx;
            // we alway start at t=0
            simTimes_.push_back(0.0);
            // we use some tolerance for checking equal dates to avoid rounding issues
            // maybe even a larger tolerance than one day could make sense
            DateType OneDay = 1.0/365.25;
            while ((sim_idx<simTimes.size()) || (obs_idx<obsTimes_.size())) {
                if  ((sim_idx<simTimes.size())&&(obs_idx<obsTimes_.size())&&
                    (abs(simTimes[sim_idx]-obsTimes_[obs_idx])<OneDay)) {
                    simTimes_.push_back(obsTimes_[obs_idx]);
                    ++sim_idx;
                    ++obs_idx;
                    continue;
                }
                if ((obs_idx>=obsTimes_.size())||((sim_idx<simTimes.size())&&(simTimes[sim_idx]<obsTimes_[obs_idx]))) {
                    simTimes_.push_back(simTimes[sim_idx]);
                    ++sim_idx;
                    continue;
                } else {
                    simTimes_.push_back(obsTimes_[obs_idx]);
                    ++obs_idx;
                    continue;
                }
            }
            // allocate memory (we don't want surprises during time-consuming simulation)
            reallocateMemory(nPaths);
            // ready to simulate...
            // but first set up adjusters
            applyNumeraireAdjuster_ = false;  // defauls
            applyZcbAdjuster_       = false;  // default
            applyAssetAdjuster_     = false;  // default
            
        }

        inline void simulate() {  // the procedure Initialise/PreEvBrownians/Simulate needs to be re-factorised
            addObsTimes_.clear();  // reset interpolated states
            Y_.clear();
            initialiseRSG();
            if (storeBrownians_)  // That's not really ok. If storeBrownians we want to avoid repeated calculation...
                for (size_t k=0; k<dW_.size(); ++k) dW_[k] = getNextBrownianIncrements();
            for (size_t k=0; k<X_.size(); ++k) simulatePath(k);
        }

        // the following two routines are for sliced simulation

        inline void prepareForSlicedSimulation() {  // this routine checks constraints and prepares for below simulate(idx) calls
            //QL_REQUIRE(storeBrownians_ == true, "TemplateMCSimulation: storeBrownians required");
            //QL_REQUIRE(richardsonExtrapolation_==false, "TemplateMCSimulation: Richardson extrapolation not supported");
            QL_REQUIRE(simTimes_.size() == obsTimes_.size(), "TemplateMCSimulation: simTimes_ == obsTimes required");
            for (size_t k=0; k<simTimes_.size(); ++k)
                QL_REQUIRE(simTimes_[k] == obsTimes_[k], "TemplateMCSimulation: simTimes_ == obsTimes required");
            if (storeBrownians_) {
                initialiseRSG();
                for (size_t k = 0; k < dW_.size(); ++k) dW_[k] = getNextBrownianIncrements();
            }
            else { // in this case we need to set up the random sequence generator differently
                rsg_ = ext::shared_ptr<PseudoRandom::rsg_type>(new PseudoRandom::rsg_type(PseudoRandom::make_sequence_generator(process_->factors(), seed_)));
                dW_[0].resize(2); // for Richardson extrapolation we need two slices ob BM's
            }
            // better set a flag that indicates that all is ready for simulation
        }

        inline void simulate(const size_t idx, const bool clearInterpolation = false) {
            // we might want to clear interpolated states, e.g. in calibration runs
            if (clearInterpolation) {
                addObsTimes_.clear();
                Y_.clear();
            }
            QL_REQUIRE(idx < obsTimes_.size(), "idx<obsTimes_.size() required");
            if (idx == 0) {
                for (size_t path = 0; path < X_.size(); ++path) X_[path][0] = process_->initialValues();
                return;
            }
            Real dt = obsTimes_[idx] - obsTimes_[idx - 1];
            
            if (richardsonExtrapolation_) {
                VecA X0T = process_->initialValues();
                VecA X1T = VecA((X0T.size()));
                VecA X12T = VecA((X0T.size()));
                VecD *dW1, *dW2;  // we need these pointers to handle storeBrownians==false
                VecD dWT = VecD((process_->factors()));
                for (size_t path = 0; path < X_.size(); ++path) {
                    if (storeBrownians_) {
                        dW1 = &dW_[path][2 * (idx - 1)];
                        dW2 = &dW_[path][2 * (idx - 1) + 1];
                    }
                    else {
                        dW_[0][0] = rsg_->nextSequence().value;
                        dW_[0][1] = rsg_->nextSequence().value;
                        dW1 = &dW_[0][0];
                        dW2 = &dW_[0][1];
                    }
                    // full Euler step
                    for (size_t k = 0; k < dWT.size(); ++k) dWT[k] = ((*dW1)[k] + (*dW2)[k]) / sqrt(2.0);
                    process_->evolve(obsTimes_[idx - 1], X_[path][idx - 1], dt, dWT, X1T);
                    // two half size Euler steps
                    process_->evolve(obsTimes_[idx - 1], X_[path][idx - 1], dt / 2.0, *dW1, X12T);
                    process_->evolve(obsTimes_[idx - 1] + dt / 2.0, X12T, dt / 2.0, *dW2, X0T);
                    // extrapolation
                    for (size_t k = 0; k < X1T.size(); ++k) X1T[k] = 2 * X0T[k] - X1T[k];
                    // extrapolation may lead to ill-defined states
                    process_->truncate(obsTimes_[idx], X1T);
                    for (size_t k = 0; k < X1T.size(); ++k) X_[path][idx][k] = X1T[k];
                }
            }
            else { // only full Euler step
                for (size_t path = 0; path < X_.size(); ++path) {
                    VecD *dW;  // we need these pointer to handle storeBrownians==false
                    if (storeBrownians_) {
                        dW = &dW_[path][idx - 1];
                    }
                    else {
                        dW_[0][0] = rsg_->nextSequence().value;
                        dW = &dW_[0][0];
                    }
                    process_->evolve(obsTimes_[idx - 1], X_[path][idx - 1], dt, *dW, X_[path][idx]);
                }
            }
        }

        // inspectors

        inline const ext::shared_ptr<ProcessType> process() { return process_; }

        inline const VecD& simTimes() { return simTimes_; }
        inline const VecD& obsTimes() { return obsTimes_; }
        inline size_t      nPaths()   { return X_.size(); }

        inline const ext::shared_ptr<Path> path(const size_t idx) {
            return ext::shared_ptr<Path>(new Path(process_,this,idx));
        }

        inline const MatA& observedPath(const size_t idx) {
            QL_REQUIRE(idx<X_.size(),"TemplateMCSimulation: path out of bounds.");
            return X_[idx]; 
        }

        // find an actual state of the model for a given observation time
        // this method is used by a path object and the result is passed on to
        // the stochastic process for payoff evaluation
        inline const VecA& state(const size_t idx, const DateType t) {
            QL_REQUIRE(idx<X_.size(),"TemplateMCSimulation: path out of bounds.");
            size_t t_idx = TemplateAuxilliaries::idx(obsTimes_,t);
            if (t==obsTimes_[t_idx]) return X_[idx][t_idx];
            QL_REQUIRE(timeInterpolation_,"TemplateMCSimulation: time interpolation not allowed");
            // allow extrapolation
            if (t<obsTimes_[0])                   return X_[idx][0];
            if (t>obsTimes_[obsTimes_.size()-1])  return X_[idx][X_[idx].size()-1];
            // first check if we already have an additional state
            if (addObsTimes_.size() > 0) {
                size_t t_addIdx = TemplateAuxilliaries::idx(addObsTimes_, t);
                if (t == addObsTimes_[t_addIdx]) return (*Y_[t_addIdx])[idx];
            }
            // if we end up here we truely need to interpolate
            ext::shared_ptr<MatA> y(new MatA(nPaths(), VecA(process_->size())));
            // linear state interpolation (this is very crude) 
            DateType rho = (t - obsTimes_[t_idx - 1]) / (obsTimes_[t_idx] - obsTimes_[t_idx - 1]);
            for (size_t p = 0; p < nPaths(); ++p)
                for (size_t k = 0; k < process_->size(); ++k)
                    (*y)[p][k] = (1.0 - rho)*X_[p][t_idx - 1][k] + rho*X_[p][t_idx][k];
            // add new state to list to re-use it for later calls
            addObsTimes_.push_back(t);
            Y_.push_back(y);
            // we need to make sure additional times and states are sorted for efficient searching
            std::vector<size_t> sorted(addObsTimes_.size());
            for (size_t k = 0; k < sorted.size(); ++k) sorted[k] = k;
            std::sort(sorted.begin(), sorted.end(), LessByVector(addObsTimes_));
            VecD tmpTimes(addObsTimes_);
            for (size_t k = 0; k < addObsTimes_.size(); ++k) addObsTimes_[k] = tmpTimes[sorted[k]];
            std::vector< ext::shared_ptr< MatA > > tmpY(Y_);
            for (size_t k = 0; k < Y_.size(); ++k) Y_[k] = tmpY[sorted[k]];
            // all done. now we only need to return the current interpolated state...
            return (*y)[idx];
        }

        inline const MatA& brownian(const size_t idx) {
            QL_REQUIRE(idx<dW_.size(),"TemplateMCSimulation: path out of bounds.");
            return dW_[idx]; 
        }

        // numeraire adjuster

        inline void calculateNumeraireAdjuster(const VecD&  numeraireObservTimes) {
            numeraireObservTimes_ = numeraireObservTimes;
            // check time grids
            QL_REQUIRE(numeraireObservTimes_.size()>1, "TemplateMCSimulation: at least two numeraireObservTimes_ required");
            QL_REQUIRE(numeraireObservTimes_[0] > 0, "TemplateMCSimulation: numeraireObservTimes_[0]>0 required");
            for (size_t k = 1; k<numeraireObservTimes_.size(); ++k) QL_REQUIRE(numeraireObservTimes_[k - 1]<numeraireObservTimes_[k], "TemplateMCSimulation: numeraireObservTimes_ in ascending order required");
            // initialise adjuster
            numeraireAdjuster_ = VecA(numeraireObservTimes_.size(), 0.0);
            for (size_t k = 0; k<nPaths(); ++k) {
                for (size_t i = 0; i<numeraireObservTimes_.size(); ++i) {
                    VecA       s = state(k, numeraireObservTimes_[i]);
                    numeraireAdjuster_[i] += 1.0 / process_->numeraire(numeraireObservTimes_[i], s);
                }
            }
            for (size_t i = 0; i < numeraireObservTimes_.size(); ++i) {
                numeraireAdjuster_[i] /= nPaths();  //  E{ 1/N(t) }
                numeraireAdjuster_[i] /= process_->zeroBond(0.0, numeraireObservTimes_[i], process_->initialValues()); // adj =  E{ 1/N(t) } / P(0,t)
                numeraireAdjuster_[i] = log(numeraireAdjuster_[i]) / numeraireObservTimes_[i];
            }

            // calculate adjuster..., tbd
            applyNumeraireAdjuster_ = true;
        }

        inline const VecA& numeraireAdjuster() { return numeraireAdjuster_; }

        inline ActiveType numeraireAdjuster(const DateType t) {
            if (!applyNumeraireAdjuster_) return 1.0; // implement adjuster interpolation
            // linear interpolation
            size_t obsIdx = TemplateAuxilliaries::idx(numeraireObservTimes_, t);
            if (obsIdx<1) obsIdx = 1;
            DateType rhoObs = (t - numeraireObservTimes_[obsIdx - 1]) / (numeraireObservTimes_[obsIdx] - numeraireObservTimes_[obsIdx - 1]);
            // flat extrapolation
            if (rhoObs<0) rhoObs = 0;
            if (rhoObs>1) rhoObs = 1;
            //
            ActiveType z = numeraireAdjuster_[obsIdx - 1] * (1 - rhoObs) + numeraireAdjuster_[obsIdx] * rhoObs;
            return exp(z*t);
        }

        // zero coupon bond adjuster

        inline void calculateZCBAdjuster( const VecD&  zcbObservTimes,
                                          const VecD&  zcbOffsetTimes ) {
            zcbObservTimes_ = zcbObservTimes;
            zcbOffsetTimes_ = zcbOffsetTimes;
            // check time grids
            QL_REQUIRE(zcbObservTimes_.size()>1,"TemplateMCSimulation: at least two zcbObservTimes_ required");
            QL_REQUIRE(zcbObservTimes_[0]>0,"TemplateMCSimulation: zcbObservTimes_[0]>0 required");
            for (size_t k=1; k<zcbObservTimes_.size(); ++k) QL_REQUIRE(zcbObservTimes_[k-1]<zcbObservTimes_[k],"TemplateMCSimulation: zcbObservTimes_ in ascending order required");
            QL_REQUIRE(zcbOffsetTimes_.size()>1,"TemplateMCSimulation: at least two zcbOffsetTimes_ required");
            QL_REQUIRE(zcbOffsetTimes_[0]>0,"TemplateMCSimulation: zcbOffsetTimes_[0]>0 required");
            for (size_t k=1; k<zcbOffsetTimes_.size(); ++k) QL_REQUIRE(zcbOffsetTimes_[k-1]<zcbOffsetTimes_[k],"TemplateMCSimulation: zcbOffsetTimes_ in ascending order required");
            // initialise zero adjuster matrix
            zcbAdjuster_ = MatA(zcbObservTimes_.size(),VecA(zcbOffsetTimes_.size(),0.0));
            MatA zcb(zcbObservTimes_.size(),VecA(zcbOffsetTimes_.size(),0.0));
            for (size_t k=0; k<nPaths(); ++k) {
                for (size_t i=0; i<zcbObservTimes_.size(); ++i) {
                    VecA       s   = state(k,zcbObservTimes_[i]);
                    ActiveType num = this->numeraireAdjuster(zcbObservTimes_[i]) * process_->numeraire(zcbObservTimes_[i],s);
                    for (size_t j=0; j<zcbOffsetTimes_.size(); ++j) {
                        zcb[i][j] += process_->zeroBond(zcbObservTimes_[i], zcbObservTimes_[i]+zcbOffsetTimes_[j], s ) / num;
                    }
                }
            }
            for (size_t i=0; i<zcbObservTimes_.size(); ++i) {
                for (size_t j=0; j<zcbOffsetTimes_.size(); ++j) {
                    ActiveType adjDF = process_->zeroBond(0.0, zcbObservTimes_[i]+zcbOffsetTimes_[j], process_->initialValues() );
                    // we ommited division by nuneraire(0) = 1 here
                    zcb[i][j] /= nPaths();
                    adjDF /= zcb[i][j];
                    zcbAdjuster_[i][j] = - log(adjDF) / (zcbOffsetTimes_[j]);
                }
            }
            applyZcbAdjuster_ = true;  // ready to use zcb adjuster
        }

        inline const MatA& zcbAdjuster() { return zcbAdjuster_; }

        inline ActiveType zcbAdjuster(const DateType t, const DateType T) {
            if (!applyZcbAdjuster_) return 1.0; // implement adjuster interpolation
            DateType dt = T - t;
            //if (dt<0) return 1.0;
            // bilinear interpolation
            size_t obsIdx = TemplateAuxilliaries::idx(zcbObservTimes_,t);
            size_t offIdx = TemplateAuxilliaries::idx(zcbOffsetTimes_,dt);
            if (obsIdx<1) obsIdx=1;
            if (offIdx<1) offIdx=1;
            DateType rhoObs = ( t-zcbObservTimes_[obsIdx-1])/(zcbObservTimes_[obsIdx]-zcbObservTimes_[obsIdx-1]);
            DateType rhoOff = (dt-zcbOffsetTimes_[offIdx-1])/(zcbOffsetTimes_[offIdx]-zcbOffsetTimes_[offIdx-1]);
            // flat extrapolation
            if (rhoObs<0) rhoObs = 0;
            if (rhoObs>1) rhoObs = 1;
            if (rhoOff<0) rhoOff = 0;
            if (rhoOff>1) rhoOff = 1;
            //
            ActiveType z = zcbAdjuster_[obsIdx-1][offIdx-1] * (1.0-rhoObs) * (1.0-rhoOff) +
                           zcbAdjuster_[obsIdx  ][offIdx-1] * (rhoObs)     * (1.0-rhoOff) +
                           zcbAdjuster_[obsIdx-1][offIdx  ] * (1.0-rhoObs) * (rhoOff)     +
                           zcbAdjuster_[obsIdx  ][offIdx  ] * (rhoObs)     * (rhoOff)     ;
            return exp(-z*(T-t));
        }

        inline ActiveType zcbAdjuster(const DateType t, const DateType T, const std::string& alias) {
            return 1.0;  // tbd. This needs to be implemented
        }

        // asset adjuster

        inline void calculateAssetAdjuster( const VecD&  assetObservTimes, const std::vector<std::string>& aliases ) {
            assetObservTimes_ = assetObservTimes;
            // check time grids
            QL_REQUIRE(assetObservTimes_.size()>1,"TemplateMCSimulation: at least two assetObservTimes_ required");
            QL_REQUIRE(assetObservTimes_[0]>=0,"TemplateMCSimulation: assetObservTimes_>=0 required");
            for (size_t k=1; k<assetObservTimes_.size(); ++k) QL_REQUIRE(assetObservTimes_[k-1]<assetObservTimes_[k],"TemplateMCSimulation: assetObservTimes_ in ascending order required");
            QL_REQUIRE(aliases.size() > 0, "TemplateMCSimulation: aliases required");
            for (size_t k = 0; k < aliases.size(); ++k) assetIndex_[aliases[k]] = k;
            QL_REQUIRE(assetIndex_.size() == aliases.size(), "TemplateMCSimulation: duplicate aliases found");

            // initialise asset adjuster vector
            assetAdjuster_ = MatA(aliases.size(),VecA(assetObservTimes_.size(), 0.0));
            for (size_t k = 0; k < assetAdjuster_.size(); ++k) {
                VecA avAsset(assetObservTimes_.size(), 0.0);
                VecA avZero(assetObservTimes_.size(), 0.0);  // numeraire
                for (size_t i = 0; i<assetObservTimes_.size(); ++i) {
                    for (size_t j = 0; j<nPaths(); ++j) {
                        VecA       s = state(j, assetObservTimes_[i]);
                        avAsset[i] += process_->asset(assetObservTimes_[i], s, aliases[k]) / (this->numeraireAdjuster(assetObservTimes_[i])*process_->numeraire(assetObservTimes_[i], s));
                        avZero[i] += 1.0 / (this->numeraireAdjuster(assetObservTimes_[i])*process_->numeraire(assetObservTimes_[i], s));
                    }
                    avAsset[i] /= nPaths();
                    avZero[i] /= nPaths();
                    // discounted expected asset (in terminal meassure)
                    assetAdjuster_[k][i] = process_->zeroBond(0, assetObservTimes_[i], process_->initialValues()) *
                        process_->forwardAsset(0, assetObservTimes_[i], process_->initialValues(), aliases[k]);
                    assetAdjuster_[k][i] = (assetAdjuster_[k][i] - avAsset[i]) / avZero[i];
                }

            }
            applyAssetAdjuster_ = true;
        }

        inline const VecA& assetAdjuster(const std::string& alias) { return assetAdjuster_[assetIndex_.at(alias)]; }

        inline ActiveType assetAdjuster(const DateType t, const std::string& alias) {
            if (!applyAssetAdjuster_) return 0.0; 
            // linear interpolation
            size_t obsIdx = TemplateAuxilliaries::idx(assetObservTimes_,t);
            if (obsIdx<1) obsIdx=1;
            DateType rhoObs = ( t-assetObservTimes_[obsIdx-1])/(assetObservTimes_[obsIdx]-assetObservTimes_[obsIdx-1]);
            // flat extrapolation
            if (rhoObs<0) rhoObs = 0;
            if (rhoObs>1) rhoObs = 1;
            ActiveType adj = (1.0-rhoObs) * assetAdjuster_[assetIndex_.at(alias)][obsIdx-1] + rhoObs * assetAdjuster_[assetIndex_.at(alias)][obsIdx];
            return adj;
        }

    };

}

#endif  /* ifndef quantlib_templatemcsimulation_hpp */
