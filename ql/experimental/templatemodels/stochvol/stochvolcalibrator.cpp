/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2015 Sebastian Schlenkrich

*/

#include <ql/pricingengines/blackformula.hpp>

#include <ql/experimental/templatemodels/stochvol/stochvolcalibrator.hpp>

namespace QuantLib {

        // count number of effective states
        size_t StochVolModelCalibrator::xDim() const {
            size_t N=0;
            if (!lambdaIsFixed_ ) ++N;
            if (!bIsFixed_		) ++N;
            if (!LIsFixed_		) ++N;
            if (!thetaIsFixed_	) ++N;
            if (!mIsFixed_		) ++N;
            if (!etaIsFixed_	) ++N;
            if (!z0IsFixed_		) ++N;
            if (!rhoIsFixed_	) ++N;
            return N;
        }

        Real StochVolModelCalibrator::direct(const Real x, const Real a, const Real b) const {
            return (b-a)*(atan(x)/M_PI + 0.5) + a;
        }

        Real StochVolModelCalibrator::inverse(const Real y, const Real a, const Real b) const {
            return tan( ((y-a)/(b-a)-0.5) * M_PI );
        }

        Array StochVolModelCalibrator::direct(const Array& X) const {
            QL_REQUIRE(X.size()==xDim(),"StochVolModelCalibrator Error: dimension mismatch.");
            Array Y(X.size());
            size_t idx = 0;
            if (!lambdaIsFixed_ ) { Y[idx] = direct(X[idx], lambdaMin_, lambdaMax_);    ++idx; }   
            if (!bIsFixed_		) { Y[idx] = direct(X[idx], bMin_     , bMax_     );    ++idx; }   
            if (!LIsFixed_		) { Y[idx] = direct(X[idx], LMin_     , LMax_     );    ++idx; }   
            if (!thetaIsFixed_	) { Y[idx] = direct(X[idx], thetaMin_ , thetaMax_ );    ++idx; }   
            if (!mIsFixed_		) { Y[idx] = direct(X[idx], mMin_     , mMax_     );    ++idx; }   
            if (!etaIsFixed_	) { Y[idx] = direct(X[idx], etaMin_   , etaMax_   );    ++idx; }   
            if (!z0IsFixed_		) { Y[idx] = direct(X[idx], z0Min_    , z0Max_    );    ++idx; }   
            if (!rhoIsFixed_	) { Y[idx] = direct(X[idx], rhoMin_   , rhoMax_   );    ++idx; }   
            return Y;
        }

        Array StochVolModelCalibrator::inverse(const Array& Y) const {
            QL_REQUIRE(Y.size()==xDim(),"StochVolModelCalibrator Error: dimension mismatch.");
            Array X(Y.size());
            size_t idx = 0;
            if (!lambdaIsFixed_ ) { X[idx] = inverse(Y[idx], lambdaMin_, lambdaMax_);    ++idx; }
            if (!bIsFixed_		) { X[idx] = inverse(Y[idx], bMin_     , bMax_     );    ++idx; }
            if (!LIsFixed_		) { X[idx] = inverse(Y[idx], LMin_     , LMax_     );    ++idx; }
            if (!thetaIsFixed_	) { X[idx] = inverse(Y[idx], thetaMin_ , thetaMax_ );    ++idx; }
            if (!mIsFixed_		) { X[idx] = inverse(Y[idx], mMin_     , mMax_     );    ++idx; }
            if (!etaIsFixed_	) { X[idx] = inverse(Y[idx], etaMin_   , etaMax_   );    ++idx; }
            if (!z0IsFixed_		) { X[idx] = inverse(Y[idx], z0Min_    , z0Max_    );    ++idx; }
            if (!rhoIsFixed_	) { X[idx] = inverse(Y[idx], rhoMin_   , rhoMax_   );    ++idx; }
            return X;
        }

        // initialize state X with model parameters and apply inverse transformation
        Array StochVolModelCalibrator::initialise() {
            Array X(xDim());
            size_t idx = 0;
            if (!lambdaIsFixed_ ) { X[idx] = lambda_;   ++idx; }
            if (!bIsFixed_		) { X[idx] = b_;		++idx; } 
            if (!LIsFixed_		) { X[idx] = L_;		++idx; } 
            if (!thetaIsFixed_	) { X[idx] = theta_;	++idx; } 
            if (!mIsFixed_		) { X[idx] = m_;		++idx; } 
            if (!etaIsFixed_	) { X[idx] = eta_;	    ++idx; }
            if (!z0IsFixed_		) { X[idx] = z0_;	    ++idx; }
            if (!rhoIsFixed_	) { X[idx] = rho_;	    ++idx; }
            return inverse(X);
        }

        // create model for given state vector
        ext::shared_ptr<RealStochVolModel> StochVolModelCalibrator::model(const Array& X) const {
            // model parameters
            Real  lambda;
            Real  b;
            Real  L;
            Real  theta;
            Real  m;
            Real  eta;
            Real  z0;
            Real  rho;
            // get parameters
            QL_REQUIRE(X.size()==xDim(),"StochVolModelCalibrator Error: dimension mismatch.");
            size_t idx = 0;
            if (!lambdaIsFixed_ ) { lambda = X[idx];  ++idx; } else { lambda = lambda_; }
            if (!bIsFixed_		) { b      = X[idx];  ++idx; } else { b      = b_;     	}
            if (!LIsFixed_		) { L      = X[idx];  ++idx; } else { L      = L_;     	}
            if (!thetaIsFixed_	) { theta  = X[idx];  ++idx; } else { theta  = theta_; 	}
            if (!mIsFixed_		) { m	   = X[idx];  ++idx; } else { m	     = m_;	   	}
            if (!etaIsFixed_	) { eta    = X[idx];  ++idx; } else { eta    = eta_;    }
            if (!z0IsFixed_		) { z0     = X[idx];  ++idx; } else { z0     = z0_;     }
            if (!rhoIsFixed_	) { rho    = X[idx];  ++idx; } else { rho    = rho_;    }
            return ext::shared_ptr<RealStochVolModel>(new RealStochVolModel(lambda,b,L,theta,m,eta,z0,rho));
        }

        Array StochVolModelCalibrator::values(const Array& x) const {
            ext::shared_ptr<RealStochVolModel> m = model(direct(x));
            std::vector<Real>          fwPrices(strikes_.size());
            std::vector<Option::Type>  callOrPut(strikes_.size());
            std::vector<Real>          normalVols(strikes_.size());
            Array                      volVariance(strikes_.size());
            for (size_t k=0; k<strikes_.size(); ++k) {
                callOrPut[k]   = (strikes_[k] < forward_) ? (Option::Put) : (Option::Call);
                fwPrices[k]    = m->vanillaOption(forward_,strikes_[k],exercTime_,callOrPut[k],glAbsAcc_,glMaxEval_);
                normalVols[k]  = bachelierBlackFormulaImpliedVol(callOrPut[k],strikes_[k],forward_,exercTime_,fwPrices[k]);
                volVariance[k] = normalVols[k] - vols_[k];
            }
            return volVariance;
        }

        Real StochVolModelCalibrator::value(const Array& x) const {
            Array y = values(x);
            Real res=0.0;
            for (size_t k=0; k<y.size(); ++k) res += y[k]*y[k];
            return res/2.0;
        }

        // update model parameters with state X
        void StochVolModelCalibrator::update(const Array& X) {
            QL_REQUIRE(X.size()==xDim(),"StochVolModelCalibrator Error: dimension mismatch.");
            size_t idx = 0;
            if (!lambdaIsFixed_ ) { lambda_ = X[idx];  ++idx; }
            if (!bIsFixed_		) { b_      = X[idx];  ++idx; } 
            if (!LIsFixed_		) { L_      = X[idx];  ++idx; } 
            if (!thetaIsFixed_	) { theta_  = X[idx];  ++idx; } 
            if (!mIsFixed_		) { m_	    = X[idx];  ++idx; } 
            if (!etaIsFixed_	) { eta_    = X[idx];  ++idx; }
            if (!z0IsFixed_		) { z0_     = X[idx];  ++idx; }
            if (!rhoIsFixed_	) { rho_    = X[idx];  ++idx; }
        }

        // construct object and calibrate
        StochVolModelCalibrator::StochVolModelCalibrator( // initial model parameters
                                 const Real  lambda,
                                 const Real  b,
                                 const Real  L,
                                 const Real  theta,
                                 const Real  m,
                                 const Real  eta,
                                 const Real  z0,
                                 const Real  rho,
                                 // calibration flags
                                 const bool  lambdaIsFixed,
                                 const bool  bIsFixed,
                                 const bool  LIsFixed,
                                 const bool  thetaIsFixed,
                                 const bool  mIsFixed,
                                 const bool  etaIsFixed,
                                 const bool  z0IsFixed,
                                 const bool  rhoIsFixed,
                                 // calibration targets
                                 const Real  exercTime,
                                 const Real  forward,
                                 const std::vector<Real>& strikes,
                                 const std::vector<Real>& vols,
                                 const std::vector<Real>& optimizationParams)  // { [min], [max], epsfcn, ftol, xtol, gtol, maxfev, glAbsAcc, glMaxEval }							 								 )
        : lambda_(lambda), b_(b), L_(L), theta_(theta), m_(m), eta_(eta), z0_(z0), rho_(rho),
          lambdaIsFixed_(lambdaIsFixed), bIsFixed_(bIsFixed), LIsFixed_(LIsFixed), thetaIsFixed_(thetaIsFixed), mIsFixed_(mIsFixed), etaIsFixed_(etaIsFixed), z0IsFixed_(z0IsFixed), rhoIsFixed_(rhoIsFixed),
          exercTime_(exercTime), forward_(forward), strikes_(strikes), vols_(vols) {
            // default optimisation parameters
            // lower boundaries
            lambdaMin_ =  0.0;
            bMin_      =  0.0;
            LMin_      =  0.0;
            thetaMin_  =  0.0;
            mMin_      =  0.0;
            etaMin_    =  0.0;
            z0Min_     =  0.0;
            rhoMin_    = -1.0;
            // upper boundaries
            lambdaMax_ =  3.0;
            bMax_      =  1.0;
            LMax_      =  1.0;
            thetaMax_  =  1.0;
            mMax_      =  3.0;
            etaMax_    =  3.0;
            z0Max_     =  3.0;
            rhoMax_    =  1.0;
            // Levenberg Marquard
            epsfcn_     =  1.0e-10;  // relative accuracy function evaluation, compare w/ glAbsAcc
            ftol_       =  1.0e-8; 
            xtol_       =  1.0e-8;
            gtol_       =  1.0e-8;
            maxfev_     =  10000;
            // Vanilla option Gauss Lobatto
            glAbsAcc_   =  1.0e-12;
            glMaxEval_  =  10000;
            // user-defined optimisation parameters
            // lower boundaries
            if (optimizationParams.size()>=8) {
                lambdaMin_ = optimizationParams[0];
                bMin_      = optimizationParams[1];
                LMin_      = optimizationParams[2];
                thetaMin_  = optimizationParams[3];
                mMin_      = optimizationParams[4];
                etaMin_    = optimizationParams[5];
                z0Min_     = optimizationParams[6];
                rhoMin_    = optimizationParams[7];
            }
            // upper boundaries
            if (optimizationParams.size()>=16) {
                lambdaMax_ = optimizationParams[8];
                bMax_      = optimizationParams[9];
                LMax_      = optimizationParams[10];
                thetaMax_  = optimizationParams[11];
                mMax_      = optimizationParams[12];
                etaMax_    = optimizationParams[13];
                z0Max_     = optimizationParams[14];
                rhoMax_    = optimizationParams[15];
            }
            // Levenberg Marquard
            if (optimizationParams.size()>=21) {
                epsfcn_    = optimizationParams[16];
                ftol_      = optimizationParams[17];
                xtol_      = optimizationParams[18];
                gtol_      = optimizationParams[19];
                maxfev_    = static_cast<Size>(optimizationParams[20]);
            }
            // Vanilla option Gauss Lobatto
            if (optimizationParams.size()>=23) {
                glAbsAcc_  =  optimizationParams[21];
                glMaxEval_ =  static_cast<Size>(optimizationParams[22]);
            }
            // set up
            NoConstraint constraint;
            Array x = initialise();
            {
                Problem problem(*this, constraint, x);
                LevenbergMarquardt optimizationMethod(epsfcn_, xtol_, gtol_);
                //EndCriteria endCriteria(60000, 100, 1e-4, 1e-4, 1e-4);
                EndCriteria endCriteria(maxfev_, 100 /* unused */, 0 /* unused */, ftol_, 0 /* unused */);
                // calibrate
                optimizationMethod.minimize(problem,endCriteria);
                update(direct(problem.currentValue()));
            }
        }


}

