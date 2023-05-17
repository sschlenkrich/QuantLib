/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2015 Sebastian Schlenkrich

*/



#ifndef quantlib_templatevanillalocalvolmodel_hpp
#define quantlib_templatevanillalocalvolmodel_hpp

#include <vector>
#include <string>

#include <stdio.h>

#include <ql/shared_ptr.hpp>
#include <ql/errors.hpp>

#include <ql/experimental/templatemodels/auxilliaries/auxilliariesT.hpp>
#include <ql/experimental/templatemodels/stochasticprocessT.hpp>


namespace QuantLib {

    template <class DateType, class PassiveType, class ActiveType>    
    class VanillaLocalVolModelT : public StochasticProcessT<DateType,PassiveType,ActiveType> {
        // from base class
        using typename StochasticProcessT<DateType,PassiveType,ActiveType>::VecD;
        using typename StochasticProcessT<DateType,PassiveType,ActiveType>::VecP;
        using typename StochasticProcessT<DateType,PassiveType,ActiveType>::VecA;
        using typename StochasticProcessT<DateType,PassiveType,ActiveType>::MatA;


    private:
        // input parameters
        DateType    T_;                    // time to expiry (in years)
        PassiveType S0_;                   // forward
        PassiveType sigmaATM_;             // ATM normal volatility as basis for straddle calculation
        std::vector<PassiveType> Sp_;      // S_i with S_i > S0
        std::vector<PassiveType> Sm_;      // S_-i with S_-i < S0
        std::vector<PassiveType> Mp_;      // slope on intervall [S_i-1, S_i)
        std::vector<PassiveType> Mm_;      // slope on intervall (S_-i, S_-[i-1]]
        // calculated parameters
        ActiveType  straddleATM_;          // ATM straddle price
        ActiveType  sigma0_;               // local vol at S0, i.e. sigma(S0)
        std::vector<ActiveType>  sigmaP_;  // sigma(Sp_[i])
        std::vector<ActiveType>  sigmaM_;  // sigma(Sm_[i])
        std::vector<PassiveType> Xp_;      // X_i with X_i > 0
        std::vector<PassiveType> Xm_;      // X_-i with X_-i < 0
        // adjusters
        PassiveType              mu_;      // in-the-model adjuster for forward
        PassiveType              alpha_;   // out-of-the-model adjuster for straddle
        PassiveType              nu_;      // out-of-the-model adjuster for forward
        // numerical accuracy parameters (maybe expose to user...)
        PassiveType              extrapolationStdevs_;  // number of stdevs used as lower and upper cutoff, default 10
        size_t                   maxCalibrationIters_;  // number of iterations for forward/sigma0 calibration
        size_t                   onlyForwardCalibrationIters_;  // we may have some initial iterations only calibrating forward, this is intended to stabilise calibration
        PassiveType              sigma0Tol_;            // tolerance for sigma convergence
        PassiveType              S0Tol_;                // tolerance for forward convergence
        bool                     adjustATM_;            // apply post-calibration ATM adjuster
        bool                     useInitialMu_;
        PassiveType              initialMu_;														
        // we may want some debug information for the calibration process
        bool                     enableLogging_;
        std::vector<std::string> logging_;
 
        // we have two constructors and want to make sure the setup is consistent
        inline void initializeDeepInTheModelParameters() {
            straddleATM_ = sigmaATM_ * sqrt(T_) * M_1_SQRTPI * M_SQRT_2 * 2.0;
            if (useInitialMu_) mu_ = initialMu_;
            else mu_ = -(Mm_[0] + Mp_[0]) / 4.0 * T_; // this should be exact for shifted log-normal models
            alpha_ = 1.0;
            nu_ = 0.0;
            extrapolationStdevs_ = 10.0;
            sigma0Tol_ = 1.0e-12;
            S0Tol_ = 1.0e-12;
        }

        // determine the lower and upper bounds for integration
        inline PassiveType lowerBoundX() { return  -extrapolationStdevs_ * sqrt(T_) + mu_; }
        inline PassiveType upperBoundX() { return   extrapolationStdevs_ * sqrt(T_) + mu_; }		

        inline PassiveType localVol(const bool isRightWing, const size_t k, const PassiveType S) {
            // this is an unsafe method specifying the vol function sigma(S) on the individual segments
            if (isRightWing) {
                QL_REQUIRE(k < sigmaP_.size(), "k < sigmaP_.size() required.");
                PassiveType sigma0 = (k > 0) ? sigmaP_[k - 1] : sigma0_;
                PassiveType S0 = (k > 0) ? Sp_[k - 1] : S0_;
                return sigma0 + Mp_[k] * (S - S0);
            } else {
                QL_REQUIRE(k < sigmaM_.size(), "k < sigmaM_.size() required.");
                PassiveType sigma0 = (k > 0) ? sigmaM_[k - 1] : sigma0_;
                PassiveType S0 = (k > 0) ? Sm_[k - 1] : S0_;
                return sigma0 + Mm_[k] * (S - S0);
            }
            return 0.0;  // this should never be reached
        }

        inline PassiveType underlyingS(const bool isRightWing, const size_t k, const PassiveType x) {
            // this is an unsafe method specifying the underlying level S(x) on the individual segments
            if (isRightWing) {
                QL_REQUIRE(k < sigmaP_.size(), "k < sigmaP_.size() required.");
                PassiveType sigma0 = (k > 0) ? sigmaP_[k - 1] : sigma0_;
                PassiveType x0 = (k > 0) ? Xp_[k - 1] : 0.0;
                PassiveType deltaS = (Mp_[k] == 0.0) ? (sigma0*(x - x0)) : (sigma0 / Mp_[k] * (exp(Mp_[k] * (x - x0)) - 1.0));
                PassiveType S0 = (k > 0) ? Sp_[k - 1] : S0_;
                return S0 + deltaS;
            }
            else {
                QL_REQUIRE(k < sigmaM_.size(), "k < sigmaM_.size() required.");
                PassiveType sigma0 = (k > 0) ? sigmaM_[k - 1] : sigma0_;
                PassiveType x0 = (k > 0) ? Xm_[k - 1] : 0.0;
                PassiveType deltaS = (Mm_[k] == 0.0) ? (sigma0*(x - x0)) : (sigma0 / Mm_[k] * (exp(Mm_[k] * (x - x0)) - 1.0));
                PassiveType S0 = (k > 0) ? Sm_[k - 1] : S0_;
                return S0 + deltaS;
            }
            return 0.0;  // this should never be reached
        }

        inline PassiveType underlyingX(const bool isRightWing, const size_t k, const PassiveType S) {
            // this is an unsafe method specifying the underlying level x(S) on the individual segments
            if (isRightWing) {
                QL_REQUIRE(k < sigmaP_.size(), "k < sigmaP_.size() required.");
                PassiveType sigma0 = (k > 0) ? sigmaP_[k - 1] : sigma0_;
                QL_REQUIRE(sigma0 > 0.0, "sigma0 > 0.0 required");
                PassiveType x0 = (k > 0) ? Xp_[k - 1] : 0.0;
                PassiveType S0 = (k > 0) ? Sp_[k - 1] : S0_;
                PassiveType deltaX = (Mp_[k] == 0.0) ? ((S - S0) / sigma0) : (log(1.0 + Mp_[k] / sigma0*(S - S0)) / Mp_[k]);
                return x0 + deltaX;
            }
            else {
                QL_REQUIRE(k < sigmaM_.size(), "k < sigmaM_.size() required.");
                PassiveType sigma0 = (k > 0) ? sigmaM_[k - 1] : sigma0_;
                QL_REQUIRE(sigma0 > 0.0, "sigma0 > 0.0 required");
                PassiveType x0 = (k > 0) ? Xm_[k - 1] : 0.0;
                PassiveType S0 = (k > 0) ? Sm_[k - 1] : S0_;
                PassiveType deltaX = (Mm_[k] == 0.0) ? ((S - S0) / sigma0) : (log(1.0 + Mm_[k] / sigma0*(S - S0)) / Mm_[k]);
                return x0 + deltaX;
            }
            return 0.0;  // this should never be reached
        }

        inline PassiveType primitiveF(const bool isRightWing, const size_t k, const PassiveType x) {
            // this is an unsafe method specifying the primitive function F(x) = \int [alpha S(x) + nu] p(x) dx
            // on the individual segments
            PassiveType sigma0, x0, S0, m0;
            if (isRightWing) {
                QL_REQUIRE(k < sigmaP_.size(), "k < sigmaP_.size() required.");
                sigma0 = (k > 0) ? sigmaP_[k - 1] : sigma0_;
                x0 = (k > 0) ? Xp_[k - 1] : 0.0;
                S0 = (k > 0) ? Sp_[k - 1] : S0_;
                m0 = Mp_[k];
            }
            else {
                QL_REQUIRE(k < sigmaM_.size(), "k < sigmaM_.size() required.");
                sigma0 = (k > 0) ? sigmaM_[k - 1] : sigma0_;
                x0 = (k > 0) ? Xm_[k - 1] : 0.0;
                S0 = (k > 0) ? Sm_[k - 1] : S0_;
                m0 = Mm_[k];
            }
            PassiveType	y0 = (x0 - mu_) / sqrt(T_);
            PassiveType y1 = (x - mu_) / sqrt(T_);
            PassiveType h = m0 * sqrt(T_);
            PassiveType Ny = TemplateAuxilliaries::Phi(y1);
            PassiveType term1, term2;
            if (m0 == 0.0) {
                term1 = (S0 + nu_ / alpha_ - sigma0 * sqrt(T_) * y0) * Ny;
                term2 = sigma0 * sqrt(T_) * TemplateAuxilliaries::phi(y1);  // use dN/dx = dN/dy / sqrt(T)
            }
            else {				
                PassiveType NyMinush = TemplateAuxilliaries::Phi(y1 - h);
                term1 = exp(h*h / 2.0 - h*y0)*sigma0 / m0 * NyMinush;
                term2 = (sigma0 / m0 - (S0 + nu_ / alpha_)) * Ny;
            }
            return alpha_ * (term1 - term2);
        }

        inline PassiveType primitiveFSquare(const bool isRightWing, const size_t k, const PassiveType x) {
            // this is an unsafe method specifying the primitive function F(x) = \int [alpha S(x) + nu]^2 p(x) dx
            // on the individual segments
            PassiveType sigma0, x0, S0, m0;
            if (isRightWing) {
                QL_REQUIRE(k < sigmaP_.size(), "k < sigmaP_.size() required.");
                sigma0 = (k > 0) ? sigmaP_[k - 1] : sigma0_;
                x0 = (k > 0) ? Xp_[k - 1] : 0.0;
                S0 = (k > 0) ? Sp_[k - 1] : S0_;
                m0 = Mp_[k];
            }
            else {
                QL_REQUIRE(k < sigmaM_.size(), "k < sigmaM_.size() required.");
                sigma0 = (k > 0) ? sigmaM_[k - 1] : sigma0_;
                x0 = (k > 0) ? Xm_[k - 1] : 0.0;
                S0 = (k > 0) ? Sm_[k - 1] : S0_;
                m0 = Mm_[k];
            }
            PassiveType	y0 = (x0 - mu_) / sqrt(T_);
            PassiveType y1 = (x - mu_) / sqrt(T_);
            PassiveType h = m0 * sqrt(T_);
            PassiveType Ny = TemplateAuxilliaries::Phi(y1);
            PassiveType sum = 0;
            if (m0 == 0.0) {
                PassiveType K3 = S0 + nu_ / alpha_ - sigma0 * sqrt(T_) * y0;
                PassiveType term1 = (K3 * K3 + sigma0 * sigma0 * T_) * Ny;
                PassiveType term2 = 2.0 * sigma0 * sqrt(T_) * K3 + sigma0 * sigma0 * T_ * y1;
                term2 *= TemplateAuxilliaries::phi(y1);  // use dN/dx = dN/dy / sqrt(T)
                sum = term1 - term2;
            }
            else {
                PassiveType NyMinush  = TemplateAuxilliaries::Phi(y1 - h);
                PassiveType NyMinus2h = TemplateAuxilliaries::Phi(y1 - 2.0*h);
                PassiveType K1 = sigma0 / m0 * exp(h*(h - y0));
                PassiveType K2 = S0 + nu_ / alpha_ - sigma0 / m0;
                PassiveType term1 = K2 * K2 * Ny;
                PassiveType term2 = 2.0 * K1 * K2 * exp(-h*h / 2.0) * NyMinush;
                PassiveType term3 = K1 * K1 * NyMinus2h;
                sum = term1 + term2 + term3;
            }
            return alpha_ * alpha_ * sum;
        }

        inline void calculateSGrid() {
            // this is an unsafe method to calculate the S-grid for a given x-grid
            // it is intended as a preprocessing step in conjunction with smile interplation
            // validity of the model is ensured by proceeding it with updateLocalVol()
            for (size_t k = 0; k < Xp_.size(); ++k) { // right wing calculations
                Sp_[k] = underlyingS(true, k, Xp_[k]);
                sigmaP_[k] = localVol(true, k, Sp_[k]);
            }
            for (size_t k = 0; k < Sm_.size(); ++k) { // left wing calculations
                PassiveType x0 = (k > 0) ? Xm_[k - 1] : 0.0;
                Sm_[k] = underlyingS(false, k, Xm_[k]);
                sigmaM_[k] = localVol(false, k, Sm_[k]);
            }
        }

        inline void updateLocalVol() {
            // use ODE solution to determine x-grid and sigma-grid taking into account constraints of
            // positive local volatility and local vol extrapolation
            for (size_t k = 0; k < Sp_.size(); ++k) { // right wing calculations
                PassiveType x0 = (k > 0) ? Xp_[k - 1] : 0.0;
                PassiveType S0 = (k > 0) ? Sp_[k - 1] : S0_;
                PassiveType sigma0 = (k > 0) ? sigmaP_[k - 1] : sigma0_;
                QL_REQUIRE(sigma0 >= 0.0, "sigma0 >= 0.0 required.");
                if ((k == Sp_.size() - 1)||(localVol(true, k, Sp_[k])<=0.0)||(underlyingX(true, k, Sp_[k])>upperBoundX()))  { // right wing extrapolation, maybe better use some epsilon here
                    PassiveType XRight = upperBoundX();  // mu might not yet be calibrated
                    QL_REQUIRE(XRight >= x0, "XRight >= x0 required.");
                    Xp_[k] = XRight;
                    Sp_[k] = underlyingS(true, k, XRight);
                    sigmaP_[k] = localVol(true, k, Sp_[k]);
                    if (k < Sp_.size() - 1) Mp_[k + 1] = Mp_[k];  // we need to make sure vol doesn't go up again
                    continue;
                }
                sigmaP_[k] = localVol(true, k, Sp_[k]);
                QL_REQUIRE(sigmaP_[k] > 0.0, "sigmaP_[k] > 0.0 required.");
                Xp_[k] = underlyingX(true, k, Sp_[k]);
            }
            for (size_t k = 0; k < Sm_.size(); ++k) { // left wing calculations
                PassiveType x0 = (k > 0) ? Xm_[k - 1] : 0.0;
                PassiveType S0 = (k > 0) ? Sm_[k - 1] : S0_;
                PassiveType sigma0 = (k > 0) ? sigmaM_[k - 1] : sigma0_;
                QL_REQUIRE(sigma0 >= 0.0, "sigma0 >= 0.0 required.");
                if ((k == Sm_.size() - 1)||(localVol(false, k, Sm_[k]) <= 0.0)||(underlyingX(false, k, Sm_[k])<lowerBoundX())) { // left wing extrapolation, maybe better use some epsilon here
                    PassiveType XLeft = lowerBoundX();  // mu might not yet be calibrated
                    QL_REQUIRE(XLeft <= x0, "XLeft <= x0 required.");
                    Xm_[k] = XLeft;
                    Sm_[k] = underlyingS(false, k, XLeft);
                    sigmaM_[k] = localVol(false, k, Sm_[k]);
                    if (k < Sm_.size() - 1) Mm_[k + 1] = Mm_[k];  // we need to make sure vol doesn't go up again
                    continue;
                }
                sigmaM_[k] = localVol(false, k, Sm_[k]);
                QL_REQUIRE(sigmaM_[k] > 0.0, "sigmaM_[k] > 0.0 required.");
                Xm_[k] = underlyingX(false, k, Sm_[k]);
            }
        }

        inline void calibrateATM() {
            PassiveType straddleVega = straddleATM_ / sigmaATM_;
            PassiveType forwardMinusStrike0, forwardMinusStrike1, straddleMinusATM0, straddleMinusATM1, dmu, dlogSigma0;
            PassiveType	dfwd_dmu, dstr_dlogSigma0, logSigma0=log(sigma0_);
            for (size_t k = 0; k < maxCalibrationIters_; ++k) {
                PassiveType call = expectation(true, S0_);
                PassiveType put = expectation(false, S0_);
                forwardMinusStrike1 = call - put;
                straddleMinusATM1 = call + put - straddleATM_;
                if (k > 0) {  // perform line search
                    PassiveType num = forwardMinusStrike0*(forwardMinusStrike1 - forwardMinusStrike0) +
                                      straddleMinusATM0*(straddleMinusATM1 - straddleMinusATM0);
                    PassiveType den = (forwardMinusStrike1 - forwardMinusStrike0)*(forwardMinusStrike1 - forwardMinusStrike0) +
                                      (straddleMinusATM1 - straddleMinusATM0)*(straddleMinusATM1 - straddleMinusATM0);
                    PassiveType lambda = -num / den;
                    PassiveType eps = 1.0e-6;  // see Griewank '86
                    if (lambda < -0.5 - eps) lambda = -0.5;
                    else if (lambda < -eps) lambda = lambda;
                    else if (lambda < 0.0) lambda = -eps;
                    else if (lambda <= eps) lambda = eps;
                    else if (lambda <= 0.5 + eps) lambda = lambda;
                    else lambda = 1.0;
                    if (lambda < 1.0) { // reject the step and calculate a new try
                        // x = x - dx + lambda dx = x + (lambda - 1.0) dx
                        mu_ += (lambda - 1.0) * dmu;
                        logSigma0 += (lambda - 1.0) * dlogSigma0;
                        dmu *= lambda;
                        dlogSigma0 *= lambda;
                        sigma0_ = exp(logSigma0);
                        updateLocalVol();
                        if (enableLogging_) logging_.push_back("k: " + std::to_string(k) +
                            "; C: " + std::to_string(call) +
                            "; P: " + std::to_string(put) +
                            "; S: " + std::to_string(straddleATM_) +
                            "; lambda: " + std::to_string(lambda) +
                            "; dmu: " + std::to_string(dmu) +
                            "; dlogSigma0: " + std::to_string(dlogSigma0));
                        continue;  // don't update derivatives and step direction for rejected steps
                    }
                }
                if (k == 0) {
                    dfwd_dmu = sigma0_;        // this is an estimate based on dS/dX at ATM
                    dstr_dlogSigma0 = straddleVega * sigma0_;   // this is an estimate based on dsigmaATM / dsigma0 =~ 1
                }
                if (k>0) { // we use secant if available
                    // only update derivative if we had a step, otherwise use from previous iteration
                    // also avoid division by zero and zero derivative
                    PassiveType eps = 1.0e-12;  // we aim at beeing a bit more robust
                    if ((fabs(forwardMinusStrike1 - forwardMinusStrike0) > eps) && (fabs(dmu) > eps)) {
                        dfwd_dmu = (forwardMinusStrike1 - forwardMinusStrike0) / dmu;
                    }
                    if ((fabs(straddleMinusATM1 - straddleMinusATM0) > eps) && (fabs(dlogSigma0) > eps)) {
                        dstr_dlogSigma0 = (straddleMinusATM1 - straddleMinusATM0) / dlogSigma0;
                    }
                }
                dmu = -forwardMinusStrike1 / dfwd_dmu;
                if (k < onlyForwardCalibrationIters_) dlogSigma0 = 0.0;  // keep sigma0 fixed and only calibrate forward
                else dlogSigma0 = -straddleMinusATM1 / dstr_dlogSigma0;
                if (dmu <= -0.9*upperBoundX()) dmu = -0.5*upperBoundX();  // make sure 0 < eps < upperBoundX() in next update
                if (dmu >= -0.9*lowerBoundX()) dmu = -0.5*lowerBoundX();  // make sure 0 > eps > lowerBoundX() in next update
                // maybe some line search could improve convergence...
                mu_ += dmu;
                logSigma0 += dlogSigma0;
                sigma0_ = exp(logSigma0);  // ensure sigma0 > 0
                updateLocalVol();
                // prepare for next iteration
                forwardMinusStrike0 = forwardMinusStrike1;
                straddleMinusATM0 = straddleMinusATM1;
                if (enableLogging_) logging_.push_back("k: " + std::to_string(k) +
                    "; C: " + std::to_string(call) +
                    "; P: " + std::to_string(put) +
                    "; S: " + std::to_string(straddleATM_) +
                    "; dfwd_dmu: " + std::to_string(dfwd_dmu) +
                    "; dstr_dlogSigma0: " + std::to_string(dstr_dlogSigma0) +
                    "; dmu: " + std::to_string(dmu) +
                    "; dlogSigma0: " + std::to_string(dlogSigma0));
                if ((fabs(forwardMinusStrike0) < S0Tol_) && (fabs(sigma0_*dlogSigma0) < sigma0Tol_)) break;
            }
        }

        inline void adjustATM() {
            // reset adjusters in case this method is invoked twice
            alpha_ = 1.0;
            nu_ = 0.0;
            PassiveType call0 = expectation(true, S0_);
            PassiveType put0 = expectation(false, S0_);
            nu_ = put0 - call0;
            if (enableLogging_) logging_.push_back("C0: " + std::to_string(call0) + "; P0: " + std::to_string(put0) + "; nu: " + std::to_string(nu_));
            PassiveType call1 = expectation(true, S0_);
            PassiveType put1 = expectation(false, S0_);
            alpha_ = straddleATM_ / (call1 + put1);
            nu_ = alpha_*nu_ + (1.0 - alpha_)*S0_;
            if (enableLogging_) logging_.push_back("C1: " + std::to_string(call1) + "; P1: " + std::to_string(put1) + "; alpha_: " + std::to_string(alpha_) + "; nu_: " + std::to_string(nu_));
        }

    public:
        // construct model based on S-grid
        VanillaLocalVolModelT(
            const DateType                   T,
            const PassiveType                S0,
            const PassiveType                sigmaATM,
            const std::vector<PassiveType>&  Sp,
            const std::vector<PassiveType>&  Sm,
            const std::vector<PassiveType>&  Mp,
            const std::vector<PassiveType>&  Mm,
            // controls for calibration
            const size_t                     maxCalibrationIters = 5,
            const size_t                     onlyForwardCalibrationIters = 0,
            const bool                       adjustATMFlag = true,
            const bool                       enableLogging = false,
            const bool                       useInitialMu  = false,
            const PassiveType                initialMu     = 0.0 )
            : T_(T), S0_(S0), sigmaATM_(sigmaATM), sigma0_(sigmaATM), Sp_(Sp), Sm_(Sm), Mp_(Mp), Mm_(Mm),
            maxCalibrationIters_(maxCalibrationIters), onlyForwardCalibrationIters_(onlyForwardCalibrationIters),
            adjustATM_(adjustATMFlag), enableLogging_(enableLogging), useInitialMu_(useInitialMu), initialMu_(initialMu) {
            // some basic sanity checks come here to avoid the need for taking care of it later on
            QL_REQUIRE(T_ > 0, "T_ > 0 required.");
            QL_REQUIRE(sigmaATM_ > 0, "sigmaATM_ > 0 required.");
            QL_REQUIRE(Sp_.size() > 0, "Sp_.size() > 0 required.");
            QL_REQUIRE(Sm_.size() > 0, "Sm_.size() > 0 required.");
            QL_REQUIRE(Mp_.size() == Sp_.size(), "Mp_.size() == Sp_.size() required.");
            QL_REQUIRE(Mm_.size() == Sm_.size(), "Mm_.size() == Sm_.size() required.");
            // check for monotonicity
            QL_REQUIRE(Sp_[0] > S0_, "Sp_[0] > S0_ required.");
            for (size_t k=1; k<Sp_.size(); ++k) QL_REQUIRE(Sp_[k] > Sp_[k-1], "Sp_[k] > Sp_[k-1] required.");
            QL_REQUIRE(Sm_[0] < S0_, "Sm_[0] < S0_ required.");
            for (size_t k = 1; k<Sm_.size(); ++k) QL_REQUIRE(Sm_[k] < Sm_[k-1], "Sm_[k] < Sm_[k-1] required.");
            // now it makes sense to allocate memory
            sigmaP_.resize(Sp_.size());
            sigmaM_.resize(Sm_.size());
            Xp_.resize(Sp_.size());
            Xm_.resize(Sm_.size());
            // initialize deep-in-the-model parameters
            initializeDeepInTheModelParameters();
            // now we may calculate local volatility
            updateLocalVol();
            calibrateATM();
            if (adjustATM_) adjustATM();
        }

        // construct model based on x-grid
        VanillaLocalVolModelT(
            const DateType                   T,
            const PassiveType                S0,
            const PassiveType                sigmaATM,
            const PassiveType                sigma0,
            const std::vector<PassiveType>&  Xp,
            const std::vector<PassiveType>&  Xm,
            const std::vector<PassiveType>&  Mp,
            const std::vector<PassiveType>&  Mm,
            // controls for calibration
            const size_t                     maxCalibrationIters = 5,
            const size_t                     onlyForwardCalibrationIters = 0,
            const bool                       adjustATMFlag = true,
            const bool                       enableLogging = false,
            const bool                       useInitialMu = false,
            const PassiveType                initialMu = 0.0)
            : T_(T), S0_(S0), sigmaATM_(sigmaATM), sigma0_(sigma0), Xp_(Xp), Xm_(Xm), Mp_(Mp), Mm_(Mm),
            maxCalibrationIters_(maxCalibrationIters), onlyForwardCalibrationIters_(onlyForwardCalibrationIters),
            adjustATM_(adjustATMFlag), enableLogging_(enableLogging), useInitialMu_(useInitialMu), initialMu_(initialMu) {
            // some basic sanity checks come here to avoid the need for taking care of it later on
            QL_REQUIRE(T_ > 0, "T_ > 0 required.");
            QL_REQUIRE(sigmaATM_ > 0, "sigmaATM_ > 0 required.");
            QL_REQUIRE(sigma0_ > 0, "sigma0_ > 0 required.");
            QL_REQUIRE(Xp_.size() > 0, "Xp_.size() > 0 required.");
            QL_REQUIRE(Xm_.size() > 0, "Xm_.size() > 0 required.");
            QL_REQUIRE(Mp_.size() == Xp_.size(), "Mp_.size() == Xp_.size() required.");
            QL_REQUIRE(Mm_.size() == Xm_.size(), "Mm_.size() == Xm_.size() required.");
            // check for monotonicity
            QL_REQUIRE(Xp_[0] > 0.0, "Xp_[0] > 0.0 required.");
            for (size_t k = 1; k<Xp_.size(); ++k) QL_REQUIRE(Xp_[k] > Xp_[k - 1], "Xp_[k] > Xp_[k-1] required.");
            QL_REQUIRE(Xm_[0] < 0.0, "Xm_[0] < 0.0 required.");
            for (size_t k = 1; k<Xm_.size(); ++k) QL_REQUIRE(Xm_[k] < Xm_[k - 1], "Xm_[k] < Xm_[k-1] required.");
            // now it makes sense to allocate memory
            sigmaP_.resize(Xp_.size());
            sigmaM_.resize(Xm_.size());
            Sp_.resize(Xp_.size());
            Sm_.resize(Xm_.size());
            // initialize deep-in-the-model parameters
            initializeDeepInTheModelParameters();
            // now we may calculate local volatility
            calculateSGrid();  // we need this preprocessing step since we only input x instead of S
            updateLocalVol();
            calibrateATM();
            if (adjustATM_) adjustATM();
        }



        // inspectors

        inline const std::vector<std::string> logging() { return logging_; }
        inline const PassiveType timeToExpiry()         { return T_;       }
        inline const PassiveType forward()              { return S0_;      }
        inline const PassiveType sigmaATM()             { return sigmaATM_; }
        inline const PassiveType alpha()                { return alpha_;   }
        inline const PassiveType mu()                   { return mu_;      }
        inline const PassiveType nu()                   { return nu_;      }
        inline const size_t      maxCalibrationIters()  { return maxCalibrationIters_; }
        inline const size_t      onlyForwardCalibrationIters() { return onlyForwardCalibrationIters_; }
        inline const bool        adjustATMFlag()        { return adjustATM_; }
        inline const bool        enableLogging()        { return enableLogging_; }
        inline const bool        useInitialMu()         { return useInitialMu_;  }
        inline const PassiveType initialMu()            { return initialMu_;     }


        // attributes in more convenient single-vector format

        const std::vector<PassiveType> underlyingX() {
            std::vector<PassiveType> X(Xm_.size() + Xp_.size() + 1);
            for (size_t k = 0; k < Xm_.size(); ++k) X[k] = Xm_[Xm_.size() - k - 1];
            X[Xm_.size()] = 0.0;
            for (size_t k = 0; k < Xp_.size(); ++k) X[Xm_.size() + 1 + k] = Xp_[k];
            return X;
        }

        const std::vector<PassiveType> underlyingS() {
            std::vector<PassiveType> S(Sm_.size() + Sp_.size() + 1);
            for (size_t k = 0; k < Sm_.size(); ++k) S[k] = Sm_[Sm_.size() - k - 1];
            S[Sm_.size()] = S0_;
            for (size_t k = 0; k < Sp_.size(); ++k) S[Sm_.size() + 1 + k] = Sp_[k];
            return S;
        }

        const std::vector<PassiveType> localVol() {
            std::vector<PassiveType> sigma(sigmaM_.size() + sigmaP_.size() + 1);
            for (size_t k = 0; k < sigmaM_.size(); ++k) sigma[k] = sigmaM_[sigmaM_.size() - k - 1];
            sigma[sigmaM_.size()] = sigma0_;
            for (size_t k = 0; k < sigmaP_.size(); ++k) sigma[sigmaM_.size() + 1 + k] = sigmaP_[k];
            return sigma;
        }

        const std::vector<PassiveType> localVolSlope() {
            std::vector<PassiveType> m(Mm_.size() + Mp_.size() + 1);
            for (size_t k = 0; k < Mm_.size(); ++k) m[k] = Mm_[Mm_.size() - k - 1];
            m[Mm_.size()] = 0.0;  // undefined
            for (size_t k = 0; k < Mp_.size(); ++k) m[Mm_.size() + 1 + k] = Mp_[k];
            return m;
        }

        // model function evaluations

        const PassiveType localVol(PassiveType S) {
            bool isRightWing = (S >= S0_) ? true : false;
            size_t idx = 0;
            if (isRightWing) while ((idx < Sp_.size() - 1) && (Sp_[idx] < S)) ++idx;
            else             while ((idx < Sm_.size() - 1) && (Sm_[idx] > S)) ++idx;
            return localVol(isRightWing, idx, S);
        }

        const PassiveType underlyingS(PassiveType x) {
            bool isRightWing = (x >= 0.0) ? true : false;
            size_t idx = 0;
            if (isRightWing) while ((idx < Xp_.size() - 1) && (Xp_[idx] < x)) ++idx;
            else             while ((idx < Xm_.size() - 1) && (Xm_[idx] > x)) ++idx;
            return underlyingS(isRightWing, idx, x);
        }

        // calculating expectations - that is the actual purpose of that model

        const PassiveType expectation(bool isRightWing, PassiveType strike) {
            // calculate the forward price of an OTM option
            size_t idx = 0;
            if (isRightWing) {
                QL_REQUIRE(strike >= S0_, "strike >= S0_ required");
                while ((idx < Sp_.size()) && (Sp_[idx] <= strike)) ++idx;  // make sure strike < Sp_[idx]
                if (idx == Sp_.size()) return 0.0;  // we are beyond exrapolation
                PassiveType strikeX = underlyingX(isRightWing, idx, strike);
                PassiveType x0 = (idx > 0) ? Xp_[idx - 1] : 0.0;
                QL_REQUIRE((x0 <= strikeX) && (strikeX <= Xp_[idx]), "(x0 <= strikeX) && (strikeX <= Xp_[idx]) required");
                PassiveType intS = 0.0;
                for (size_t k = idx; k < Sp_.size(); ++k) {
                    PassiveType xStart = (k == idx) ? strikeX : Xp_[k - 1];
                    intS += (primitiveF(isRightWing, k, Xp_[k]) - primitiveF(isRightWing, k, xStart));
                }
                // we need to adjust for the strike integral
                PassiveType xEnd = Xp_.back();
                PassiveType intK = TemplateAuxilliaries::Phi((xEnd - mu_) / sqrt(T_)) - TemplateAuxilliaries::Phi((strikeX - mu_) / sqrt(T_));
                return intS - strike * intK;
            }
            else {
                QL_REQUIRE(strike <= S0_, "strike <= S0_ required");
                while ((idx < Sm_.size()) && (Sm_[idx] >= strike)) ++idx;  // make sure Sm_[idx] < strke
                if (idx == Sm_.size()) return 0.0;  // we are beyond exrapolation
                PassiveType strikeX = underlyingX(isRightWing, idx, strike);
                PassiveType x0 = (idx > 0) ? Xm_[idx - 1] : 0.0;
                QL_REQUIRE((x0 >= strikeX) && (strikeX >= Xm_[idx]), "(x0 >= strikeX) && (strikeX >= Xm_[idx]) required");
                PassiveType intS = 0.0;
                for (size_t k = idx; k < Sm_.size(); ++k) {
                    PassiveType xStart = (k == idx) ? strikeX : Xm_[k - 1];
                    intS += (primitiveF(isRightWing, k, Xm_[k]) - primitiveF(isRightWing, k, xStart));
                }
                // we need to adjust for the strike integral
                PassiveType xEnd = Xm_.back();
                PassiveType intK = TemplateAuxilliaries::Phi((xEnd - mu_) / sqrt(T_)) - TemplateAuxilliaries::Phi((strikeX - mu_) / sqrt(T_));
                return intS - strike * intK;
            }
        }

        const PassiveType variance(bool isRightWing, PassiveType strike) {
            // calculate the forward price of an OTM power option with payoff 1_{S>K}(S-K)^2
            size_t idx = 0;
            if (isRightWing) {
                QL_REQUIRE(strike >= S0_, "strike >= S0_ required");
                while ((idx < Sp_.size()) && (Sp_[idx] <= strike)) ++idx;  // make sure strike < Sp_[idx]
                if (idx == Sp_.size()) return 0.0;  // we are beyond exrapolation
                PassiveType strikeX = underlyingX(isRightWing, idx, strike);
                PassiveType x0 = (idx > 0) ? Xp_[idx - 1] : 0.0;
                QL_REQUIRE((x0 <= strikeX) && (strikeX <= Xp_[idx]), "(x0 <= strikeX) && (strikeX <= Xp_[idx]) required");
                PassiveType intS=0.0, intS2 = 0.0;
                for (size_t k = idx; k < Sp_.size(); ++k) {
                    PassiveType xStart = (k == idx) ? strikeX : Xp_[k - 1];
                    intS  += (primitiveF(isRightWing, k, Xp_[k]) - primitiveF(isRightWing, k, xStart));
                    intS2 += (primitiveFSquare(isRightWing, k, Xp_[k]) - primitiveFSquare(isRightWing, k, xStart));
                }
                // we need to adjust for the Vanilla and strike integral
                PassiveType xEnd = Xp_.back();
                PassiveType intK = TemplateAuxilliaries::Phi((xEnd - mu_) / sqrt(T_)) - TemplateAuxilliaries::Phi((strikeX - mu_) / sqrt(T_));
                return intS2 - 2.0 * strike * intS + strike * strike * intK;
            }
            else {
                QL_REQUIRE(strike <= S0_, "strike <= S0_ required");
                while ((idx < Sm_.size()) && (Sm_[idx] >= strike)) ++idx;  // make sure Sm_[idx] < strke
                if (idx == Sm_.size()) return 0.0;  // we are beyond exrapolation
                PassiveType strikeX = underlyingX(isRightWing, idx, strike);
                PassiveType x0 = (idx > 0) ? Xm_[idx - 1] : 0.0;
                QL_REQUIRE((x0 >= strikeX) && (strikeX >= Xm_[idx]), "(x0 >= strikeX) && (strikeX >= Xm_[idx]) required");
                PassiveType intS = 0.0, intS2 = 0.0;
                for (size_t k = idx; k < Sm_.size(); ++k) {
                    PassiveType xStart = (k == idx) ? strikeX : Xm_[k - 1];
                    intS  += (primitiveF(isRightWing, k, Xm_[k]) - primitiveF(isRightWing, k, xStart));
                    intS2 += (primitiveFSquare(isRightWing, k, Xm_[k]) - primitiveFSquare(isRightWing, k, xStart));
                }
                // we need to adjust for the strike integral
                PassiveType xEnd = Xm_.back();
                PassiveType intK = TemplateAuxilliaries::Phi((xEnd - mu_) / sqrt(T_)) - TemplateAuxilliaries::Phi((strikeX - mu_) / sqrt(T_));
                return -(intS2 - 2.0 * strike * intS + strike * strike * intK);
            }
        }

        // we need to implement the stochastic process interface to use the model in MC simulation

        // dimension of X
        inline virtual size_t size() { return 1; }
        // stochastic factors (underlying, volatilities and spreads)
        inline virtual size_t factors() { return 1; }
        // initial values for simulation
        inline virtual VecP initialValues() { return VecP(1, S0_); }
        // a[t,X(t)]
        inline virtual VecA drift(const DateType t, const VecA& X) { return VecP(1, 0.0); }
        // b[t,X(t)]
        inline virtual MatA diffusion(const DateType t, const VecA& X) { return MatA(1, VecA(1, localVol(X[0]))); }
        // we need to implement evolution based on local shifted-lognormal dynamics
        inline virtual void evolve(const DateType t0, const VecA& X0, const DateType dt, const VecD& dW, VecA& X1) {
            // we do not want to use low-level implementation details because these might change
            // instead we "assume" unknown local shifted lognormal dynamics
            PassiveType S = X0[0], epsilon = 1.0e-6;
            PassiveType sigma = localVol(S), sigma2 = localVol(S + epsilon);
            PassiveType m = (sigma2 - sigma) / epsilon;
            PassiveType dS = sigma * dW[0] * sqrt(dt);  // default for normal model
            if (m != 0.0) dS = (exp(m*dW[0]*sqrt(dt) - m*m/2.0*dt) - 1.0)*sigma/m;
            X1[0] = S + dS;
            truncate(t0 + dt, X1);
            return;
        }
        // the numeraire in the domestic currency used for discounting future payoffs
        inline virtual ActiveType numeraire(const DateType t, const VecA& X) { return 1.0; }
        // an asset with (individual) drift and volatility
        inline virtual ActiveType asset(const DateType t, const VecA& X, const std::string& alias) { return X[0]; }
        inline virtual ActiveType forwardAsset(const DateType t, const DateType T, const VecA& X, const std::string& alias) { return X[0]; }
        inline virtual ActiveType zeroBond(const DateType t, const DateType T, const VecA& X) { return 1.0; }


    };


}


#endif  /* ifndef quantlib_templatevanillalocalvolmodel_hpp_hpp */
