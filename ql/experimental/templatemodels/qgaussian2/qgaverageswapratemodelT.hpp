/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2018 Sebastian Schlenkrich

*/

/*! \file qgaverageswapratemodelT.hpp
    \brief (Approximate) swaption pricing for multi-factor quasi-Gaussian model with stochastic vol
           derive effective time-homogenous parameters via averaging
               
               dS(t) = [ sigma + slope (S - S0) ] sqrt(z) dW
               dz(t) = theta [z0 - z(t)]dt + eta sqrt(z) dZ
                z(0) = z0 = 1,  dW dZ = 0

*/


#ifndef quantlib_templateqgaverageswapratemodelT_hpp
#define quantlib_templateqgaverageswapratemodelT_hpp

#include <ql/shared_ptr.hpp>

#include <ql/experimental/templatemodels/stochvol/hestonmodelT.hpp>
#include <ql/experimental/templatemodels/qgaussian2/qgswapratemodelT.hpp>



namespace QuantLib {

    // Declaration of the quasi-Gaussian model class
    template <class DateType, class PassiveType, class ActiveType>
    class QGAverageSwaprateModelT : public QGSwaprateModelT<DateType, PassiveType, ActiveType> {
    protected:
        // average parameters
        ActiveType sigma_;
        ActiveType slope_;
        ActiveType eta_;

        // Vanilla option pricing
        enum { Heston, ShiftedLogNormal, Normal, StochVolNormal } type_;
        ext::shared_ptr< HestonModelT<DateType, PassiveType, ActiveType> > hestonModel_;
        ActiveType S0_;
        ActiveType shift_;

        inline ActiveType slopeOverSigma(const DateType t) { return QGSwaprateModelT<DateType,PassiveType,ActiveType>::slope(t) / QGSwaprateModelT<DateType,PassiveType,ActiveType>::sigma(t); }

        // helper functions for vol averaging, Piterbarg, 10.2.4
        inline static ActiveType A_CIR(ActiveType c1, ActiveType c2, ActiveType z0, ActiveType theta, ActiveType eta, DateType dt) {
            ActiveType gamma = sqrt((theta*theta + 2.0 * eta*eta * c2));
            ActiveType t1 = theta*z0 / eta / eta * (theta + gamma)*dt;
            ActiveType t2 = 1.0 + (theta + gamma + c1 * eta*eta) * (exp(gamma * dt) - 1.0) / 2.0 / gamma;
            QL_REQUIRE(t2>0, "QGAverageSwaprateModelT: A_CIR: require positive log()-argument");
            return t1 - 2.0*theta*z0 / eta / eta*log(t2);
        }

        inline static ActiveType B_CIR(ActiveType c1, ActiveType c2, ActiveType z0, ActiveType theta, ActiveType eta, DateType dt) {
            ActiveType gamma = sqrt((theta*theta + 2.0 * eta*eta * c2));
            ActiveType emGdt = exp(-gamma * dt);
            ActiveType numer = (2.0*c2 - theta*c1)*(1.0 - emGdt) + gamma*c1*(1.0 + emGdt);
            ActiveType denum = (theta + gamma + c1*eta*eta) * (1.0 - emGdt) + 2.0*gamma*emGdt;
            return numer / denum;
        }

        // functor for volatility averaging solver
        class AverageLambdaObjective {
        protected:
            ActiveType z0_;
            ActiveType theta_;
            ActiveType eta_;
            DateType   dt_;
            ActiveType target_;
        public:
            AverageLambdaObjective(const ActiveType z0,
                const ActiveType theta,
                const ActiveType eta,
                const DateType   dt,
                const ActiveType target)
                : z0_(z0), theta_(theta), eta_(eta), dt_(dt), target_(target) {}
            ActiveType operator() (const ActiveType avLambda2c) {
                ActiveType A = A_CIR(0, -avLambda2c, z0_, theta_, eta_, dt_);
                ActiveType B = B_CIR(0, -avLambda2c, z0_, theta_, eta_, dt_);
                ActiveType res = A - B*z0_ - target_;
                return res;
            }
        };


    public:

        // constructor
        QGAverageSwaprateModelT(
            const ext::shared_ptr< QGSwaprateModelT<DateType,PassiveType,ActiveType> >&   model)
            : QGSwaprateModelT<DateType,PassiveType,ActiveType>(*model) {
            // set up averaging
            eta_ = averageEta();
            QL_REQUIRE(eta_>=0.0, "QGAverageSwaprateModelT: eta >= 0 required.");
            ActiveType b = averageSlopeOverSigma();
            sigma_ = averageSigma(eta_, b);
            QL_REQUIRE(sigma_ >= 0.0, "QGAverageSwaprateModelT: sigma >= 0 required.");  // maybe we need strict > 0
            slope_ = sigma_*b;
            // set up Vanilla option pricing
            if (eta_ < 1.0e-2) {  // vol-of-vol below 1% is effectively deterministic
                if (slope_ < 1.0e-4) type_ = Normal;             // we can not really handle negative slopes, therefore also normal if negative
                else                 type_ = ShiftedLogNormal;
            }
            else {
                if (slope_ < 1.0e-6) type_ = StochVolNormal;     // we are a bit more agressive here because this is not implemented yet
                else                 type_ = Heston;
            }
            QL_REQUIRE(type_ != StochVolNormal, "QGAverageSwaprateModelT: StochVolNormal not implemented.");
            S0_ = QGSwaprateModelT<DateType,PassiveType,ActiveType>::S0();
            if ((type_ == ShiftedLogNormal) || (type_ == Heston)) shift_ = sigma_ / slope_ - S0_;
            if (type_ == Heston) {
                hestonModel_ = ext::shared_ptr< HestonModelT<DateType, PassiveType, ActiveType> >(
                    new HestonModelT<DateType, PassiveType, ActiveType>(
                        // state transformations ~S(t) = S(t) + shift, v(t) = z(t) slope^2
                        QGSwaprateModelT<DateType,PassiveType,ActiveType>::theta(),               // kappa
                        QGSwaprateModelT<DateType,PassiveType,ActiveType>::z0()*slope_*slope_,    // theta
                        eta_*slope_,                                                              // sigma
                        QGSwaprateModelT<DateType,PassiveType,ActiveType>::rho(),                 // rho
                        QGSwaprateModelT<DateType,PassiveType,ActiveType>::z0()*slope_*slope_     // v0
                        ));
            }
        }

        virtual ~QGAverageSwaprateModelT() = default;

        inline ActiveType averageEta() {
            std::vector<DateType> times(QGSwaprateModelT<DateType,PassiveType,ActiveType>::modelTimes());
            std::vector<ActiveType> f(times.size());
            std::vector<ActiveType> w(times.size() - 1);
            f[f.size() - 1] = 0.0;
            for (size_t k = f.size() - 1; k>0; --k) {
                ActiveType sigma = QGSwaprateModelT<DateType,PassiveType,ActiveType>::sigma((times[k - 1] + times[k]) / 2.0);
                ActiveType tmp = exp(-QGSwaprateModelT<DateType,PassiveType,ActiveType>::theta()*(times[k] - times[k - 1]));
                f[k - 1] = sigma*sigma / QGSwaprateModelT<DateType,PassiveType,ActiveType>::theta()*(1.0 - tmp) + tmp*f[k];
            }
            ActiveType sum = 0.0;
            for (size_t k = 0; k<w.size(); ++k) {
                ActiveType sigma = QGSwaprateModelT<DateType,PassiveType,ActiveType>::sigma((times[k] + times[k + 1]) / 2.0);
                ActiveType sigma2 = sigma*sigma;
                ActiveType theta = QGSwaprateModelT<DateType,PassiveType,ActiveType>::theta();
                ActiveType theta2 = theta*theta;
                w[k] = (f[k + 1] * f[k + 1] - f[k] * f[k]) / 2.0 / theta +
                    (f[k + 1] - f[k])*sigma2 / theta2 +
                    (times[k + 1] - times[k])*sigma2*sigma2 / theta2;
                w[k] *= 0.5;
                sum += w[k];
            }
            ActiveType eta2 = 0.0;
            for (size_t k = 0; k<w.size(); ++k) {
                ActiveType eta = QGSwaprateModelT<DateType,PassiveType,ActiveType>::eta((times[k] + times[k + 1]) / 2.0);
                eta2 += w[k] * eta * eta;
            }
            eta2 = eta2 / sum;
            return sqrt(eta2);
        }

        ActiveType averageSlopeOverSigma() {
            std::vector<DateType> times(QGSwaprateModelT<DateType,PassiveType,ActiveType>::modelTimes());
            std::vector<ActiveType> w(times.size() - 1);
            ActiveType theta = QGSwaprateModelT<DateType,PassiveType,ActiveType>::theta();
            ActiveType z0 = QGSwaprateModelT<DateType,PassiveType,ActiveType>::z0();
            ActiveType sumSigma2dT = 0.0;
            ActiveType S1 = 0.0, S2 = 0.0, S3 = 0.0;
            ActiveType sum = 0.0;
            for (size_t k = 0; k<w.size(); ++k) {
                // v1
                ActiveType sigma = QGSwaprateModelT<DateType,PassiveType,ActiveType>::sigma((times[k] + times[k + 1]) / 2.0);
                ActiveType sigma2 = sigma*sigma;
                ActiveType eta = QGSwaprateModelT<DateType,PassiveType,ActiveType>::eta((times[k] + times[k + 1]) / 2.0);
                ActiveType sigma2dT = sigma2*(times[k + 1] - times[k]);
                ActiveType v1 = z0*z0*(times[k + 1] - times[k])*(sigma2dT / 2.0 + sumSigma2dT);
                sumSigma2dT += sigma2dT;
                // v3, v4, v5
                ActiveType expmThdT = exp(-theta*(times[k + 1] - times[k]));
                ActiveType v3 = z0 / theta*(S2 + S3)*(1.0 - expmThdT);
                ActiveType v4 = (times[k + 1] - times[k]) - (1.0 - expmThdT) / theta -
                    (1.0 - expmThdT)*(1.0 - expmThdT) / 2.0 / theta;
                ActiveType siEtaTh = sigma*eta / theta;
                v4 *= z0*siEtaTh*siEtaTh / 2.0;
                ActiveType theta2 = theta*theta;
                ActiveType v5 = z0*sigma2 / theta2 / 2.0*S1*(1.0 - expmThdT)*(1.0 - expmThdT);
                // updating S1, S2, S3
                S3 = expmThdT * (S3 + S1*sigma2 / theta*(1.0 - expmThdT));
                S2 = expmThdT*S2 + siEtaTh*siEtaTh / 2.0*(1.0 - expmThdT)*(1.0 - expmThdT);
                S1 = expmThdT*expmThdT*S1 + eta*eta / theta / 2.0*(1.0 - expmThdT*expmThdT);
                // gathering things together...
                w[k] = sigma2 * (v1 + v3 + v4 + v5);
                sum += w[k];
            }
            ActiveType b = 0;
            for (size_t k = 0; k<w.size(); ++k) {
                b += w[k] * slopeOverSigma((times[k] + times[k + 1]) / 2.0);
            }
            b = b / sum;
            return b;
        }

        ActiveType averageSigma(ActiveType eta, ActiveType b) {
            // ActiveType b = averageSlopeOverSigma(T);  // better avoid calculation twice
            // ActiveType eta = averageEta(T);           // 
            std::vector<DateType> times(QGSwaprateModelT<DateType,PassiveType,ActiveType>::modelTimes());
            // c = h''(zeta) / h'(zeta)
            ActiveType zeta = 0.0;
            for (size_t k = 0; k<times.size() - 1; ++k) {
                ActiveType sigma = QGSwaprateModelT<DateType,PassiveType,ActiveType>::sigma((times[k] + times[k + 1]) / 2.0);
                zeta += sigma*sigma*(times[k + 1] - times[k]);
            }
            ActiveType avSigma2 = zeta / (times[times.size() - 1] - times[0]);
            zeta *= QGSwaprateModelT<DateType,PassiveType,ActiveType>::z0();
            // if there is no skew or no stoch vol we are save to use average sigma
            if ((b == 0.0)||(eta==0.0)) return sqrt(avSigma2);  // maybe define a sutable threshold
            ActiveType c = -(b*b / 4.0 + 1.0 / zeta) / 2.0; 
            // Psi_{z lambda^2}
            ActiveType A = 0.0, B = 0.0;
            for (size_t k = times.size() - 1; k>0; --k) {
                DateType  t = (times[k] + times[k - 1]) / 2.0;
                DateType dt = (times[k] - times[k - 1]);
                ActiveType sigma = QGSwaprateModelT<DateType,PassiveType,ActiveType>::sigma(t);
                A = A + A_CIR(B, -c*sigma*sigma, QGSwaprateModelT<DateType,PassiveType,ActiveType>::z0(), QGSwaprateModelT<DateType,PassiveType,ActiveType>::theta(), eta, dt);
                B = B_CIR(B, -c*sigma*sigma, QGSwaprateModelT<DateType,PassiveType,ActiveType>::z0(), QGSwaprateModelT<DateType,PassiveType,ActiveType>::theta(), eta, dt);
            }
            ActiveType target = A - B*QGSwaprateModelT<DateType,PassiveType,ActiveType>::z0();
            AverageLambdaObjective f(QGSwaprateModelT<DateType,PassiveType,ActiveType>::z0(), QGSwaprateModelT<DateType,PassiveType,ActiveType>::theta(), eta, times[times.size() - 1] - times[0], target);
            ActiveType avSigma2c = avSigma2 * c;
            avSigma2c = TemplateAuxilliaries::solve1d<ActiveType>(f, 1.0e-8, avSigma2c, avSigma2c, 10);
            ActiveType avSigma = sqrt(avSigma2c / c);
            return avSigma;
        }

        // the time-indexed inspectors are unsafe in the sense that they don't check if t > floatTimes[0]

        // overwrite with averaging
        inline virtual ActiveType sigma(const DateType t) { return sigma();  }
            
        // overwrite with averaging
        inline virtual ActiveType slope(const DateType t) { return slope(); }

        // overwrite with averaging
        inline virtual ActiveType eta(const DateType t)   { return eta(); }

        inline virtual ActiveType sigma() { return sigma_; }
        inline virtual ActiveType slope() { return slope_; }
        inline virtual ActiveType eta()   { return eta_;   }

        // undiscounted expectation of vanilla payoff
        inline ActiveType vanillaOption(const PassiveType strikePrice, const int callOrPut, const PassiveType accuracy = 1.0-6, const size_t maxEvaluations = 1000) {
            const std::vector<DateType>& times = QGSwaprateModelT<DateType,PassiveType,ActiveType>::modelTimes();
            DateType term = times[times.size() - 1] - times[0];
            if (type_ == Heston)
                return hestonModel_->vanillaOption(S0_ + shift_, strikePrice + shift_, term, callOrPut, accuracy, maxEvaluations);
            if (type_ == ShiftedLogNormal)
                return TemplateAuxilliaries::Black76(S0_ + shift_, strikePrice + shift_, slope_, term, callOrPut);
            if (type_ == Normal)
                return TemplateAuxilliaries::Bachelier(S0_, strikePrice, sigma_, term, callOrPut);
            QL_REQUIRE(false, "QGAverageSwaprateModelT: unknown model type.");
            return 0;
        }

    };


}

#endif  /* ifndef quantlib_templateqgaverageswapratemodelT_hpp */
