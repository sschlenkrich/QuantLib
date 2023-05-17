/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2015 Sebastian Schlenkrich

*/



#ifndef quantlib_templatetdstochvolmodel_hpp
#define quantlib_templatetdstochvolmodel_hpp

#include <complex>
#include <ql/shared_ptr.hpp>
#include <boost/function.hpp>
#include <ql/errors.hpp>
#include <ql/experimental/templatemodels/auxilliaries/auxilliariesT.hpp>
#include <ql/experimental/templatemodels/auxilliaries/gausslobattoT.hpp>
#include <ql/experimental/templatemodels/auxilliaries/complexT.hpp>
#include <ql/experimental/templatemodels/auxilliaries/solver1dT.hpp>
#include <ql/experimental/templatemodels/stochasticprocessT.hpp>
#include <ql/experimental/templatemodels/stochvol/hestonmodelT.hpp>



#define _MIN_( a, b ) ( (a) < (b) ? (a) : (b) )
#define _MAX_( a, b ) ( (a) > (b) ? (a) : (b) )

namespace QuantLib {

    // general stochastic volatility model interface with time-dependent parameters and parameter averaging pricing formula
    //
    //    dS(t) = lambda(t) [ b(t) S(t) + (1-b(t)) L ] sqrt[z(t)] dW(t)
    //    dz(t) = theta [ m - z(t) ] dt + eta(t) sqrt[z(t)] dZ(t)
    //    dW(t) dZ(t) = rho dt
    //
    template <class DateType, class PassiveType, class ActiveType>
    class TimeDependentStochVolModelT : public StochasticProcessT<DateType,PassiveType,ActiveType> {
    public:
        // abstract inspectors
        virtual ActiveType  lambda( const DateType t) = 0;
        virtual ActiveType  b(      const DateType t) = 0;
        virtual ActiveType  L()                       = 0;
        virtual ActiveType  theta()                   = 0;
        virtual ActiveType  m()                       = 0;
        virtual ActiveType  eta(    const DateType t) = 0;
        virtual ActiveType  z0()                      = 0;
        virtual ActiveType  rho()                     = 0;
        virtual ActiveType  S0()                      = 0;

        // from base class
        using typename StochasticProcessT<DateType, PassiveType, ActiveType>::VecD;
        using typename StochasticProcessT<DateType, PassiveType, ActiveType>::VecP;
        using typename StochasticProcessT<DateType, PassiveType, ActiveType>::VecA;
        using typename StochasticProcessT<DateType, PassiveType, ActiveType>::MatA;
        using typename StochasticProcessT<DateType, PassiveType, ActiveType>::VolEvolv;
        using StochasticProcessT<DateType, PassiveType, ActiveType>::volEvolv;


        // helper functions for vol averaging, Piterbarg, 10.2.4
        inline static ActiveType A_CIR ( ActiveType c1, ActiveType c2, ActiveType z0, ActiveType theta, ActiveType eta, DateType dt ) {
            ActiveType gamma = sqrt((theta*theta + 2.0 * eta*eta * c2));
            ActiveType t1 = theta*z0/eta/eta * (theta + gamma)*dt;
            ActiveType t2 = 1.0 + (theta + gamma + c1 * eta*eta) * (exp(gamma * dt) - 1.0) / 2.0 / gamma;
            QL_REQUIRE(t2>0,"A_CIR: require positive log()-argument");
            return t1 - 2.0*theta*z0/eta/eta*log(t2);
        }

        inline static ActiveType B_CIR ( ActiveType c1, ActiveType c2, ActiveType z0, ActiveType theta, ActiveType eta, DateType dt ) {
            ActiveType gamma = sqrt((theta*theta + 2.0 * eta*eta * c2));
            ActiveType emGdt = exp(-gamma * dt);
            ActiveType numer = (2.0*c2 - theta*c1)*(1.0 - emGdt) + gamma*c1*(1.0 + emGdt);
            ActiveType denum = (theta + gamma + c1*eta*eta) * (1.0 - emGdt) + 2.0*gamma*emGdt;
            return numer / denum;
        }		

        // abstract averaging formula definitions
        inline virtual ActiveType  averageLambda( const DateType T) = 0;
        inline virtual ActiveType  averageB     ( const DateType T) = 0;
        inline virtual ActiveType  averageEta   ( const DateType T) = 0;

        // undiscounted expectation of vanilla payoff
        inline ActiveType vanillaOption(const PassiveType forwardPrice,
                                        const PassiveType strikePrice,
                                        const DateType    term,
                                        const int         callOrPut,
                                        const PassiveType accuracy,
                                        const size_t      maxEvaluations) {
                StochVolModelT<DateType,PassiveType,ActiveType> model(averageLambda(term), averageB(term),L(),theta(),m(),averageEta(term),z0(),rho());
                return model.vanillaOption( forwardPrice, strikePrice, term, callOrPut, accuracy, maxEvaluations );
        }

        // conditional moments of vol process used for z-integration, Piterbarg, 8.3.3.
        // E[ z(T) | z(t) ]
        inline virtual ActiveType expectationZ( DateType t, ActiveType zt, DateType dT ) {
            return z0() + (zt - z0())*exp(-theta()*dT);
        }
        // Var[ z(T) | z(t) ]
        inline virtual ActiveType varianceZ( DateType t, ActiveType zt, DateType dT ) {
            ActiveType expmThDT = exp(-theta()*dT);
            ActiveType onemETDT = 1 - expmThDT;
            ActiveType eta2oThe = eta(t+dT/2.0)*eta(t+dT/2.0)/theta();  // approx eta(t)=eta for s \in [t, t+dT]
            return zt*eta2oThe*expmThDT*onemETDT + z0()*eta2oThe/2.0*onemETDT*onemETDT; 
        }


        // stochastic process interface
        // dimension of X
        inline virtual size_t size() { return 2; }
        // stochastic factors of x and z (maybe distinguish if trivially eta=0)
        inline virtual size_t factors() { return 2; }
        // initial values for simulation
        inline virtual VecP initialValues() {
            VecP X(2);
            X[0] = S0();
            X[1] = z0();
            return X;
        }
        // a[t,X(t)]
        inline virtual VecA drift( const DateType t, const VecA& X) {
            VecA a(2);
            // S-variable drift-less
            a[0] = 0.0;
            // z-variable theta [ m - z(t)^+ ]  (full truncation)
            a[1] = theta()*(m() - ((X[1]>0)?(X[1]):(0.0)));
            return a;
        }
        // b[t,X(t)]
        inline virtual MatA diffusion( const DateType t, const VecA& X) {
            MatA B(2);
            B[0].resize(2);
            B[1].resize(2);
            ActiveType sqrtz = ( (X[1]>0) ? (sqrt(X[1])) : (0.0) );   // full truncation
            // S-variable lambda(t) [ b(t) S(t) + (1-b(t)) L ] sqrt[z(t)] dW(t)
            B[0][0] = lambda(t) * (b(t) * X[0] + (1.0-b(t)) * L()) * sqrtz;
            B[0][1] = 0.0;
            // z-variable
            B[1][0] = rho()*eta(t)*sqrtz;
            B[1][1] = sqrt(1-rho()*rho())*eta(t)*sqrtz;
            // finished
            return B;
        }

        // integrate X1 = X0 + drift()*dt + diffusion()*dW*sqrt(dt)
        inline virtual void evolve( const DateType t0, const VecA& X0, const DateType dt, const VecD& dW, VecA& X1 ) {
            // ensure X1 has size of X0
            VecA a = drift(t0, X0);
            MatA b = diffusion(t0, X0);
            // S-variable
            X1[0] = X0[0] + a[0]*dt + b[0][0]*dW[0]*sqrt(dt);
            // z-variable
            if (volEvolv()==VolEvolv::FullTruncation) {
                X1[1] = X0[1] + a[1]*dt + (b[1][0]*dW[0]+b[1][1]*dW[1])*sqrt(dt);
                return;
            }
            if (volEvolv()==VolEvolv::LogNormalApproximation) {
                ActiveType e = expectationZ(t0, X0[1], dt);
                ActiveType v = varianceZ(t0, X0[1], dt);
                ActiveType dZ = rho()*dW[0] + sqrt(1-rho()*rho())*dW[1];
                ActiveType si = sqrt(log(1.0 + v/e/e));
                ActiveType mu = log(e) - si*si/2.0;
                X1[1] = exp(mu + si*dZ);
            }
            return;
        }



        // embeded classes

        class PieceWiseConstant  {
        protected:
            // check for dimensions
            bool isConsistent_;
            // time grid
            std::vector<DateType>     times_;
            // time-dependent model parameters
            std::vector<ActiveType>   lambda_;
            std::vector<ActiveType>   b_;
            std::vector<ActiveType>   eta_;
            // time-homogeneous model parameters
            ActiveType                L_;
            ActiveType                theta_;
            ActiveType                m_;
            ActiveType                z0_;
            ActiveType                rho_;
            ActiveType                S0_;
            // helper
            inline size_t idx( const DateType t ) { return TemplateAuxilliaries::idx(times_,t); }
        public:
            // constructor
            PieceWiseConstant( const std::vector<DateType>&    times,
                               const std::vector<ActiveType>&  lambda,
                               const std::vector<ActiveType>&  b,
                               const std::vector<ActiveType>&  eta,
                               const ActiveType                L,
                               const ActiveType                theta,
                               const ActiveType                m,
                               const ActiveType                z0,
                               const ActiveType                rho,
                               const ActiveType                S0 )
                : times_(times), lambda_(lambda), b_(b), eta_(eta), L_(L), theta_(theta), m_(m), z0_(z0), rho_(rho), S0_(S0) {
                isConsistent_ = true;
                QL_REQUIRE(times_.size()>0,"TemplateTimeDependentStochVolModel::PieceWiseConstant: non-empty times required");
                for (size_t k=0; k<times_.size()-1; ++k) QL_REQUIRE(times_[k]<times_[k+1],"TemplateTimeDependentStochVolModel::PieceWiseConstant: ascending time-grid required");
                QL_REQUIRE(lambda_.size()==times_.size(),"TemplateTimeDependentStochVolModel::PieceWiseConstant: lambda dimension mismatch");
                QL_REQUIRE(b_.size()     ==times_.size(),"TemplateTimeDependentStochVolModel::PieceWiseConstant: b dimension mismatch");
                QL_REQUIRE(eta_.size()   ==times_.size(),"TemplateTimeDependentStochVolModel::PieceWiseConstant: eta dimension mismatch");
            }
            // inspectors
            virtual ActiveType  lambda( const DateType t)  { return lambda_[idx(t)]; }
            virtual ActiveType  b(      const DateType t)  { return b_[idx(t)];      }
            virtual ActiveType  eta(    const DateType t)  { return eta_[idx(t)];    }
            virtual ActiveType  L()                        { return L_;              }
            virtual ActiveType  theta()                    { return theta_;          }
            virtual ActiveType  m()                        { return m_;              }
            virtual ActiveType  z0()                       { return z0_;             }
            virtual ActiveType  rho()                      { return rho_;            }
            virtual ActiveType  S0()                       { return S0_;            }

        };

        // averaging assuming piecewise constant parameters on ( T_k-1, T_k ]
        class MidPointIntegration {
        protected:
            // reference to model
            TimeDependentStochVolModelT* model_;
            // time grid for integration
            std::vector<DateType>               times_;
            inline std::vector<DateType> getTimes( const DateType T ) {
                std::vector<DateType> times;
                times.push_back(0.0);
                for (size_t k=0; k<times_.size(); ++k) {
                    if (times_[k]>0.0 && times_[k]<T) times.push_back(times_[k]);
                    if (times_[k]>=T) break;
                }
                times.push_back(T);
                return times;
            }

            class AverageLambdaObjective {
            protected:
                ActiveType z0_;
                ActiveType theta_;
                ActiveType eta_;
                DateType   dt_;
                ActiveType target_;
            public:
                AverageLambdaObjective ( const ActiveType z0,
                                         const ActiveType theta,
                                         const ActiveType eta,
                                         const DateType   dt,
                                         const ActiveType target )
                                         : z0_(z0), theta_(theta), eta_(eta), dt_(dt), target_(target) {}
                ActiveType operator() ( const ActiveType avLambda2c ) {
                    ActiveType A = A_CIR(0, -avLambda2c, z0_, theta_, eta_, dt_);
                    ActiveType B = B_CIR(0, -avLambda2c, z0_, theta_, eta_, dt_);
                    ActiveType res = A - B*z0_ - target_;
                    return res;
                }
            };

            // Prop. 9.1.2
            class RiccatiODE {
                TimeDependentStochVolModelT* model_;
                ActiveType v_;
                ActiveType u_;
                ActiveType averageB_;
            public:
                RiccatiODE (TimeDependentStochVolModelT* model, ActiveType v, ActiveType u, ActiveType averageB)
                    : model_(model), v_(v), u_(u), averageB_(averageB) {}
                void operator() (const DateType t, const std::vector<ActiveType>& y, std::vector<ActiveType>& fy) {
                    fy[0] = - model_->theta() * model_->z0() * y[1];
                    fy[1] = (model_->theta() - model_->rho()*model_->eta(t)*averageB_*u_*model_->lambda(t))*y[1]
                            - model_->eta(t)*model_->eta(t)/2.0*y[1]*y[1] - v_*model_->lambda(t)*model_->lambda(t);
                }
            };

        public:
            // constructor
            MidPointIntegration ( TimeDependentStochVolModelT* model,
                                  std::vector<DateType>        times )
                                  : model_(model), times_(times) { }

            virtual ~MidPointIntegration() = default;

            // averaging formula implementations
            virtual ActiveType  averageEta   ( const DateType T) {
                std::vector<DateType> times(getTimes(T));
                std::vector<ActiveType> f(times.size());
                std::vector<ActiveType> w(times.size()-1);
                f[f.size()-1] = 0.0;
                for (size_t k=f.size()-1; k>0; --k) {
                    ActiveType lambda = model_->lambda((times[k-1]+times[k])/2.0);
                    ActiveType tmp    = exp(-model_->theta()*(times[k]-times[k-1]));
                    f[k-1] = lambda*lambda/model_->theta()*(1.0 - tmp) + tmp*f[k];
                }
                ActiveType sum = 0.0;
                for (size_t k=0; k<w.size(); ++k) {
                    ActiveType lambda  = model_->lambda((times[k]+times[k+1])/2.0);
                    ActiveType lambda2 = lambda*lambda;
                    ActiveType theta   = model_->theta();
                    ActiveType theta2  = theta*theta;
                    w[k] = (f[k+1]*f[k+1] - f[k]*f[k])/2.0/theta +
                           (f[k+1] - f[k])*lambda2/theta2 +
                           (times[k+1]-times[k])*lambda2*lambda2/theta2;
                    w[k] *= 0.5;
                    sum  += w[k];
                }
                ActiveType eta2=0.0;
                for (size_t k=0; k<w.size(); ++k) {
                    ActiveType eta = model_->eta((times[k]+times[k+1])/2.0);
                    eta2 += w[k] * eta * eta;
                }
                eta2 = eta2 / sum;
                return sqrt(eta2);
            }

            virtual ActiveType  averageB( const DateType T) {
                std::vector<DateType> times(getTimes(T));
                std::vector<ActiveType> w(times.size()-1);
                ActiveType theta = model_->theta();
                ActiveType z0 = model_->z0();
                ActiveType sumLambda2dT = 0.0;
                ActiveType S1=0.0, S2=0.0, S3=0.0;
                ActiveType sum = 0.0;
                for (size_t k=0; k<w.size(); ++k) {
                    // v1
                    ActiveType lambda    = model_->lambda((times[k]+times[k+1])/2.0);
                    ActiveType lambda2   = lambda*lambda;
                    ActiveType eta       = model_->eta((times[k]+times[k+1])/2.0);
                    ActiveType lambda2dT = lambda2*(times[k+1]-times[k]);
                    ActiveType v1        = z0*z0*(times[k+1]-times[k])*(lambda2dT/2.0 + sumLambda2dT);
                    sumLambda2dT        += lambda2dT;
                    // v3, v4, v5
                    ActiveType expmThdT  = exp(-theta*(times[k+1]-times[k]));
                    ActiveType v3        = z0/theta*(S2+S3)*(1.0-expmThdT);
                    ActiveType v4        = (times[k+1]-times[k]) - (1.0-expmThdT)/theta -
                                           (1.0-expmThdT)*(1.0-expmThdT)/2.0/theta;
                    ActiveType laEtaTh   = lambda*eta/theta;
                    v4                  *= z0*laEtaTh*laEtaTh/2.0;
                    ActiveType theta2    = theta*theta;
                    ActiveType v5        = z0*lambda2/theta2/2.0*S1*(1.0-expmThdT)*(1.0-expmThdT);
                    // updating S1, S2, S3
                    S3                   = expmThdT * ( S3 + S1*lambda2/theta*(1.0-expmThdT) );
                    S2                   = expmThdT*S2 + laEtaTh*laEtaTh/2.0*(1.0-expmThdT)*(1.0-expmThdT);
                    S1                   = expmThdT*expmThdT*S1 + eta*eta/theta/2.0*(1.0 - expmThdT*expmThdT);
                    // gathering things together...
                    w[k]                 = lambda2 * ( v1 + v3 + v4 + v5 );
                    sum                 += w[k];
                }
                ActiveType b=0;
                for (size_t k=0; k<w.size(); ++k) {
                    b += w[k] * model_->b((times[k]+times[k+1])/2.0);
                }
                b = b / sum;
                return b;
            }
            
            virtual ActiveType  averageLambda     ( const DateType T) {
                ActiveType b   = averageB(T);
                ActiveType eta = averageEta(T);  // maybe better use time-dep eta
                std::vector<DateType> times(getTimes(T));
                // c = h''(zeta) / h'(zeta)
                ActiveType zeta = 0.0;
                for (size_t k=0; k<times.size()-1; ++k) {
                    ActiveType lambda    = model_->lambda((times[k]+times[k+1])/2.0);
                    zeta += lambda*lambda*(times[k+1]-times[k]);
                }
                ActiveType avLambda2 = zeta / (times[times.size()-1] - times[0]);
                zeta *= model_->z0();
                ActiveType c = -(b*b/4.0 + 1.0/zeta)/2.0;
                // Psi_{z lambda^2}
                ActiveType A = 0.0, B = 0.0;
                //std::vector<ActiveType> y1(2,0.0), y0(2,0.0);
                //RiccatiODE ode(model_,c,0.0,b);
                for (size_t k=times.size()-1; k>0; --k) {
                    DateType  t = (times[k]+times[k-1])/2.0;
                    DateType dt = (times[k]-times[k-1]);
                    ActiveType lambda = model_->lambda(t);
                    A = A + A_CIR(B, -c*lambda*lambda, model_->z0(), model_->theta(), eta, dt);
                    B = B_CIR(B, -c*lambda*lambda, model_->z0(), model_->theta(), eta, dt);
                    // Riccati ODE via Runge Kutta Method
                    //y1 = y0;
                    //TemplateAuxilliaries::rungeKuttaStep<DateType,ActiveType>(y1,times[k],ode,times[k-1]-times[k],y0);
                }
                ActiveType target = A - B*model_->z0();
                //ActiveType target = y0[0] + y0[1]*model_->z0();
                AverageLambdaObjective f(model_->z0(), model_->theta(), eta, T, target);
                ActiveType avLambda2c = avLambda2 * c;
                avLambda2c = TemplateAuxilliaries::solve1d<ActiveType>(f, 1.0e-8, avLambda2c, avLambda2c, 10);
                ActiveType avLambda = sqrt(avLambda2c / c);
                return avLambda;
            };
        };
    };

    // piecewise constant parameters and numerical integration
    template <class DateType, class PassiveType, class ActiveType>
    class PWCAnalytical : public TimeDependentStochVolModelT<DateType,PassiveType,ActiveType> {
        // from base class
        using typename TimeDependentStochVolModelT<DateType, PassiveType, ActiveType>::PieceWiseConstant;
        using typename TimeDependentStochVolModelT<DateType, PassiveType, ActiveType>::MidPointIntegration;
        using typename TimeDependentStochVolModelT<DateType, PassiveType, ActiveType>::VolEvolv;

    private:
        ext::shared_ptr<PieceWiseConstant>    pwc_;
        ext::shared_ptr<MidPointIntegration>  mp_;
        VolEvolv volEvolv_;
    public:
        PWCAnalytical(const std::vector<DateType>&    times,
                        const std::vector<ActiveType>&  lambda,
                        const std::vector<ActiveType>&  b,
                        const std::vector<ActiveType>&  eta,
                        const ActiveType                L,
                        const ActiveType                theta,
                        const ActiveType                m,
                        const ActiveType                z0,
                        const ActiveType                rho,
                        const ActiveType                S0,
                        const VolEvolv                  volEvolv = VolEvolv::FullTruncation )
            :    pwc_(new PieceWiseConstant(times,lambda,b,eta,L,theta,m,z0,rho,S0)),
                    mp_(new MidPointIntegration(this,times)), volEvolv_(volEvolv) {}						
        // inspectors
        virtual ActiveType  lambda( const DateType t)  { return pwc_->lambda(t);    }
        virtual ActiveType  b(      const DateType t)  { return pwc_->b(t) ;        }
        virtual ActiveType  eta(    const DateType t)  { return pwc_->eta(t);       }
        virtual ActiveType  L()                        { return pwc_->L();          }
        virtual ActiveType  theta()                    { return pwc_->theta();      }
        virtual ActiveType  m()                        { return pwc_->m();          }
        virtual ActiveType  z0()                       { return pwc_->z0();         }
        virtual ActiveType  rho()                      { return pwc_->rho();        }
        virtual ActiveType  S0()                       { return pwc_->S0();         }
        virtual VolEvolv    volEvolv()                 { return volEvolv_;          }
        // averaging formula implementations
        virtual ActiveType  averageLambda( const DateType T) { return mp_->averageLambda( T ); }
        virtual ActiveType  averageB     ( const DateType T) { return mp_->averageB( T );      }
        virtual ActiveType  averageEta   ( const DateType T) { return mp_->averageEta( T );    }
    };

}

#undef _MIN_
#undef _MAX_

#endif  /* ifndef quantlib_templatetdstochvolmodel_hpp */
