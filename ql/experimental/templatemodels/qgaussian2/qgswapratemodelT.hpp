/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2018 Sebastian Schlenkrich

*/

/*! \file qgswapratemodelT.hpp
    \brief (Approximate) swaption pricing for multi-factor quasi-Gaussian model with stochastic vol
               
               dS(t) = [ sigma + slope (S - S0) ] sqrt(z) dW
               dz(t) = theta [z0 - z(t)]dt + eta sqrt(z) dZ
                z(0) = z0 = 1,  dW dZ = 0

               generalised swap rate

               S(t)  =  [ sum_j u_j P(t,T_j) ] / [ sum_i w_i P(t,T_i) ]

               single-curve, annual fixed payment use case is
               u_1 = 1, u_2 = -1 (float leg), w_i = 1 (annuity with annual year fractions)

               multi-curve, semi-annual fixed payment use case is
               u_1 = 1+spr_1, u_j = spr_j (j=2..M-1), u_M = -1 (float leg with deterministic tenor basis spreads)
               w_i = 0.5 (annuity with semi-annual year fractions)
                           
           All methods are template based to allow incorporation of Automatic Differentiation
           tools
*/


#ifndef quantlib_templateqgswapratemodelT_hpp
#define quantlib_templateqgswapratemodelT_hpp

#include <ql/shared_ptr.hpp>

#include <ql/experimental/templatemodels/qgaussian2/quasigaussianmodel2T.hpp>


namespace QuantLib {

    // Declaration of the quasi-Gaussian model class
    template <class DateType, class PassiveType, class ActiveType>
    class QGSwaprateModelT : public StochasticProcessT<DateType, PassiveType, ActiveType> {
    protected:

        // container class definitions
        typedef std::vector<DateType>                      VecD;
        typedef std::vector<PassiveType>                   VecP; 
        typedef std::vector<ActiveType>                    VecA;
        typedef std::vector< std::vector<DateType> >       MatD;
        typedef std::vector< std::vector<PassiveType> >    MatP;
        typedef std::vector< std::vector<ActiveType> >     MatA;

        // swap spec container
        struct Swap {
            // float leg
            VecD   floatTimes;    // T[0], ..., T[M]
            VecD   floatWeights;  // u[0], ..., u[M]
            // fixed (annuity) leg
            VecD   fixedTimes;    // T[0], ..., T[N]
            VecD   fixedWeights;  // w[0], ..., w[N]
            Swap ( const VecD& floatT, const VecD& floatW, const VecD& fixedT, const VecD& fixedW)
                : floatTimes(floatT), floatWeights(floatW), fixedTimes(fixedT), fixedWeights(fixedW) {}
        };

        // reference to QG modelspecs
        ext::shared_ptr< QuasiGaussianModel2T<DateType,PassiveType,ActiveType> >   model_;

        // vanilla swap specs for this model
        Swap               swap_;

        VecD               times_;    // time grid for discretisation and time-domain specification, t_0=0, t_N=T_exercise
        std::vector<VecA>  barX_;     // E^A [ x(t) ] 
        std::vector<MatA>  barY_;     // E^A [ y(t) ] 
        VecA               x0_;       // = VecA(model_->factors()-1,0.0);
        MatA               y0_;       // = MatA(model_->factors()-1,x0_);
        bool               useExpectedXY_;
        // we cache sigma and slope evaluation
        VecA               sigma_;
        VecA               slope_;

        // helper
        inline size_t idx(const DateType t) {
            if ((t <= times_[0]) | (times_.size()<2)) return 0;
            if (t >  times_[times_.size() - 2])       return times_.size() - 1;
            // bisection search
            size_t a = 0, b = times_.size() - 2;
            while (b - a>1) {
                size_t s = (a + b) / 2;
                if (t <= times_[s]) b = s;
                else                a = s;
            }
            return b;
        }

        // float leg
        inline ActiveType floatLeg( const Swap& s, DateType t, const VecA& x, const MatA& y) {
            ActiveType flg = 0.0;
            for (size_t k=0; k<s.floatTimes.size(); ++k) flg += s.floatWeights[k]*model_->ZeroBond(t,s.floatTimes[k],x,y);
            return flg;
        }

        // annuity
        inline ActiveType annuity( const Swap& s, DateType t, const VecA& x, const MatA& y) {
            ActiveType ann = 0.0;
            for (size_t k=0; k<s.fixedTimes.size(); ++k) ann += s.fixedWeights[k]*model_->ZeroBond(t,s.fixedTimes[k],x,y);
            return ann;
        }

        // forward swaprate
        inline ActiveType swapRate( const Swap& s, DateType t, const VecA& x, const MatA& y) {
            return floatLeg(s,t,x,y) / annuity(s,t,x,y);
        }

        // d P(t,T,x,y) / dx  = -P(t,T,x,y) G(t,T) 
        inline VecA zcbGradient( DateType t, DateType T, const VecA& x, const MatA& y) {
            VecA grad(x.size());
            ActiveType DF = model_->ZeroBond(t,T,x,y);
            for (size_t k=0; k<grad.size(); ++k) grad[k] = -DF * model_->G(k,t,T);
            return grad;
        }

        // d^2 P(t,T,x,y) / dx dx_i  = P(t,T,x,y) G(t,T) G_i(t,T) 
        inline VecA zcbHessian( size_t i, DateType t, DateType T, const VecA& x, const MatA& y) {
            VecA hess(x.size());
            ActiveType DF = model_->ZeroBond(t,T,x,y);
            for (size_t k=0; k<hess.size(); ++k) hess[k] = DF * model_->G(k,t,T) * model_->G(i,t,T);
            return hess;
        }

        // float leg gradient
        inline VecA floatLegGradient( const Swap& s, DateType t, const VecA& x, const MatA& y) {
            VecA grad(x.size(),0.0);
            for (size_t k=0; k<s.floatTimes.size(); ++k) {
                VecA zcbGrad = zcbGradient(t,s.floatTimes[k],x,y);
                for (size_t i=0; i<x.size(); ++i) grad[i] += s.floatWeights[k] * zcbGrad[i];
            }
            return grad;
        }		

        // annuity gradient
        inline VecA annuityGradient( const Swap& s, DateType t, const VecA& x, const MatA& y) {
            VecA grad(x.size(),0.0);
            for (size_t k=0; k<s.fixedTimes.size(); ++k) {
                VecA zcbGrad = zcbGradient(t,s.fixedTimes[k],x,y);
                for (size_t i=0; i<x.size(); ++i) grad[i] += s.fixedWeights[k] * zcbGrad[i];
            }
            return grad;
        }

        // float leg Hessian
        inline MatA floatLegHessian(const Swap& s, DateType t, const VecA& x, const MatA& y) {
            MatA hess(x.size(), VecA(x.size(), 0.0));
            for (size_t k = 0; k<s.floatTimes.size(); ++k) {
                ActiveType DF = model_->ZeroBond(t, s.floatTimes[k], x, y);
                for (size_t i = 0; i < x.size(); ++i)
                    for (size_t j = 0; j < x.size(); ++j)
                        hess[i][j] += s.floatWeights[k] * DF * model_->G(i, t, s.floatTimes[k]) * model_->G(j, t, s.floatTimes[k]);
            }
            return hess;
        }

        // annuity Hessian
        inline MatA annuityHessian(const Swap& s, DateType t, const VecA& x, const MatA& y) {
            MatA hess(x.size(), VecA(x.size(), 0.0));
            for (size_t k = 0; k<s.fixedTimes.size(); ++k) {
                ActiveType DF = model_->ZeroBond(t, s.fixedTimes[k], x, y);
                for (size_t i = 0; i < x.size(); ++i)
                    for (size_t j = 0; j < x.size(); ++j)
                        hess[i][j] += s.fixedWeights[k] * DF * model_->G(i, t, s.fixedTimes[k]) * model_->G(j, t, s.fixedTimes[k]);
            }
            return hess;
        }


        // gradient
        VecA swapGradient( const Swap& s, DateType t, const VecA& x, const MatA& y) {
            VecA grad = floatLegGradient(s,t,x,y);
            VecA angr = annuityGradient(s,t,x,y);
            ActiveType annuit = annuity(s,t,x,y);
            ActiveType dSdAnnuity = - swapRate(s,t,x,y) / annuit;
            for (size_t k=0; k<grad.size(); ++k) grad[k] = grad[k]/annuit + angr[k]*dSdAnnuity;
            return grad;
        }

        // hessian
        MatA swapHessian(const Swap& s, DateType t, const VecA& x, const MatA& y) {
            // we set f = u v, v = 1/z, u...floatleg, z...annuity
            // f'' = u''v + u'v' + [u'v']^T + uv''
            // v'  = -1/z/z z',  v'' = 1/z/z [2/z z'z' - z'']
            ActiveType floatL = floatLeg(s, t, x, y);       // u
            ActiveType annuit = annuity(s, t, x, y);        // z
            VecA floatGrad = floatLegGradient(s, t, x, y);  // u'
            VecA annuGrad  = annuityGradient(s, t, x, y);   // z'
            MatA floatHess = floatLegHessian(s, t, x, y);   // u''
            MatA annuHess  = annuityHessian(s, t, x, y);    // z''
            // we calculate u''v
            MatA hess(floatHess);
            for (size_t i = 0; i < x.size(); ++i)
                for (size_t j = 0; j < x.size(); ++j)
                    hess[i][j] /= annuit;
            // we add u'v' + [u'v']^T
            VecA vPrime(annuGrad);
            for (size_t i = 0; i < x.size(); ++i) vPrime[i] *= (-1.0) / annuit / annuit;
            for (size_t i = 0; i < x.size(); ++i)
                for (size_t j = 0; j < x.size(); ++j)
                    hess[i][j] += floatGrad[i] * vPrime[j] + floatGrad[j] * vPrime[i];
            // we add uv''
            for (size_t i = 0; i < x.size(); ++i)
                for (size_t j = 0; j < x.size(); ++j)
                    hess[i][j] += floatL / annuit / annuit * (2.0 / annuit * annuGrad[i] * annuGrad[j] - annuHess[i][j]);
            return hess;
        }
            

        // E^A [ x(T) ] = H(T)H(t)^-1 E^A [ x(t) ] + int_t^T H(T)H(s)^-1 [ E^A[y(s)]*1 + sigma_x^T sigma_x G_A(s) ] ds
        inline VecA expectationAx( const Swap& s, DateType t, DateType T, const MatA& barYtT, const VecA& barXt ) {
            VecA barXT = barXt;
            // H(T)H(t)^-1 E^A [ x(t) ]
            for (size_t i=0; i<model_->factors()-1; ++i) barXT[i] *= exp(-model_->chi()[i] * (T-t) );
            // int_t^T H(T)H(s)^-1 [ E^A[y(s)]*1 + sigma_x^T sigma_x G_A(s) ] ds
            // G_A((t+T)/2)
            VecA GA(model_->factors()-1);
            for (size_t i=0; i<model_->factors()-1; ++i) {
                GA[i] = 0.0;
                for (size_t j=0; j<s.fixedTimes.size(); ++j) GA[i] += s.fixedWeights[j] * model_->ZeroBond(0.0,s.fixedTimes[j],x0_,y0_) * model_->G(i,(t+T)/2.0,s.fixedTimes[j]);
                GA[i] *= 1.0 / annuity(s,0.0,x0_,y0_);
            }
            // v = sigma_x^T [sigma_x G_A(s)]  | s=(t+T)/2
            MatA sig_xT = model_->sigma_xT((t+T)/2.0,x0_,y0_);
            VecA u(model_->factors()-1);
            for (size_t i=0; i<model_->factors()-1; ++i) {
                u[i] = 0.0;
                for (size_t j=0; j<model_->factors()-1; ++j) u[i] += sig_xT[j][i]*GA[j];
            }
            VecA v(model_->factors()-1);
            for (size_t i=0; i<model_->factors()-1; ++i) {
                v[i] = 0.0;
                for (size_t j=0; j<model_->factors()-1; ++j) v[i] += sig_xT[i][j]*u[j];
                // E^A[y(s)]*1 - v = -v + E^A[y(s)]*1  | s=(t+T)/2
                v[i] *= -1.0;
                for (size_t j=0; j<model_->factors()-1; ++j) v[i] += barYtT[i][j];
                // H(T)H(s)^-1 (v - E^A[y(s)]*1)   | s=(t+T)/2
                v[i] *= exp( -model_->chi()[i] * (T-t)/2.0 );
                // mid-point rule
                barXT[i] += v[i] * (T - t);
            }
            return barXT;
        }

        // E^A [ y(T) ] = H(T) H(t)^-1 E^A [ y(t) ] H(t)^-1 H(T) + int_t^T H(T) H(s)^-1 sigma_x^T sigma_x H(s)^-1 H(T) ds
        //              = H(T) H(t)^-1 E^A [ y(t) ] H(t)^-1 H(T) + int_t^T [sigma_x^T sigma_x]_i,j exp{-(chi_i + chi_j)(t-s)} |_i,j=1,..,d
        // via mid-point rule
        inline MatA expectationAy( DateType t, DateType T, const MatA& barYt ) {
            MatA barYT = barYt;
            // H(T) H(t)^-1 E^A [ y(T) ] H(t)^-1 H(T)
            for (size_t i=0; i<model_->factors()-1; ++i) {
                for (size_t j=0; j<model_->factors()-1; ++j) {
                    barYT[i][j] *= exp(-(model_->chi()[i] + model_->chi()[j]) * (T-t) );
                }
            }
            // int_t^T H(T) H(s)^-1 sigma_x^T sigma_x H(s)^-1 H(T) ds
            MatA sig_xT = model_->sigma_xT((t+T)/2.0,x0_,y0_);
            for (size_t i=0; i<model_->factors()-1; ++i) {
                for (size_t j=0; j<model_->factors()-1; ++j) {
                    ActiveType tmp = 0.0;  // sigma^2
                    for (size_t k=0; k<model_->factors()-1; ++k) tmp += sig_xT[i][k] * sig_xT[j][k];
                    tmp *= (1.0 - exp(-(model_->chi()[i] + model_->chi()[j]) * (T - t))) / (model_->chi()[i] + model_->chi()[j]);
                    barYT[i][j] += tmp;
                    //tmp *= exp( -(model_->chi()[i] + model_->chi()[j]) * (T-t)/2.0 );
                    // mid-point rule
                    //barYT[i][j] += tmp * (T - t);
                }
            }
            return barYT;
        }

        // sigma for given x,y
        ActiveType sigma( const Swap& s, DateType t, const VecA& x, const MatA& y) {
            VecA grad  = swapGradient(s,t,x,y);
            MatA sigxT = model_->sigma_xT(t,x,y);
            VecA u(x.size());
            // grad * sigma_x^T
            for (size_t j=0; j<x.size(); ++j) {
                u[j] = 0.0;
                for (size_t i=0; i<x.size(); ++i) u[j] += grad[i] * sigxT[i][j];
            }
            // [ grad sigma_x^T sigma_x grad^T ]^{1/2}
            ActiveType num = 0.0;
            for (size_t j=0; j<x.size(); ++j) num += u[j]*u[j];
            num = sqrt(num);
            return num;
        }

        // slope for given x,y
        ActiveType slope( const Swap& s, DateType t, const VecA& x, const MatA& y) {

            // u = grad * sigma_x^T
            VecA u(x.size());
            VecA grad  = swapGradient(s,t,x,y);
            MatA sigxT = model_->sigma_xT(t,x,y);
            for (size_t j=0; j<x.size(); ++j) {
                u[j] = 0.0;
                for (size_t i=0; i<x.size(); ++i) u[j] += grad[i] * sigxT[i][j];
            }
            // den = [ grad sigma_x^T sigma_x grad^T ]^2
            ActiveType den = 0.0;
            for (size_t i=0; i<x.size(); ++i) den += u[i]*u[i];
            ActiveType sigma = sqrt(den);  // maybe better re-use sigma calculation from above...
            den = den * den;

            // d_x ...

            // v = grad c_x = grad sigma_x^T sigma_x = u * sigma_x = sigma_x^T u^T
            VecA v(x.size());
            for (size_t i=0; i<x.size(); ++i) {
                v[i] = 0.0;
                for (size_t j=0; j<x.size(); ++j) v[i] += sigxT[i][j] * u[j];
            }

            // inputs...
            MatP HHfInv = model_->HHfInv();
            MatP DfT    = model_->DfT();
            VecP delta  = model_->delta();
            VecP chi    = model_->chi();
            VecA slope_f(x.size());
            for (size_t k=0; k<x.size(); ++k) { 
                slope_f[k] = model_->slope(k,t);  // this makes it model-specific; maybe better delegate to model
            }

            // d sigma_x^T / d x_k
            MatA dx(x.size());
            MatA dsigmaxT_dxk(x.size());
            VecA diag(x.size());
            MatA tmp(x.size());
            for (size_t k=0; k<x.size(); ++k) {
                dx[k].resize(x.size(),(ActiveType)0.0);
                dsigmaxT_dxk[k].resize(x.size());
                tmp[k].resize(x.size());
            }

            // d_x = sum...
            for (size_t k=0; k<x.size(); ++k) {
                // diag^k = diag_i [ slope_{f,i} h_k(t+delta_i)/h_k(t) ]
                //        = diag_i [ slope_{f,i} exp(-chi_k delta_i)   ]
                for (size_t i=0; i<x.size(); ++i) diag[i] = slope_f[i] * exp(-chi[k]*delta[i]);
                // tmp = diag * DfT
                for (size_t i=0; i<x.size(); ++i) {
                    for (size_t j=0; j<x.size(); ++j) tmp[i][j] = diag[i] * DfT[i][j];
                }
                // d sigma_x^T / d x_k = HHf^-1 * diag * DfT = HHf^-1 * tmp
                for (size_t i=0; i<x.size(); ++i) {
                    for (size_t j=0; j<x.size(); ++j) {
                        dsigmaxT_dxk[i][j] = 0.0;
                        for (size_t l=0; l<x.size(); ++l) dsigmaxT_dxk[i][j] += HHfInv[i][l] * tmp[l][j];
                    }
                }
                // tmp = [d sigma_x^T / d x_k] sigma_x = [d sigma_x^T / d x_k] [sigma_x^T]^T
                for (size_t i=0; i<x.size(); ++i) {
                    for (size_t j=0; j<x.size(); ++j) {
                        tmp[i][j] = 0.0;
                        for (size_t l=0; l<x.size(); ++l) tmp[i][j] += dsigmaxT_dxk[i][l] * sigxT[j][l];
                    }
                }
                // d_x = sum_k v_k [tmp + tmp^T]
                for (size_t i=0; i<x.size(); ++i) {
                    for (size_t j=0; j<x.size(); ++j) dx[i][j] += v[k] * (tmp[i][j] + tmp[j][i]); 
                }

            }  // for k=0..d

            // num1 and num2...

            // num1 = grad * (dx grad^T)  (overwriting vector u)
            for (size_t i=0; i<x.size(); ++i) {
                u[i] = 0.0;
                for (size_t j=0; j<x.size(); ++j) u[i] += dx[i][j] * grad[j];
            }
            ActiveType num1 = 0.0;
            for (size_t i=0; i<x.size(); ++i)  num1 += grad[i] * u[i];

            // num2 = grad c_x grad^2 c_x grad^T = v (grad^2 v^T)
            tmp = swapHessian(s,t,x,y);
            for (size_t i=0; i<x.size(); ++i) {
                u[i] = 0.0;
                for (size_t j=0; j<x.size(); ++j) u[i] += tmp[i][j] * v[j];
            }
            ActiveType num2 = 0.0;
            for (size_t i=0; i<x.size(); ++i)  num2 += v[i] * u[i];

            // finishing up
            ActiveType res = sigma * ( num1/den/2.0 + num2/den );
            return res;
        }

    public:

        // constructor
        QGSwaprateModelT (
            const ext::shared_ptr< QuasiGaussianModel2T<DateType,PassiveType,ActiveType> >&           model,
            const VecD&                                                                               floatTimes,    // T[1], ..., T[M]
            const VecD&                                                                               floatWeights,  // u[1], ..., u[M]
            const VecD&                                                                               fixedTimes,    // T[1], ..., T[N]
            const VecD&                                                                               fixedWeights,  // w[1], ..., w[N-1]
            const VecD&                                                                               modelTimes,    // time grid for numerical integration
            const bool                                                                                useExpectedXY  // evaluate E^A [ x(t) ], E^A [ y(t) ] as expansion points
            ) 	: model_(model), swap_(floatTimes,floatWeights,fixedTimes,fixedWeights), times_(modelTimes), useExpectedXY_(useExpectedXY) {
            // check consistency of swap
            // float leg
            QL_REQUIRE(swap_.floatWeights.size()>0,"TemplateQGSwaptionModel: empty float weights.");
            QL_REQUIRE(swap_.floatTimes.size()==swap_.floatWeights.size(),"TemplateQGSwaptionModel: float sizes mismatch.");
            QL_REQUIRE(swap_.floatTimes[0]>0,"TemplateQGSwaptionModel: future float times required");
            for (size_t k=1; k<swap_.floatTimes.size(); ++k) QL_REQUIRE(swap_.floatTimes[k]>=swap_.floatTimes[k-1],"TemplateQGSwaptionModel: ascending float times required");
            // fixed leg
            QL_REQUIRE(swap_.fixedWeights.size()>0,"TemplateQGSwaptionModel: empty fixed weights.");
            QL_REQUIRE(swap_.fixedTimes.size()==swap_.fixedWeights.size(),"TemplateQGSwaptionModel: fixed sizes mismatch.");
            QL_REQUIRE(swap_.fixedTimes[0]>0,"TemplateQGSwaptionModel: future fixed times required");
            for (size_t k=1; k<swap_.fixedTimes.size(); ++k) QL_REQUIRE(swap_.fixedTimes[k]>=swap_.fixedTimes[k-1],"TemplateQGSwaptionModel: ascending fixed times required");
            // check consistency of times
            QL_REQUIRE(times_.size()>1,"TemplateQGSwaptionModel: at least two model times required");
            QL_REQUIRE(times_[0]==0,"TemplateQGSwaptionModel: model times T[0]=0 required");
            for (size_t k=1; k<times_.size(); ++k) QL_REQUIRE(times_[k]>times_[k-1],"TemplateQGSwaptionModel: ascending model times required");
            QL_REQUIRE(times_[times_.size()-1]<=swap_.floatTimes[0],"TemplateQGSwaptionModel: model times T[n] <= T_float[0] required");
            QL_REQUIRE(times_[times_.size()-1]<=swap_.fixedTimes[0],"TemplateQGSwaptionModel: model times T[n] <= T_fixed[0] required");
            // initialise x0_, y0_
            x0_ = VecA(model_->factors()-1,0.0);
            y0_ = MatA(model_->factors()-1,x0_);
            // initialise barX_, barY_
            barX_ = std::vector<VecA>(times_.size(),x0_);
            barY_ = std::vector<MatA>(times_.size(),y0_);
            // evaluate expectations...
            if (useExpectedXY_) evaluateExpXY();
            // cache sigma and slope
            sigma_.resize(times_.size());
            slope_.resize(times_.size());
            for (size_t k = 0; k < times_.size(); ++k) {
                sigma_[k] = sigma(swap_, times_[k], barX(times_[k]), barY(times_[k]));
                slope_[k] = slope(swap_, times_[k], barX(times_[k]), barY(times_[k]));
            }
        }

        virtual ~QGSwaprateModelT() = default;

        inline const ext::shared_ptr< QuasiGaussianModel2T<DateType, PassiveType, ActiveType> >&  model() { return model_; }
        inline const VecD& floatTimes()    { return swap_.floatTimes;   }
        inline const VecD& floatWeights()  { return swap_.floatWeights; }
        inline const VecD& fixedTimes()    { return swap_.fixedTimes;   }
        inline const VecD& fixedWeights()  { return swap_.fixedWeights; }
        inline const VecD& modelTimes()    { return times_;             }
        inline const bool  useExpectedXY() { return useExpectedXY_;     }



        inline void evaluateExpXY() {
            // first we need y as it inputs to x
            for (size_t k=1; k<times_.size(); ++k) {
                barY_[k] = expectationAy(times_[k-1],times_[k],barY_[k-1]);
                barX_[k] = expectationAx(swap_,times_[k-1],times_[k],barY((times_[k-1]+times_[k])/2.0),barX_[k-1]);
            }
        }

        // inspect model parameters

        // linear interpolation
        inline VecA barX( DateType t ) {
            if (t<=times_[0])               return barX_[0];
            if (t>=times_[times_.size()-1]) return barX_[barX_.size()-1];
            size_t k = idx(t);
            VecA x(barX_[k-1]);
            DateType rho = (t - times_[k-1]) / (times_[k] - times_[k-1]);
            for (size_t i=0; i<x.size(); ++i) x[i] += rho * (barX_[k][i] - barX_[k-1][i]);
            return x;
        }

        // linear interpolation
        inline MatA barY( DateType t ) {
            if (t<=times_[0])               return barY_[0];
            if (t>=times_[times_.size()-1]) return barY_[barY_.size()-1];
            size_t k = idx(t);
            MatA y(barY_[k-1]);
            DateType rho = (t - times_[k-1]) / (times_[k] - times_[k-1]);
            for (size_t i=0; i<y.size(); ++i) {
                for (size_t j=0; j<y[i].size(); ++j) {
                    y[i][j] += rho * (barY_[k][i][j] - barY_[k-1][i][j]);
                }
            }
            return y;
        }

        // the time-indexed inspectors are unsafe in the sense that they don't check if t > floatTimes[0]

        inline ActiveType annuity(DateType t, const VecA& x, const MatA& y)      { return annuity(swap_,t,x,y);      }
        inline ActiveType swapRate(DateType t, const VecA& x, const MatA& y)     { return swapRate(swap_, t, x, y);  }
        inline VecA       swapGradient(DateType t, const VecA& x, const MatA& y) { return swapGradient(swap_,t,x,y); }
        inline MatA       swapHessian(DateType t, const VecA& x, const MatA& y)  { return swapHessian(swap_,t,x,y);  }

        // inspectors 
        inline virtual ActiveType sigma(const DateType t) {
            // return sigma(swap_,t,barX(t),barY(t));
            // linear interpolation of cached values
            if (t <= times_[0])                 return sigma_[0];
            if (t >= times_[times_.size() - 1]) return sigma_[times_.size() - 1];
            size_t k = idx(t);
            ActiveType x = sigma_[k - 1];
            DateType rho = (t - times_[k - 1]) / (times_[k] - times_[k - 1]);
            x += rho * (sigma_[k] - sigma_[k - 1]);
            return x;
        }
            
        inline virtual ActiveType slope ( const DateType t ) {
            // return slope(swap_,t,barX(t),barY(t)); 
            // linear interpolation of cached values
            if (t <= times_[0])                 return slope_[0];
            if (t >= times_[times_.size() - 1]) return slope_[times_.size() - 1];
            size_t k = idx(t);
            ActiveType x = slope_[k - 1];
            DateType rho = (t - times_[k - 1]) / (times_[k] - times_[k - 1]);
            x += rho * (slope_[k] - slope_[k - 1]);
            return x;
        }

        inline virtual  ActiveType eta( const DateType t )   { return model_->eta(t);                 }

        inline virtual ActiveType  theta()                   { return model_->theta();                 }
        inline virtual ActiveType  z0()                      { return model_->z0();                    }
        inline virtual ActiveType  rho()                     { return 0.0;                             }
        inline virtual ActiveType  S0()                      { return swapRate(swap_,0.0,x0_,y0_);     }

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
        inline virtual VecA drift(const DateType t, const VecA& X) {
            VecA a(2);
            // S-variable drift-less
            a[0] = 0.0;
            // z-variable theta [ m - z(t)^+ ]  (full truncation)
            a[1] = theta()*(z0() - ((X[1]>0) ? (X[1]) : (0.0)));
            return a;
        }
        // b[t,X(t)]
        inline virtual MatA diffusion(const DateType t, const VecA& X) {
            MatA B(2);
            B[0].resize(2);
            B[1].resize(2);
            ActiveType sqrtz = ((X[1]>0) ? (sqrt(X[1])) : (0.0));   // full truncation
                                                                    // S-variable lambda(t) [ b(t) S(t) + (1-b(t)) L ] sqrt[z(t)] dW(t)
            B[0][0] = (sigma(t) + slope(t)*(X[0] - S0())) * sqrtz;
            B[0][1] = 0.0;
            // z-variable
            B[1][0] = rho()*eta(t)*sqrtz;
            B[1][1] = sqrt(1 - rho()*rho())*eta(t)*sqrtz;
            // finished
            return B;
        }

        // integrate X1 = X0 + drift()*dt + diffusion()*dW*sqrt(dt)
        inline virtual void evolve(const DateType t0, const VecA& X0, const DateType dt, const VecD& dW, VecA& X1) {
            // first we simulate stochastic vol via lognormal approximation...
            ActiveType e = model_->expectationZ(t0, X0[1], dt);
            ActiveType v = model_->varianceZ(t0, X0[1], dt);
            ActiveType dZ = dW[1];  // last risk factor is for vol process
            ActiveType si = sqrt(log(1.0 + v / e / e));
            ActiveType mu = log(e) - si*si / 2.0;
            X1[1] = exp(mu + si*dZ);
            ActiveType averageZ = 0.5*(X0[1] + X1[1]);  // we freeze z for subsequent calculation
            // S-variable
            ActiveType sqrtz = (averageZ > 0.0) ? (sqrt(averageZ)) : (0.0);  // full truncation as safe guard
            ActiveType localVol = sigma(t0) + slope(t0)*(X0[0] - S0());  // maybe average sigma and slope...
            localVol = (localVol > 0.0) ? (localVol) : (0.0);  // truncate negative local volatility
            X1[0] = X0[0] + localVol * sqrtz * dW[0] * sqrt(dt);
            return;
        }

        // the numeraire in the domestic currency used for discounting future payoffs
        inline virtual ActiveType numeraire(const DateType t, const VecA& X) { return 1.0; }

        // a domestic currency zero coupon bond
        inline virtual ActiveType zeroBond(const DateType t, const DateType T, const VecA& X) { return 1.0; }

        // default implementation for single-asset models
        inline virtual ActiveType asset(const DateType t, const VecA& X, const std::string& alias) { return X[0]; }

        // the expectation E^T in the domestic currency terminal meassure
        // this is required to calculate asset adjusters without knowing the implementation of the model
        inline virtual ActiveType forwardAsset(const DateType t, const DateType T, const VecA& X, const std::string& alias) { return asset(t, X, alias); }

    };


}

#endif  /* ifndef quantlib_templateqgswapratemodelT_hpp */
