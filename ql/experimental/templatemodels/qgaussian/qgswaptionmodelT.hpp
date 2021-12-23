/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2015 Sebastian Schlenkrich

*/

/*! \file qgswaptionmodelT.hpp
    \brief (Approximate) swaption pricing for multi-factor quasi-Gaussian model with stochastic vol
               
               dS(t) = lambda [(1-b)L +bS(0)] sqrt(z) dW
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


#ifndef quantlib_templateqgswaptionmodel_hpp
#define quantlib_templateqgswaptionmodel_hpp

#include <ql/shared_ptr.hpp>
#include <ql/types.hpp>

#include <ql/experimental/templatemodels/auxilliaries/auxilliariesT.hpp>
#include <ql/experimental/templatemodels/qgaussian/quasigaussianmodelT.hpp>
#include <ql/experimental/templatemodels/stochvol/hestonmodelT.hpp>
#include <ql/experimental/templatemodels/stochvol/tdstochvolmodelT.hpp>


//#include <ql/experimental/template/auxilliaries/gausslobatto.hpp>



namespace QuantLib {

    // Declaration of the quasi-Gaussian model class
    template <class DateType, class PassiveType, class ActiveType>
    class QGSwaptionModelT : public TimeDependentStochVolModelT<DateType,PassiveType,ActiveType> {
        // from base class
        using typename TimeDependentStochVolModelT<DateType, PassiveType, ActiveType>::MidPointIntegration;


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
        ext::shared_ptr< QuasiGaussianModelT<DateType,PassiveType,ActiveType> >   model_;

        // vanilla swap specs for this model
        Swap               swap_;

        // time grid for barX, barY evaluation
        VecD               times_;

        // reference to averaging algorithm
        ext::shared_ptr< MidPointIntegration >  mp_;

        std::vector<VecA>  barX_;     // E^A [ x(t) ] 
        std::vector<MatA>  barY_;     // E^A [ y(t) ] 
        VecA               x0_;       // = VecA(model_->factors()-1,0.0);
        MatA               y0_;       // = MatA(model_->factors()-1,x0_);

        // helper
        inline size_t idx( const DateType t ) { return TemplateAuxilliaries::idx(times_,t); }

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

        // gradient
        VecA swapGradient( const Swap& s, DateType t, const VecA& x, const MatA& y) {
            //VecA grad(x.size());
            //VecA gZCB(x.size());
            // numerator
            //grad = zcbGradient(t,s.times[0],x,y);
            //gZCB = zcbGradient(t,s.times[s.N],x,y);
            VecA grad = floatLegGradient(s,t,x,y);
            VecA angr = annuityGradient(s,t,x,y);
            ActiveType annuit = annuity(s,t,x,y);
            ActiveType dSdAnnuity = - swapRate(s,t,x,y) / annuit;
            for (size_t k=0; k<grad.size(); ++k) grad[k] = grad[k]/annuit + angr[k]*dSdAnnuity;
            return grad;
        }

        // hessian
        MatA swapHessian( const Swap& s, DateType t, const VecA& x, const MatA& y) {
            MatA hess(x.size());
            VecA grad(x.size(),0.0);
            VecA gZCB(x.size(),0.0);
            VecA dGrad(x.size(),0.0);
            VecA dGZCB(x.size(),0.0);
            ActiveType swpRate = swapRate(s,t,x,y);
            ActiveType annuit  = annuity(s,t,x,y);
            VecA swapGrad = swapGradient(s,t,x,y);
            VecA annuGrad = annuityGradient(s,t,x,y);
            for (size_t i=0; i<x.size(); ++i) {
                // numerator
                //grad  = zcbGradient(t,s.times[0],x,y);
                //gZCB  = zcbGradient(t,s.times[s.N],x,y);
                //dGrad = zcbHessian(i,t,s.times[0],x,y);
                //dGZCB = zcbHessian(i,t,s.times[s.N],x,y);
                for (size_t j=0; j<s.floatTimes.size(); ++j) {
                    gZCB  = zcbGradient(t,s.floatTimes[j],x,y);
                    dGZCB = zcbHessian(i,t,s.floatTimes[j],x,y);
                    for (size_t k=0; k<grad.size(); ++k) {
                        grad[k]  += s.floatWeights[j] * gZCB[k] / annuit;                                       // z  = u / v
                        dGrad[k] += s.floatWeights[j] * (dGZCB[k]*annuit - gZCB[k]*annuGrad[i])/annuit/annuit;  // z' = (u'v - uv')/v^2
                    }
                    //for (size_t k=0; k<grad.size(); ++k) {
                    //    grad[k]  = (grad[k] - gZCB[k]) / annuit;                        // z  = u / v
                    //    dGrad[k] = (dGrad[k] - dGZCB[k] - grad[k]*annuGrad[i])/annuit;  // z' = (u' - z v') / v
                    // }
                }
                // denumerator
                ActiveType dSdAnnuity  = - swpRate / annuit;
                ActiveType ddSdAnnuity = - (swapGrad[i] + dSdAnnuity*annuGrad[i])/annuit;   // check + vs -
                for (size_t j=0; j<s.fixedTimes.size(); ++j) {
                    gZCB  = zcbGradient(t,s.fixedTimes[j],x,y);
                    dGZCB = zcbHessian(i,t,s.fixedTimes[j],x,y);
                    for (size_t k=0; k<grad.size(); ++k) {
                        grad[k]  += s.fixedWeights[j] * dSdAnnuity * gZCB[k];
                        dGrad[k] += s.fixedWeights[j] * (ddSdAnnuity*gZCB[k] + dSdAnnuity*dGZCB[k]);
                    }
                }
                // allocate and copy
                hess[i] = dGrad;
            }
            return hess;
        }


        /*
        VecA swapGradient( const Swap& s, DateType t, const VecA& x, const MatA& y) {
            PassiveType eps = 1.0e-7;
            VecA grad(x.size());
            VecA xm(x), xp(x);
            for (size_t k=0; k<x.size(); ++k) {
                // bump
                xm[k] -= eps;
                xp[k] += eps;
                // recalculate
                ActiveType ym = swapRate(s,t,xm,y);
                ActiveType yp = swapRate(s,t,xp,y);
                // first derivs
                grad[k] = yp/2.0/eps - ym/2.0/eps;
                // reset
                xm[k] = x[k];
                xp[k] = x[k];
            }
            return grad;
        }
        */
        /*
        // hessian via central differences to avoid errors in tedious calculations
        MatA swapHessian( const Swap& s, DateType t, const VecA& x, const MatA& y) {
            PassiveType eps = 1.0e-7;
            MatA hess(x.size());
            VecA xm(x), xp(x), ym(x.size()), yp(x.size());
            for (size_t k=0; k<x.size(); ++k) {
                // bump
                xm[k] -= eps;
                xp[k] += eps;
                // recalculate
                ym = swapGradient(s,t,xm,y);
                yp = swapGradient(s,t,xp,y);
                // second derivs
                hess[k].resize(x.size());
                for (size_t i=0; i<x.size(); ++i) hess[k][i] = yp[i]/2.0/eps - ym[i]/2.0/eps;
                // reset
                xm[k] = x[k];
                xp[k] = x[k];
            }
            // maybe better check symetry
            return hess;
        }
        */

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
                // E^A[y(s)]*1 + v  | s=(t+T)/2
                for (size_t j=0; j<model_->factors()-1; ++j) v[i] += barYtT[i][j];
                // H(T)H(s)^-1 E^A[y(s)]*1 + v  | s=(t+T)/2
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
                    // integrand at (t+T)/s
                    ActiveType tmp = 0.0;
                    for (size_t k=0; k<model_->factors()-1; ++k) tmp += sig_xT[i][k] * sig_xT[j][k];
                    tmp *= exp( -(model_->chi()[i] + model_->chi()[j]) * (T-t)/2.0 );
                    // mid-point rule
                    barYT[i][j] += tmp * (T - t);
                }
            }
            return barYT;
        }

        // lambda_S for given x,y
        ActiveType lambda( const Swap& s, DateType t, const VecA& x, const MatA& y) {
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
            // swapt rate at (t,x,y)=0
            //VecA x0(x.size(),(ActiveType)0.0);
            //MatA y0(x.size());
            //for (size_t i=0; i<x.size(); ++i) y0[i].resize(x.size(),(ActiveType)0.0);
            ActiveType den = swapRate(s,0.0,x0_,y0_);
            return num / den;
        }

        // b_S for given x,y
        ActiveType b( const Swap& s, DateType t, const VecA& x, const MatA& y) {

            // swapt rate at (t,x,y)=0
            VecA x0(x.size(),(ActiveType)0.0);
            MatA y0(x.size());
            for (size_t i=0; i<x.size(); ++i) y0[i].resize(x.size(),(ActiveType)0.0);
            ActiveType S0 = swapRate(s,0.0,x0,y0);

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
            VecA lambda_f(x.size());
            VecA b_f(x.size());
            for (size_t k=0; k<x.size(); ++k) { 
                lambda_f[k] = model_->lambda(k,t);
                b_f[k]      = model_->b(k,t);
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
                // diag^k = diag_i [ lambda_{f,i} b_{f,i} h_k(t+delta_i)/h_k(t) ]
                //        = diag_i [ lambda_{f,i} b_{f,i} exp(-chi_k delta_i)   ]
                for (size_t i=0; i<x.size(); ++i) diag[i] = lambda_f[i] * b_f[i] * exp(-chi[k]*delta[i]);
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
            ActiveType res = S0 * ( num1/den/2.0 + num2/den );
            return res;
        }

    public:

        // constructor
        QGSwaptionModelT (
            const ext::shared_ptr< QuasiGaussianModelT<DateType,PassiveType,ActiveType> >&          model,
            const VecD&                                                                               floatTimes,    // T[1], ..., T[M]
            const VecD&                                                                               floatWeights,  // u[1], ..., u[M]
            const VecD&                                                                               fixedTimes,    // T[1], ..., T[N]
            const VecD&                                                                               fixedWeights,  // w[1], ..., w[N-1]
            const VecD&                                                                               modelTimes,    // time grid for numerical integration
            const bool                                                                                useExpectedXY  // evaluate E^A [ x(t) ], E^A [ y(t) ] as expansion points
            ) 	: model_(model), swap_(floatTimes,floatWeights,fixedTimes,fixedWeights), times_(modelTimes), mp_(new MidPointIntegration(this,modelTimes))   {
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
            QL_REQUIRE(times_[times_.size()-1]>=swap_.floatTimes[0],"TemplateQGSwaptionModel: model times T[n] >= T_float[0] required");
            QL_REQUIRE(times_[times_.size()-1]<=swap_.fixedTimes[0],"TemplateQGSwaptionModel: model times T[n] <= T_fixed[0] required");
            // initialise x0_, y0_
            x0_ = VecA(model_->factors()-1,0.0);
            y0_ = MatA(model_->factors()-1,x0_);
            // initialise barX_, barY_
            barX_ = std::vector<VecA>(times_.size(),x0_);
            barY_ = std::vector<MatA>(times_.size(),y0_);
            // evaluate expectations...
            if (useExpectedXY) evaluateExpXY();
        }

        virtual ~QGSwaptionModelT() = default;

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

        // inspectors 
        virtual ActiveType lambda( const DateType t ) { return lambda(swap_,t,barX(t),barY(t)); }
        virtual ActiveType b     ( const DateType t ) { return b(swap_,t,barX(t),barY(t));      }
        virtual ActiveType eta   ( const DateType t ) { return model_->eta(t);                   }

        virtual ActiveType  L()                       { return swapRate(swap_,0.0,x0_,y0_);     }
        virtual ActiveType  theta()                   { return model_->theta();                 }
        virtual ActiveType  m()                       { return model_->z0();                    }
        virtual ActiveType  z0()                      { return model_->z0();                    }
        virtual ActiveType  rho()                     { return 0.0;                             }
        virtual ActiveType  S0()                      { return swapRate(swap_,0.0,x0_,y0_);     }

        // averaging formula implementations
        virtual ActiveType  averageLambda( const DateType T) { return mp_->averageLambda( T ); }
        virtual ActiveType  averageB     ( const DateType T) { return mp_->averageB( T );      }
        virtual ActiveType  averageEta   ( const DateType T) { return mp_->averageEta( T );    }

        // averaging formula implementations for current swap
        virtual ActiveType  averageLambda()                  { return mp_->averageLambda( swap_.floatTimes[0] ); }
        virtual ActiveType  averageB     ()                  { return mp_->averageB(      swap_.floatTimes[0] ); }
        virtual ActiveType  averageEta   ()                  { return mp_->averageEta(    swap_.floatTimes[0] ); }

    };


}

#endif  /* ifndef quantlib_templateqgswaptionmodel_hpp */
