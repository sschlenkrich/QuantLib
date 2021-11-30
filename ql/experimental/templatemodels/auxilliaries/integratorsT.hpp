/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2010 Sebastian Schlenkrich

*/

/*! \file integratorsT.hpp
    \brief provide template functions for numerical integration
*/


#ifndef quantlib_templateintegrators_hpp
#define quantlib_templateintegrators_hpp

#include <boost/function.hpp>
//#include <boost/math/special_functions/erf.hpp>
//#include <ql/experimental/template/auxilliaries/MinimADVariable2.hpp>


namespace TemplateAuxilliaries {

    // evaluate \int_a^b v(t) f(t) dt = \sum v_i [F(t_i) - F(t_i-1)] with
    // v(t) piece-wise left-constant,
    // F'(t) = f(t)
    template <typename PassiveType, typename ActiveType, typename FuncType>
    class PieceWiseConstantIntegral {
    private:
        std::vector<PassiveType> t_;
        std::vector<ActiveType> v_;
        FuncType F_;
    public:
        PieceWiseConstantIntegral(const std::vector<PassiveType>& t, const std::vector<ActiveType>& v, const FuncType& F) : t_(t), v_(v), F_(F) {}
        ActiveType operator()(PassiveType startTime, PassiveType endTime) {
            int sgn = 1;
            if (startTime>endTime) {  // we want to ensure startTime <= endTime
                PassiveType t = startTime;
                startTime = endTime;
                endTime = t;
                sgn = -1;
            }
            // organising indices
            size_t idx_min  = 0;
            size_t idx_max  = std::min(t_.size(),v_.size())-1;
            size_t idx_last = idx_max;
            // enforce a < t_min <= t_max < b or special treatment
            while ((startTime>=t_[idx_min])&&(idx_min<idx_last)) ++idx_min;
            while ((endTime  <=t_[idx_max])&&(idx_max>0       )) --idx_max;
            ActiveType tmp = sgn * ( F_(endTime) - F_(startTime) );
            if (endTime<=t_[0])    return v_[0]        * tmp;  // short end
            if (idx_min==idx_last) return v_[idx_last] * tmp;  // long end
            if (idx_min> idx_max)  return v_[idx_min]  * tmp;  // integration within grid interval
            // integral a ... x_min
            tmp = v_[idx_min] * ( F_(t_[idx_min]) - F_(startTime) );
            // integral x_min ... x_max
            for (size_t i=idx_min; i<idx_max; ++i) tmp += v_[i+1] * ( F_(t_[i+1]) - F_(t_[i]) );
            // integral x_max ... b
            if (idx_max<idx_last) tmp += v_[idx_max+1] * ( F_(endTime) - F_(t_[idx_max]) );
            else                  tmp += v_[idx_max]   * ( F_(endTime) - F_(t_[idx_max]) );
            // finished
            return sgn * tmp;
        }
    };

    // evaluate \int_x(0)^x(n) v(x) f(x) dx via trapezoidal rule
    // v(x) interpolated values on variable x-grid 
    // f(x) scalar function as functor
    class TrapezoidalIntegral {
    public:
        template <typename PassiveType, typename ActiveType, typename FuncType>
        ActiveType operator()(const std::vector<PassiveType>& x, const std::vector<ActiveType>& v, const FuncType& f) {
            size_t n=std::min(x.size(),v.size());
            if (n<2) return (ActiveType)0.0;
            ActiveType sum=0;
            for (size_t i=0; i<n-1; ++i) sum += 0.5*(v[i]*f(x[i]) + v[i+1]*f(x[i+1]))*(x[i+1] - x[i]);
            return sum;
        }
    };

    // evaluate \int_x[0]^x[n-1] v(x) f(x) dx via Gau�-Tschebyschow-Integration
    // x[0] left boundary, x[n-1] right boundary
    // x[1], ..., x[n-2] Gauss-Tschebyschow grid points
    // v(x) interpolated values on x-grid 
    // f(x) scalar function as functor
    class GaussTschebyschowIntegral {
    public:
        template <typename PassiveType>
        std::vector<PassiveType> getGrid(PassiveType a, PassiveType b, size_t n) {
            std::vector<PassiveType> x(n);
            if (n==0) return x;
            if (n==1) { x[0] = 0.5*(a+b); return x; }
            x[0] = a; 
            x[n-1] = b; 
            if (n==2) return x; 
            // n>2
            for (size_t k=1; k<n-1; ++k) {
                x[k] = -cos( (2.0*k-1.0) / (2.0*(n-2.0)) * M_PI );      // x \in (-1, 1)
                x[k] = 0.5*(b-a)*(x[k]+1) + a;                 // x \in ( a, b)
            }
            return x;
        }
        template <typename PassiveType, typename ActiveType, typename FuncType>
        ActiveType operator()(const std::vector<PassiveType>& x, const std::vector<ActiveType>& v, const FuncType& f) {
            size_t n=std::min(x.size(),v.size());
            if (n<3) return TrapezoidalIntegral()(x,v,f);
            ActiveType sum=0;
            for (size_t k=1; k<n-1; ++k) {
                sum += v[k] * f(x[k]) * sin( (2*k-1)/(2*(n-2))*M_PI );
            }
            sum *= M_PI / (n-2);
            return sum;
        }
    };


    // 4th order Runge Kutta step for y' = f(t,y)
    // via y1   = y0 + b^T k dt
    //     k_i  = f(t+c_i dt, y0 + a_i^T k dt)
    template <typename DateType, typename ActiveType>
    void rungeKuttaStep( const std::vector<ActiveType>&                                                                      y0, 
                         const DateType                                                                                      t, 
                         const boost::function< void (const DateType, const std::vector<ActiveType>&, std::vector<ActiveType>&) >& f,
                         const DateType                                                                                      dt,
                         std::vector<ActiveType>&                                                                      y1  ) {
        std::vector<ActiveType> k1(y0.size()), k2(y0.size()), k3(y0.size()), k4(y0.size());
        y1 = y0;
        // we add an epsilon in case we are at a boundary of piecewise constant parameters
        DateType eps = 1.0e-8*dt;
        // k1 = f(t,y0)
        f( t+eps,    y0, k1);
        for (size_t i=0; i<y0.size(); ++i) y1[i] += 1.0/6.0 * dt * k1[i];
        // k2 = f(t + 0.5dt, y0 + 0.5 k1)
        for (size_t i=0; i<y0.size(); ++i) k1[i] = y0[i] + 0.5 * dt * k1[i];
        f( t+0.5*dt, k1, k2);
        for (size_t i=0; i<y0.size(); ++i) y1[i] += 1.0/3.0 * dt * k2[i];
        // k3 = f(t + 0.5dt, y0 + 0.5 k2)
        for (size_t i=0; i<y0.size(); ++i) k2[i] = y0[i] + 0.5 * dt * k2[i];
        f( t+0.5*dt, k2, k3);
        for (size_t i=0; i<y0.size(); ++i) y1[i] += 1.0/3.0 * dt * k3[i];
        // k4 = f(t + dt, y0 + k3)
        for (size_t i=0; i<y0.size(); ++i) k3[i] = y0[i] + dt * k3[i];
        f( t+1.0*dt-eps, k3, k4);
        for (size_t i=0; i<y0.size(); ++i) y1[i] += 1.0/6.0 * dt * k4[i];
        return;
    }

}

#endif  /* quantlib_templateintegrators_hpp */
