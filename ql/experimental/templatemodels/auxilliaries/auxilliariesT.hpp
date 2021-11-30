/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2010 Sebastian Schlenkrich

*/

/*! \file auxilliariesT.hpp
    \brief provide template functions for required computations with active data types
*/


#ifndef quantlib_templateauxilliaries_hpp
#define quantlib_templateauxilliaries_hpp

#include <boost/math/special_functions/erf.hpp>
#include <ql/experimental/templatemodels/auxilliaries/minimADVariable2T.hpp>
#include <ql/experimental/templatemodels/auxilliaries/integratorsT.hpp>



namespace TemplateAuxilliaries {

    inline double erf(const double x)     { return boost::math::erf(x); }
    inline double erf_inv(const double x) { return boost::math::erf_inv(x); }
    inline double DBL(const double x)     { return x; }
//	inline double DBL(const ADTAGEO::daglad x) { return x.val(); }
    inline double DBL(const MinimAD::Variable<QuantLib::Real> x) { return x.value(); }


    // transformation (-inf, +inf) -> (a, b)
    template <typename Type> inline Type direct(const Type x, const Type a, const Type b) {
        Type y = (x<0.0) ? (-1.0/(x-1.0)) : (-1.0/(x+1.0)+2.0);
        y = 0.5*(b-a)*y + a;
        return y;
    }

    // transformation (a, b) -> (-inf, +inf)
    template <typename Type> inline Type inverse(const Type y, const Type a, const Type b) {
        Type x = 2.0*(y-a)/(b-a);
        x = (x<1.0) ? (-1.0/x+1.0) : (-1.0/(x-2.0)-1.0);
        return x;
    }


    // square of vector elements
    template <typename Type> inline
    std::vector<Type> sqr(std::vector<Type> x) {
        std::vector<Type> x2(x.size());
        for (size_t i=0; i<x.size(); ++i) x2[i] = x[i]*x[i];
        return x2;
    }


    // find index in ascending vector
    // evaluate n s.t. t[n-1] < t <= t[n]
    template <typename Type> inline
    size_t idx(const std::vector<Type>& times, const Type t) {
        if ((t <= times[0]) | (times.size()<2)) return 0;
        if (t >  times[times.size()-2 ])        return times.size()-1;
        // bisection search
        size_t a = 0, b = times.size()-2;
        while (b-a>1) {
            size_t s = (a + b) / 2;
            if (t <= times[s]) b = s;
            else                a = s;
        }
        return b;
    }

    // union of two vectors
    template <typename Type> inline
    std::vector<Type> unionVector(const std::vector<Type>& v1, const std::vector<Type>& v2) {
        std::vector<Type> res;
        std::set_union(v1.begin(), v1.end(), v2.begin(), v2.end(), std::back_inserter(res));
        return res;
    }

    //  normal distribution and Black76 functions

    template <typename Type> inline
    Type Phi( Type x) {
        return 0.5*(erf(x*M_SQRT1_2)+1.0);
    }

    template <typename Type> inline
    Type phi( Type x) {
        return M_SQRT1_2 * M_1_SQRTPI * exp(-x*x/2.0);
    }

    template <typename Type> inline
    Type PhiInv(const Type x) {
        return M_SQRT2*erf_inv(2*x-1);
    }

    template <typename Type> inline
    Type Black76(Type F,        // forward price
                 Type K,        // strike
                 Type sigma,    // Black76 volatility
                 Type T,        // Time to expiry
                 int  cop       // call (1) or put (-1) option
                ) {
        Type d1 = log(F/K)/sigma/sqrt(T) + sigma*sqrt(T)/2.0;
        Type d2 = d1 - sigma*sqrt(T);
        return cop * (F*Phi(cop*d1) - K*Phi(cop*d2));
    }

    template <typename Type> inline
    Type Black76Vega(Type F,        // forward price
                    Type K,        // strike
                    Type sigma,    // Black76 volatility
                    Type T,        // Time to expiry
                    int    cop     // call (1) or put (-1) option
                    ) {
        Type d1 = log(F/K)/sigma/sqrt(T) + sigma*sqrt(T)/2.0;
        return F * exp(-d1*d1/2.0) * sqrt( T / 2.0 / M_PI);
    }

    template <typename Type> inline
    Type Bachelier(Type F,        // forward price
                   Type K,        // strike
                   Type sigma,    // Normal volatility
                   Type T,        // Time to expiry
                   int  cop       // call (1) or put (-1) option
                ) {
        Type d = cop*(F-K);
        Type h = d/sigma/sqrt(T);
        return d*Phi(h) + phi(h)*sigma*sqrt(T);
    }

    //  cubic interpolation and integration of expectation
    
    //  solve A X = Y by LU decomposition where A = tridiag { a, b, c }
    template <typename DateType, typename ValueType> inline void
    solveTridiagLinearSystem ( std::vector<DateType>&   a,  // input: a_1 to a_{dim-1}, output: l_1 to l_{dim-1} of L
                               std::vector<DateType>&   b,  // input: b_0 to b_{dim-1}, output: u_0 to u_{dim-1} of U
                               std::vector<DateType>&   c,  // input: c_0 to c_{dim-2}, output: v_0 to v_{dim-2} of U
                               std::vector<ValueType>&  y,  // input: right hand sides y, output: solutions x
                               std::vector<ValueType>&  z   // intermediates
                               ) {
        size_t dim = b.size();
        // in place LU decomposition; no error handling if LU decomposition does not exist
        for (size_t i=1; i<dim; ++i) {
            a[i] /= b[i-1];
            b[i] -= c[i-1]*a[i];
        }
        // forward substitution
        z[0] = y[0];
        for (size_t i=1; i<dim; ++i) z[i] = y[i] - a[i]*z[i-1];
        // backward substitution, eliminate input
        y[dim-1] = z[dim-1]/b[dim-1];
        //for (long i=dim-2; i>=0; --i) y[i] = (z[i] - c[i]*y[i+1])/b[i];
        for (size_t i=dim-1; i>0; --i) y[i-1] = (z[i-1] - c[i-1]*y[i])/b[i-1];
    }

    //  (Log) Cubic interpolation requires precomputed derivatives g[i] = y'[i]
    template <typename DateType, typename ValueType> inline
    void c2splineDerivatives( std::vector<DateType>&  x,  // input parameter, strictly increasing grid point 
                              std::vector<ValueType>& y,  // input parameter, function values at x[i] 
                              std::vector<ValueType>&       g,  // output parameter, derivatives at x[i] 
                              std::vector<ValueType>&       z,  // intermediates
                              int       logInterpolation  = 0,  // (1) interpolate log(y[i]), (0) standard   
                              int       boundaryCondition = 1   // (0) g[0] = g[dim-1] = 0,
                                                                // (1) g[0], g[dim-1] by secants   
                        ) {
        ValueType dy1, dy2;
        size_t i, dim = std::min(x.size(),y.size()); 
        std::vector<DateType> a(dim), b(dim), c(dim);
        // constant extrapolation
        if (dim==1) {
            g[0] = 0.0;
        }
        // see Wikipedia 'spline interpolation'
        for (i=1; i<dim-1; ++i) {
            a[i] = x[i+1] - x[i];
            b[i] = 2.0*(x[i+1] - x[i-1]);
            c[i] = x[i] - x[i-1];
            dy1  = (logInterpolation) ? log(y[i]/y[i-1]) : y[i] - y[i-1];
            dy2  = (logInterpolation) ? log(y[i+1]/y[i]) : y[i+1] - y[i];
            g[i] = 3.0*(a[i]/c[i]*dy1 + c[i]/a[i]*dy2);
            dy1 = 0.0; dy2 = 0.0; // elimination
        }
        // boundary conditions
        a[0] = 0.0; c[0] = 0.0; g[0] = 0.0; a[dim-1] = 0.0; c[dim-1] = 0.0; g[dim-1] = 0.0;
        b[0] = 1.0; b[dim-1] = 1.0;
        if (boundaryCondition && (dim>1)) { // fix boundaries
            g[0] =  (logInterpolation) ? log(y[1]/y[0]) : (y[1]-y[0]);
            g[0] /= (x[1]-x[0]);
            g[dim-1] =  (logInterpolation) ? log(y[dim-1]/y[dim-2]) : (y[dim-1]-y[dim-2]);
            g[dim-1] /= (x[dim-1]-x[dim-2]);
        }
        // solve tridiag [ a, b, c ] x = g and x -> g
        solveTridiagLinearSystem( a, b, c, g, z );
    }

    template <typename DateType, typename ValueType> inline
    ValueType interpolCSpline ( DateType                      xi,
                                std::vector<DateType>&  x,
                                std::vector<ValueType>& y,
                                std::vector<ValueType>& g,
                                int    logInterpolation = 0) {
        size_t prevIdx, nextIdx, idx;
        size_t dim = std::min(x.size(),y.size());
        dim = std::min(dim,g.size());
        ValueType h, u, v, p, q, yPrev, yNext, res;
        // linear extrapolation
        if (xi<=x[0]) {
            res = (logInterpolation) ? log(y[0]) : y[0];
            res += g[0]*(xi-x[0]);
            res = (logInterpolation) ? exp(res) : res;
            return res;
        }
        if (xi>=x[dim-1]) {
            res = (logInterpolation) ? log(y[dim-1]) : y[dim-1];
            res += g[dim-1]*(xi-x[dim-1]);
            res = (logInterpolation) ? exp(res) : res;
            return res;
        }
        // find prevIdx and nextIdx s.t. x[prevIdx] < xi <= x[nextIdx]
        prevIdx = 0;
        nextIdx = dim-1;
        while (nextIdx>prevIdx+1) {
            idx = (prevIdx + nextIdx)/2;
            if (xi<=x[idx]) nextIdx = idx;
            else            prevIdx = idx;
        }
        // auxilliary variables
        h = x[nextIdx] - x[prevIdx];
        u = (xi - x[prevIdx])/h;
        v = 1.0 - u;
        yPrev = (logInterpolation) ? log(y[prevIdx]) : y[prevIdx];
        yNext = (logInterpolation) ? log(y[nextIdx]) : y[nextIdx];
        p = 3.0*yNext - h*g[nextIdx];
        q = 3.0*yPrev + h*g[prevIdx];
        // final interpolation
        res = yNext*u*u*u + p*u*u*v + q*u*v*v + yPrev*v*v*v;
        res = (logInterpolation) ? exp(res) : res;
        return res;
    }


    // integration x_0, ..., x_N plus extrapolation via Bachelier formula
    template <typename PassiveType, typename ActiveType> inline
    ActiveType normalExpectation(
                std::vector<PassiveType>&   x,    //  grid points of payoff
                std::vector<ActiveType>&    v,    //  payoff
                std::vector<ActiveType>&    g,  //  derivatives of payoff for interpolation
                ActiveType                  mu,   //  expectation of normal distribution
                ActiveType                  var,  //  variance sigma^2 of normal distribution
                std::string                 method="" //  integration method  
                  ) {
        ActiveType res=0.0;
        if (x.size()==0) return res;

        // low rate extrapolation
        ActiveType Q0 = Phi((x[0]-mu)/sqrt(var));
        res = v[0]*Q0;
        if (x.size()==1) return res;  
        res -= (v[1]-v[0])/(x[1]-x[0])*Bachelier(mu,(ActiveType)x[0],sqrt(var),(ActiveType)(1.0),-1);
        //res -= g[0]*Bachelier(mu,(ActiveType)x[0],sqrt(var),(ActiveType)(1.0),-1);

        // high rate extrapolation
        ActiveType QN = Phi((x[x.size()-1]-mu)/sqrt(var));
        res += v[v.size()-1]*(1.0-QN);
        res += (v[v.size()-1]-v[v.size()-2])/(x[x.size()-1]-x[x.size()-2])*Bachelier(mu,(ActiveType)x[x.size()-1],sqrt(var),(ActiveType)(1.0),+1);
        //res += g[x.size()-1]*Bachelier(mu,(ActiveType)x[x.size()-1],sqrt(var),(ActiveType)(1.0),+1);

        // switch methods...

        // default intervall integrations
        // boundaries
        res += v[v.size()-1]*QN - v[0]*Q0;
        // replication via Put-Spreads
        ActiveType B1, B2 = Bachelier(mu,(ActiveType)x[0],sqrt(var),(ActiveType)(1.0),-1);
        for (size_t i=0; i<x.size()-1; ++i) {
            B1  = B2;
            B2  = Bachelier(mu,(ActiveType)x[i+1],sqrt(var),(ActiveType)(1.0),-1);
            res -= (v[i+1] - v[i]) / (x[i+1] - x[i]) * (B2 - B1);
            //res -= 0.5 * (g[i] + g[i+1]) * (B2 - B1);
        }
        return res;
    }


    template <typename PassiveType, typename ActiveType> inline
    ActiveType normalExpectation(
                std::vector<PassiveType>&   x,  //  grid points of payoff
                std::vector<ActiveType>&    v,  //  payoff
                std::vector<ActiveType>&    g,  //  derivatives of payoff for interpolation
                ActiveType                       mu,  //  expectation of normal distribution
                ActiveType                      var,  //  variance sigma^2 of normal distribution
                PassiveType                     tol   //  tolerance for accuracy
                  ) {

        if ((std::min(x.size(),v.size())<2)||(x.size()!=v.size())) return 0;
        if ((g.size()<2)||(g.size()!=x.size())) return 0;
        if (tol<=0) return normalExpectation( x, v, g, mu, var );
        PassiveType x0, x1, x2, h, err, tmp, lambda, sum1;
        ActiveType  v0, v1, v2;
        ActiveType  y0, y1, y2, sum2, res;
        std::vector<ActiveType> sums(x.size());
        sums.push_back((ActiveType)0.0);
        for (int i=0; i<2; ++i) {
            x0 = DBL(mu);
            v0 = interpolCSpline(x0,x,v,g);
            y0 = v0 * exp(-(x0-mu)*(x0-mu)/2.0/var);
            //h = (1-2*i)*tol; // first stupid guess
            h = (1-2*i)*x[x.size()-1]/x.size()*2;
            x2 = x0 + h;
            v2 = interpolCSpline(x2,x,v,g);
            y2 = v2 * exp(-(x2-mu)*(x2-mu)/2.0/var);				      
            // integrate mu to +infty
            tmp = Phi(DBL((x0-mu)/sqrt(var)));
            tmp = (tmp<0.5) ? tmp : (1.0-tmp);
            while ((tmp>tol*tol)||(fabs(DBL(v0))*tmp>tol*tol)) {
                x1 = x0 + 0.5*h;
                v1 = interpolCSpline(x1,x,v,g);
                y1 = v1 * exp(-(x1-mu)*(x1-mu)/2.0/var);				      
                sum1 = DBL(y0 + y2)*h/2.0;           // order-3
                sum2 = (y0 + 4.0*y1 + y2)*h/6.0;  // order-5
                err = fabs((sum1-DBL(sum2))/h);
                if (err>tol) {  // reject the step
                    h *= 0.5;
                    x2 = x1;
                    v2 = v1;
                    y2 = y1;
                    continue;  // try again with half the step size
                }
                //res += (1-2*i)*sum2;  // use the more accurate estimate
                sums.push_back(sums.back() + (1-2*i)*sum2);
                x0 = x2;
                v0 = v2;
                y0 = y2;	
                lambda = sqrt(tol/(err+1.0e-32));
                lambda = (lambda>2.0) ? 2.0 : lambda; // cap step size increase
                lambda = ((lambda>1.0)&(lambda<1.3)) ? 1.0 : lambda; // avoid step size oscillation
                h *= lambda;
                x2 = x0 + h;
                v2 = interpolCSpline(x2,x,v,g);
                y2 = v2 * exp(-(x2-mu)*(x2-mu)/2.0/var);				      
                tmp = Phi(DBL((x0-mu)/sqrt(var)));
                tmp = (tmp<0.5) ? tmp : (1.0-tmp);
            }
        }
        // don't forget the scaling factor...
        res = sums.back() / sqrt(2.0 * M_PI * var);
        // reverse elimination
        for (size_t i=sums.size(); i>0; --i) sums[i-1] = 0.0;
        return res;
    }

}

#endif  /* quantlib_templateauxilliaries_hpp */
