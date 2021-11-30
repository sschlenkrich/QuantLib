/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2010 Sebastian Schlenkrich

*/



#ifndef quantlib_templateauxilliaries_solver1d_hpp
#define quantlib_templateauxilliaries_solver1d_hpp

#include <ql/types.hpp>
#include <boost/function.hpp>


namespace TemplateAuxilliaries {

    //! Template for 1-D solution f(x) = 0, s.t. x \in [a,b] via secant method
    template <class Type>
    Type solve1d( const boost::function<Type (Type)>& f, Type xTol, Type a, Type b, size_t nTrials = 10 ) {
        Type fa = f(a);
        Type fb = f(b);
        if (fa*fb>0) {  // we need new arguments enclosing the solution
            b = (a + b);
            a = b / 4.0;
            for (size_t k=0; k<nTrials; ++k) {
                fa = f(a);
                fb = f(b);
                if (fa*fb<=0) break;
                a = a/2.0;
                b = b*2.0;
            }
        }
        QL_REQUIRE(fa*fb<=0,"Solve1d: Can't find intervall enclosing a solution");
        if (a>b) { // swap a <-> b
            Type tmp = a; a = b; b = tmp;
            tmp = fa; fa = fb; fb = tmp;
        }
        Type m = (fb - fa)/(b-a);
        Type x1 = a;
        Type y1 = fa;
        Type s = - y1/m;
        while (fabs(s)>xTol) {
            Type x0=x1, y0=y1;
            // find a new solution
            x1 = x0 + s;
            if ((x1<a)||(x1>b)) x1 = (a + b)/2.0;
            y1 = f(x1);
            m  = (y1-y0)/(x1-x0);
            s  = -y1/m;
            // update intervalls
            if (fa*y1>=0) {
                a  = x1;
                fa = y1;
            } else {
                b  = x1;
                fb = y1;
            }
        }
        return x1;
    }
    
}

#endif  /* ifndef quantlib_solve1d_hpp */
