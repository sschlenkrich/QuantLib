/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2010 Sebastian Schlenkrich

*/



#ifndef quantlib_gausslobatto_hpp
#define quantlib_gausslobatto_hpp

#include <ql/types.hpp>
#include <boost/function.hpp>


namespace TemplateAuxilliaries {

    //! Template for Gauss-Lobatto integration, see gausslobattointegral.hpp/.cpp
    template <class Type>
    class GaussLobatto {
    private:
        //! local types
        typedef double Real;
        //! attributes
        mutable size_t maxEvaluations_;
        mutable Type absAccuracy_;
        mutable Type relAccuracy_;
        mutable Type absError_;
        mutable size_t evaluations_;
        const bool useConvergenceEstimate_;
        //! constants
        const static double alpha_, beta_, x1_, x2_, x3_;
        //! specific integration routines
        Type calculateAbsTolerance(const boost::function<Type (Type)>& f, Type a, Type b) const;
        Type adaptivGaussLobattoStep(const boost::function<Type (Type)>& f, Type a, Type b, Type fa, Type fb, Type acc) const;
    public:
        //! constructor
        GaussLobatto( size_t maxEvaluations,
                      Type absAccuracy,
                      Type relAccuracy = 0,
                      bool useConvergenceEstimate = true )
            : maxEvaluations_(maxEvaluations), absAccuracy_(absAccuracy), relAccuracy_(relAccuracy),
              absError_(0), evaluations_(0), useConvergenceEstimate_(useConvergenceEstimate) {}
        //! inspectors
        Type absoluteAccuracy() const { return absAccuracy_;      }
        Type absoluteError()    const { return absError_;         }
        size_t maxEvaluations()   const { return maxEvaluations_;   }
        size_t evaluations()      const { return evaluations_;      }
        //! modifiers
        void setNumberOfEvaluations(size_t number)         const { evaluations_ = number; }
        void increaseNumberOfEvaluations(size_t increase)  const { evaluations_ += increase; }
        //! integrate function interface
        Type integrate(const boost::function<Type (Type)>& f, Type a, Type b) const {
            //setNumberOfEvaluations(0);
            const Type calcAbsTolerance = calculateAbsTolerance(f, a, b);
            //increaseNumberOfEvaluations(2);
            return adaptivGaussLobattoStep(f, a, b, f(a), f(b), calcAbsTolerance);
        }
    };

    template <class Type> const double GaussLobatto<Type>::alpha_ = std::sqrt(2.0/3.0); 
    template <class Type> const double GaussLobatto<Type>::beta_  = 1.0/std::sqrt(5.0);
    template <class Type> const double GaussLobatto<Type>::x1_    = 0.94288241569547971906; 
    template <class Type> const double GaussLobatto<Type>::x2_    = 0.64185334234578130578;
    template <class Type> const double GaussLobatto<Type>::x3_    = 0.23638319966214988028;

    template <class Type> Type
    GaussLobatto<Type>::calculateAbsTolerance( const boost::function<Type (Type)>& f, Type a, Type b) const {
        Type relTol = (relAccuracy_ > QL_EPSILON) ? relAccuracy_ : QL_EPSILON;
        
        const Type m = (a+b)/2; 
        const Type h = (b-a)/2;
        const Type y1 = f(a);
        const Type y3 = f(m-alpha_*h);
        const Type y5 = f(m-beta_*h);
        const Type y7 = f(m);
        const Type y9 = f(m+beta_*h);
        const Type y11= f(m+alpha_*h);
        const Type y13= f(b);

        Type acc=h*(0.0158271919734801831*(y1+y13)
                  +0.0942738402188500455*(f(m-x1_*h)+f(m+x1_*h))
                  +0.1550719873365853963*(y3+y11)
                  +0.1888215739601824544*(f(m-x2_*h)+ f(m+x2_*h))
                  +0.1997734052268585268*(y5+y9) 
                  +0.2249264653333395270*(f(m-x3_*h)+f(m+x3_*h))
                  +0.2426110719014077338*y7);  
        
        increaseNumberOfEvaluations(13);
        QL_REQUIRE(acc != 0.0, "can not calculate absolute accuracy from "
                               "relative accuracy");

        Type r = 1.0;
        if (useConvergenceEstimate_) {
            const Type integral2 = (h/6)*(y1+y13+5*(y5+y9));
            const Type integral1 = (h/1470)*(77*(y1+y13)+432*(y3+y11)+
                                             625*(y5+y9)+672*y7);
        
            if (fabs(integral2-acc) != 0.0) 
                r = fabs(integral1-acc)/fabs(integral2-acc);
            if (r == 0.0 || r > 1.0)
                r = 1.0;
        }

        if (relAccuracy_ != 0)
            return  ((absoluteAccuracy() < acc*relTol) ? absoluteAccuracy() : acc*relTol)/(r*QL_EPSILON);
                    //min(absoluteAccuracy(), acc*relTol)/(r*QL_EPSILON);
        else {
            return absoluteAccuracy()/(r*QL_EPSILON);
        }
    }

    template <class Type> Type
    GaussLobatto<Type>::adaptivGaussLobattoStep(const boost::function<Type (Type)>& f, Type a, Type b, Type fa, Type fb, Type acc) const {
        QL_REQUIRE(evaluations() < maxEvaluations(),
                   "max number of iterations reached");
        
        const Type h=(b-a)/2; 
        const Type m=(a+b)/2;
        
        const Type mll=m-alpha_*h; 
        const Type ml =m-beta_*h; 
        const Type mr =m+beta_*h; 
        const Type mrr=m+alpha_*h;
        
        const Type fmll= f(mll);
        const Type fml = f(ml);
        const Type fm  = f(m);
        const Type fmr = f(mr);
        const Type fmrr= f(mrr);
        increaseNumberOfEvaluations(5);
        
        const Type integral2=(h/6)*(fa+fb+5*(fml+fmr));
        const Type integral1=(h/1470)*(77*(fa+fb)
                                       +432*(fmll+fmrr)+625*(fml+fmr)+672*fm);
        
        // avoid 80 bit logic on x86 cpu
        Type dist = acc + (integral1-integral2);
        if(dist==acc || mll<=a || b<=mrr) {
            QL_REQUIRE(m>a && b>m,"Interval contains no more machine number");
            return integral1;
        }
        else {
            return  adaptivGaussLobattoStep(f,a,mll,fa,fmll,acc)  
                  + adaptivGaussLobattoStep(f,mll,ml,fmll,fml,acc)
                  + adaptivGaussLobattoStep(f,ml,m,fml,fm,acc)
                  + adaptivGaussLobattoStep(f,m,mr,fm,fmr,acc)
                  + adaptivGaussLobattoStep(f,mr,mrr,fmr,fmrr,acc)
                  + adaptivGaussLobattoStep(f,mrr,b,fmrr,fb,acc);
        }
    }
    
}

#endif  /* ifndef quantlib_gausslobatto_hpp */
