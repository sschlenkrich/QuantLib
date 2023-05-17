/*!
     \brief Complex a template implementation for complex arithmetics
        Copyright (C) 2012 Sebastian Schlenkrich
        Version 1.0
*/

#ifndef Cpx_Complex_hpp
#define Cpx_Complex_hpp


//#define _USE_MATH_DEFINES // for Visual Studio

#include <iostream>
#include <iomanip>

#include <map>
#include <vector>
#include <iterator>
//#include <cmath>
#include <complex>

//#define DEBUG              // low-level io debugging
//#define EXECUTABLE               // compile example and main function
#define MS_VISUAL_STUDIO     // math functions not in global scope


#ifdef MS_VISUAL_STUDIO
    #include <boost/math/special_functions/erf.hpp>
#endif


namespace Cpx {

#ifdef DEBUG
    #define MSG(msg_stream) std::cerr << msg_stream << std::endl ;
#else
    #define MSG(msg_stream)
#endif

#ifdef MS_VISUAL_STUDIO
    //! Some compilers don't find intrinsic double functions
    inline double exp(double x) { return std::exp(x); }
    inline double log(double x) { return std::log(x); }
    inline double sqrt(double x) { return std::sqrt(x); }
    inline double sin(double x) { return std::sin(x); }
    inline double cos(double x) { return std::cos(x); }
    inline double tan(double x) { return std::tan(x); }
    inline double atan2(double y, double x) { return std::atan2(y,x); }
    inline double erf(double x)     { return boost::math::erf(x); }
    inline double erf_inv(double x) { return boost::math::erf_inv(x); }
#endif

    //! Variable defines the class template for derivative evaluations
    template <class Type = double>
    class Complex {
    protected:
        Type real_;
        Type imag_;
    public:
        //! "standard" constructor
        Complex( const Type& real, const Type& imag) : real_(real), imag_(imag) {}
        //! Type conversion constructor
        Complex( const Type& real) : real_(real), imag_(0) {}
        //! no standard constructor and default copy constructor
        //! inspectors
        const Type& real()  const { return real_;  }
        const Type& imag() const { return imag_; }
        //! IO steaming
        //! Apply IO to the real and imaginary component
        inline friend std::ostream& operator << (std::ostream &output, const Complex<Type> &x) {
            return output << "(" << x.real_ << "," << x.imag_ << ")";
        }
        //! Intrinsic functions
        inline friend Complex<Type> exp(const Complex<Type> &x) {
            Type r = exp(x.real_);
            return Complex<Type>( r*cos(x.imag_), r*sin(x.imag_) );
        }
        inline friend Complex<Type> log(const Complex<Type> &x) {
            return Complex<Type>( 0.5*log(x.real_*x.real_ + x.imag_*x.imag_), atan2(x.imag_,x.real_) );
        }
        inline friend Complex<Type> sqrt(const Complex<Type> &x) {
            if (x.imag_==0) {
                if (x.real_<0) return Complex<Type>( 0, sqrt(-x.real_) );
                else           return Complex<Type>( sqrt(x.real_), 0 );
            }
            if (x.real_==0) {
                if (x.imag_<0) {
                    Type tmp  = sqrt(-x.imag_)*M_SQRT1_2;
                    return Complex<Type>( tmp, -tmp );
                } else {
                    Type tmp  = sqrt(x.imag_)*M_SQRT1_2;
                    return Complex<Type>( tmp, tmp );
                }
            }
            Type sqr_r  = sqrt(sqrt(x.real_*x.real_ + x.imag_*x.imag_));
            Type phi_2  = atan2(x.imag_,x.real_) / 2;
            return Complex<Type>( sqr_r*cos(phi_2), sqr_r*sin(phi_2) );
        }

        //! Unary operators
        inline friend Complex<Type> operator + (const Complex<Type> &x ) {
            return x;
        }
        inline friend Complex<Type> operator - (const Complex<Type> &x ) {
            return Complex<Type>( -x.real_, -x.imag_ );
        }
        //! Binary operators
        //! Complex x Complex
        inline friend Complex<Type> operator + (const Complex<Type> &x, const Complex<Type> &y ) {
            return Complex<Type>( x.real_ + y.real_, x.imag_ + y.imag_ );
        }
        inline friend Complex<Type> operator - (const Complex<Type> &x, const Complex<Type> &y ) {
            return Complex<Type>( x.real_ - y.real_, x.imag_ - y.imag_ );
        }
        inline friend Complex<Type> operator * (const Complex<Type> &x, const Complex<Type> &y ) {
            return Complex<Type>( x.real_*y.real_ - x.imag_*y.imag_, x.real_*y.imag_ + x.imag_*y.real_ );
        }
        inline friend Complex<Type> operator / (const Complex<Type> &x, const Complex<Type> &y ) {
            Type den = y.real_*y.real_ + y.imag_*y.imag_;
            return Complex<Type>( (x.real_*y.real_ + x.imag_*y.imag_)/den,
                          (x.imag_*y.real_ - x.real_*y.imag_)/den );
        }
        //! Complex x Type
        inline friend Complex<Type> operator + (const Complex<Type> &x, const Type &y ) {
            return Complex<Type>( x.real_ + y, x.imag_ );
        }
        inline friend Complex<Type> operator - (const Complex<Type> &x, const Type &y ) {
            return Complex<Type>( x.real_ - y, x.imag_ );
        }
        inline friend Complex<Type> operator * (const Complex<Type> &x, const Type &y ) {
            return Complex<Type>( x.real_*y, x.imag_*y );
        }
        inline friend Complex<Type> operator / (const Complex<Type> &x, const Type &y ) {
            return Complex<Type>( x.real_/y, x.imag_/y );
        }
        //! Type x Complex
        inline friend Complex<Type> operator + (const Type &x, const Complex<Type> &y ) {
            return Complex<Type>( x + y.real_, y.imag_ );
        }
        inline friend Complex<Type> operator - (const Type &x, const Complex<Type> &y ) {
            return Complex<Type>( x - y.real_, -y.imag_ );
        }
        inline friend Complex<Type> operator * (const Type &x, const Complex<Type> &y ) {
            return Complex<Type>( x*y.real_, x*y.imag_ );
        }
        inline friend Complex<Type> operator / (const Type &x, const Complex<Type> &y ) {
            Type den = y.real_*y.real_ + y.imag_*y.imag_;
            return Complex<Type>( x*y.real_/den, -x*y.imag_/den );
        }
    };
    
    //! Type x Complex
    template <class Type>
    inline Complex<Type> operator + (const double &x, const Complex<Type> &y ) {
        return Complex<Type>( x + y.real(), y.imag() );
    }
    template <class Type>
    inline Complex<Type> operator - (const double &x, const Complex<Type> &y ) {
        return Complex<Type>( x - y.real(), -y.imag() );
    }
    template <class Type>
    inline Complex<Type> operator * (const double &x, const Complex<Type> &y ) {
        return Complex<Type>( x*y.real(), x*y.imag() );
    }
    template <class Type>
    inline Complex<Type> operator / (const double &x, const Complex<Type> &y ) {
        Type den = y.real()*y.real() + y.imag()*y.imag();
        return Complex<Type>( x*y.real()/den, -x*y.imag()/den );
    }

    
}   // namespace Cpx

#endif   // Cpx_Complex_hpp

#ifdef EXECUTABLE

int main() {
    typedef Cpx::Complex<double> complex;
    //typedef std::complex<double> complex;
    
    std::cout << "Cpx Complex (C) Sebastian Schlenktich (2012)"<< std::endl;
    std::cout << std::endl;
    
    complex a(0,4), b(1,2);
    complex c(0); // = a/b;
    c = a + complex(5.0);
    std::cout << "c = " << (c = sqrt(a)) << std::endl;
    double z;
    std::cout << (z = atan2(-1,-1)) << std::endl;
    
    return 0;
}

#endif   // EXECUTABLE

/*  Example output


*/
