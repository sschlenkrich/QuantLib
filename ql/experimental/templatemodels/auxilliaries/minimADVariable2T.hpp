/*!
\brief MinimAD a minimum purpose AD tool
Copyright (C) 2012 Sebastian Schlenkrich
Version 2.0

MinimAD provides an Algorithmic Differentiation framework
based on C++ Template Meta Programming. It is inspired by
the instant elimination approach proposed by J. Riehme
and A. Griewank (2009).

The key feature of MinimAD is the definition of a class
template Variable with a generic type (ValueType) for
its value component. All function and derivative
evaluations are defined in the template without
specifying the concrete type.

If the class template Variable is instantiated, for
example, with double then the compiler generates the
concrete source code. By inlining the generated code
further optimizations can be performed by the compiler.
Thus, the template based AD approach is not pure
operator overloading but also a kind of source
transformation by the compiler.

An important consequence of the template definition is
the fact that the class template can be instantiated with
any object type that provides arithmetic operations. Thus,
ValueType can be an active AD enabled type or, in
particular an instantiation of the template class Variable
itself. For example, by declaring objects of type
Variable< Variable<double> > second order derivatives can
easily be evaluated.

The approach for nested template instanstantiations can
be extended for higher order derivatives if required. A
further applications could, for example, be the use of
intervall types as underlying ValueType.

Version 2.0 improves performance by incorporating static
arrays for edges in addition to std::map<>. The first
edges are realized by static arrays. This avoids STL
containers in particular for temporary variables.

*/

#ifndef MinimAD_MinimADVariable_hpp
#define MinimAD_MinimADVariable_hpp


//#define _USE_MATH_DEFINES // for Visual Studio

#include <iostream>
#include <iomanip>

#include <map>
#include <vector>
#include <iterator>

//#include <cmath>


//#define DEBUG              // low-level io debugging
//#define EXECUTABLE           // compile example and main function
#define MS_VISUAL_STUDIO     // math functions not in global scope


#ifdef MS_VISUAL_STUDIO
#include <boost/math/special_functions/erf.hpp>
#endif


namespace MinimAD {

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

    //! A generic exception that can be thrown if something goes wrong
    class MinimADException : public std::exception {
    protected:
        std::string what_;
    public:
        MinimADException(const std::string& what) { what_ = what; }
        virtual const char* what() const throw()  { return what_.c_str(); }
        ~MinimADException() throw() {}
    };

    //! Control if static/dynamic containers are used
#define STATIC_REFERENCES
#define DYNAMIC_REFERENCES
    //! We abstract the template and Variable type to ease readability
#define TEMPLATE_CLASS_VARIABLE  template <class ValueType = double, class DerivativeType = ValueType, size_t N_ARG = 3, size_t N_RES = 3>
#define TEMPLATE_FUNCT_VARIABLE  template <class ValueType, class DerivativeType, size_t N_ARG, size_t N_RES>
#define VARIABLE_TYPE            Variable<ValueType,DerivativeType,N_ARG,N_RES>

    //! Variable defines the class template for derivative evaluations
    TEMPLATE_CLASS_VARIABLE
    class Variable {	
        //! ContainerType is the unique definition of the pointer maps used
        typedef std::map< VARIABLE_TYPE* ,DerivativeType* > ContainerType;
    protected:
#ifdef STATIC_REFERENCES
        VARIABLE_TYPE  *argRef_[N_ARG];
        VARIABLE_TYPE  *resRef_[N_RES];
        DerivativeType *argDrv_[N_ARG];
        DerivativeType *resDrv_[N_RES];
#endif
#ifdef DYNAMIC_REFERENCES
        ContainerType *arguments_;                // dynamic edges to preceeding/succeeding variables
        ContainerType *results_;
#endif    
        // attributes
        ValueType value_;                         // function value
        static DerivativeType zeroDerivative_;    // a zero object is required as return value
        //! This is function is called by constructors
        inline void initialise() {
#ifdef STATIC_REFERENCES
            for (size_t k=0; k<N_ARG; ++k) argRef_[k] = 0;
            for (size_t k=0; k<N_RES; ++k) resRef_[k] = 0;
#endif		
#ifdef DYNAMIC_REFERENCES
            arguments_ = 0;
            results_   = 0;
#endif		
        }

        //! Modifiers of the graph representation that add and remove edges
        //! method addResult is called only by method addArgument to ensure symetric references
        inline void addResult(const VARIABLE_TYPE &result, DerivativeType *p_derivative) {
#ifdef STATIC_REFERENCES
            //! first try static container
            for (size_t k=0; k<N_RES; ++k) {
                if (!resRef_[k]) {
                    resRef_[k] = &(const_cast<VARIABLE_TYPE &>(result));
                    resDrv_[k] = p_derivative;
                    return;
                }
            }
#endif
#ifdef DYNAMIC_REFERENCES
            //! if this does not work use dynamic container
            if (!results_) results_ = new ContainerType;
            //(const_cast<VARIABLE_TYPE &>(argument)).results_[this] = p_derivative;
            results_->insert(std::pair< VARIABLE_TYPE* ,DerivativeType* >(&(const_cast<VARIABLE_TYPE &>(result)),p_derivative));
            return;
#endif
            //! In case we don't have a dynamic container and we could not insert something is wrong
            MSG( "Variable Add Result Error: " << this << " " << &result )
                throw MinimADException("ERROR! Variable result can not be added.");	    			    
        } 

        //! method addArgumwnts is the interface for graph creation
        //! this is called by constructors, assignment operators, and elimitate() method
        inline void addArgument(const VARIABLE_TYPE &argument, const DerivativeType &derivative) {
#ifdef STATIC_REFERENCES
            //! find argument in static container and update derivative
            size_t idx = N_ARG;
            for (size_t k=0; k<N_ARG; ++k) {
                if (argRef_[k]==&argument) {
                    idx = k;
                    break;
                }
            }
            if (idx<N_ARG) {
                (*(argDrv_[idx])) += derivative;
                return;
            }
#endif
#ifdef DYNAMIC_REFERENCES
            //! also try the dynamic container
            if (arguments_) {
                typename ContainerType::iterator it = arguments_->find(&(const_cast<VARIABLE_TYPE &>(argument)));
                if (it!=arguments_->end()) {
                    (*it->second) += derivative;
                    return;
                }
            }
#endif
            //! now we have to add a new argument and its derivative to the graph
            DerivativeType *p_derivative = new DerivativeType(derivative);
#ifdef STATIC_REFERENCES    
            for (size_t k=0; k<N_ARG; ++k) {
                if (!argRef_[k]) {
                    argRef_[k] = &(const_cast<VARIABLE_TYPE &>(argument));
                    argDrv_[k] = p_derivative;
                    (const_cast<VARIABLE_TYPE &>(argument)).addResult(*this,p_derivative);
                    return;
                }
            }
#endif
#ifdef DYNAMIC_REFERENCES
            if (!arguments_) arguments_ = new ContainerType;
            //(*arguments_)[&(const_cast<VARIABLE_TYPE &>(argument))] = p_derivative;
            arguments_->insert( std::pair< VARIABLE_TYPE* ,DerivativeType* > (&(const_cast<VARIABLE_TYPE &>(argument)), p_derivative));
            (const_cast<VARIABLE_TYPE &>(argument)).addResult(*this,p_derivative);
            return;
#endif
            //! In case we don't have a dynamic container and we could not insert something is wrong
            MSG( "Variable Add Argument Error: " << this << " " << &argument )
                throw MinimADException("ERROR! Variable argument can not be added.");	    			    
        }
        //! method removeResults encapsulates the removal process of results called by elliminate()

        inline void removeResult(const VARIABLE_TYPE &result) {
#ifdef STATIC_REFERENCES
            for (size_t k=0; k<N_RES; ++k) {
                if (resRef_[k]==&result) {
                    resRef_[k] = 0;
                    return;
                }
            }
#endif
#ifdef DYNAMIC_REFERENCES
            if (results_) {
                typename ContainerType::iterator it = results_->find(&(const_cast<VARIABLE_TYPE &>(result)));
                if (it!=results_->end()){
                    results_->erase(it);
                    return;
                }
            }
#endif
            //! \warning if there is nothing to remove sonething is wrong with the function call
            MSG( "Variable Remove result Error: " << this << " " << &result )
                throw MinimADException("ERROR! Variable Result can not be removed.");	    		
        }

        //! method removeArguement encapsulates the removal process of arguments called by elliminate()
        inline void removeArgument(const VARIABLE_TYPE &argument) {
#ifdef STATIC_REFERENCES
            for (size_t k=0; k<N_ARG; ++k) {
                if (argRef_[k]==&argument) {
                    argRef_[k] = 0;
                    return;
                }
            }
#endif
#ifdef DYNAMIC_REFERENCES
            if (arguments_) {	    
                typename ContainerType::iterator it = arguments_->find(&(const_cast<VARIABLE_TYPE &>(argument)));
                if (it!=arguments_->end()) {
                    arguments_->erase(it);
                    return;
                }
            }
#endif	    
            //! \warning if there is nothing to remove sonething is wrong with the function call
            MSG( "Variable Remove argument Error: " << this << " " << &argument )
                throw MinimADException("ERROR! Variable Argument can not be removed.");	    		
        }

        inline void eliminate() {
#ifdef STATIC_REFERENCES
            //! first remove current node from static graph
            for (size_t k=0; k<N_RES; ++k) if (resRef_[k]) resRef_[k]->removeArgument( *this );
            for (size_t k=0; k<N_ARG; ++k) if (argRef_[k]) argRef_[k]->removeResult( *this );
            //! in a second step conect incoming and outgoing nodes applying the chain rule
            for (size_t i=0; i<N_RES; ++i) {
                if (resRef_[i]) {
                    for (size_t j=0; j<N_ARG; ++j) {
                        if (argRef_[j]) {
                            resRef_[i]->addArgument( *argRef_[j], (*resDrv_[i]) * (*argDrv_[j]) );
                            VARIABLE_TYPE::eliminations++;
                        }
                    }
                }
            }
#endif
#ifdef DYNAMIC_REFERENCES
            //! first remove current node from dynamic graph
            if (results_) for (typename ContainerType::iterator it=results_->begin(); it!=results_->end(); ++it) it->first->removeArgument( *this );
            if (arguments_) for (typename ContainerType::iterator it=arguments_->begin(); it!=arguments_->end(); ++it) it->first->removeResult( *this );
            //! in a second step connect incoming and outgoing nodes applying the chain rule
            if (arguments_ && results_) {
                for (typename ContainerType::iterator i=results_->begin(); i!=results_->end(); ++i) {
                    for (typename ContainerType::iterator j=arguments_->begin(); j!=arguments_->end(); ++j) {
                        i->first->addArgument( *(j->first), (*i->second) * (*j->second) ); 
                        VARIABLE_TYPE::eliminations++;
                    }
                }
            }
#endif
#ifdef STATIC_REFERENCES
#ifdef DYNAMIC_REFERENCES
            //! connect static results with dynamic arguments
            if (arguments_) {
                for (size_t i=0; i<N_RES; ++i) {
                    if (resRef_[i]) {
                        for (typename ContainerType::iterator j=arguments_->begin(); j!=arguments_->end(); ++j) {
                            resRef_[i]->addArgument( *(j->first), (*resDrv_[i]) * (*j->second) );
                            VARIABLE_TYPE::eliminations++;
                        }
                    }
                }
            }
            //! connect dynamic results with static arguments
            if (results_) {
                for (typename ContainerType::iterator i=results_->begin(); i!=results_->end(); ++i) {
                    for (size_t j=0; j<N_ARG; ++j) {
                        if (argRef_[j]) {
                            i->first->addArgument( *argRef_[j], (*i->second) * (*argDrv_[j]) ); 
                            VARIABLE_TYPE::eliminations++;			    
                        }
                    }
                }
            }
#endif // STATIC_REFERENCES
#endif // DYNAMIC_REFERENCES
            //! we have to clean up memory manually which is allocated in addArguments() and addResults()
            //! this is the place to modify the higher order elimination sequence
#ifdef STATIC_REFERENCES
            for (size_t k=0; k<N_RES; ++k) {
                if (resRef_[k]) {
                    delete resDrv_[k];
                    resRef_[k] = 0;
                }
            }
            for (size_t k=0; k<N_ARG; ++k) {
                if (argRef_[k]) {
                    delete argDrv_[k];
                    argRef_[k] = 0;
                }
            }
#endif
#ifdef DYNAMIC_REFERENCES
            if (results_) {
                for (typename ContainerType::iterator it=results_->begin(); it!=results_->end(); ++it) delete (it->second);
                results_->clear();
                delete results_;
                results_ = 0;
            }
            if (arguments_) {
                for (typename ContainerType::iterator it=arguments_->begin(); it!=arguments_->end(); ++it) delete (it->second);
                arguments_->clear();
                delete arguments_;
                arguments_ = 0;
            }
#endif
        } // eliminate()       

    public:
        static double eliminations;
        //! Constructors for Variable objects
        //! Internal unary constructor
        Variable( const ValueType&                            value,
                  const VARIABLE_TYPE&                        argument,
                  const DerivativeType&                       derivative )
            : value_(value) {
            MSG( "UC: " << this << " " << value << " " << &argument )
            initialise();
            this->addArgument( argument, derivative );		    
        }
        //! Internal binary constructor
        Variable( const ValueType&                            value,
                  const VARIABLE_TYPE&   			          argument1,
                  const VARIABLE_TYPE&                        argument2,
                  const DerivativeType&                       derivative1,
                  const DerivativeType&                       derivative2  )
            : value_(value) {
            MSG( "BC: " << this << " " << value << " " << &argument1 << " " << &argument2 )
            initialise();
            this->addArgument( argument1, derivative1 );		    
            this->addArgument( argument2, derivative2 );
        }	        	    
        //! Standard constructor
        Variable() : value_(0) {
            MSG( "SC: " << this )
            initialise();
        }
        //! ValueType constructor
        Variable( const ValueType &argument) 
            : value_(argument) {
            MSG( "VC: " << this << " " << argument )
            initialise();
        }
        //! Copy constructor
        Variable( const VARIABLE_TYPE &argument ) 
            : value_(argument.value_) {
            MSG( "CC: " << this << " " << argument.value_ << " " << &argument )
            initialise();
            this->addArgument( argument, (DerivativeType)1);
        }            
        //! Destructor
        ~Variable() { 
            MSG( "DR: " << this )
            eliminate();
        }	
        //! Inspectors
        //! Returning a reference to the value component
        inline const ValueType& value() const { return value_; }
        inline const ValueType& val() const { return value_; }	
        //! Returning a reference to the derivative object of this w.r.t. y
        inline const DerivativeType& operator % (const VARIABLE_TYPE &y ) const {
#ifdef STATIC_REFERENCES
            //! search for edge in static container
            for (size_t k=0; k<N_ARG; ++k) {
                if (argRef_[k]==&y) {
                    return *argDrv_[k];
                }
            }
#endif            
#ifdef DYNAMIC_REFERENCES
            if (arguments_) {
                typename ContainerType::iterator arg_it = (const_cast< VARIABLE_TYPE* >(this))->arguments_->find(&(const_cast< VARIABLE_TYPE& >(y)));
                if (arg_it!=this->arguments_->end()) return (*arg_it->second);
            }
#endif
            /*
            typename ContainerType::iterator arg_it = (const_cast< VARIABLE_TYPE* >(this))->arguments_->find(&(const_cast< VARIABLE_TYPE& >(y)));
            if (arg_it!=this->arguments_->end()) {
            MSG("arg: " << this << " derivative: " << (*arg_it->second) )
            }
            typename ContainerType::iterator res_it = (const_cast< VARIABLE_TYPE& >(y)).results_->find((const_cast< VARIABLE_TYPE* >(this))); 
            if (res_it!=y.results_->end()) {
            MSG("res: " << &y << " derivative: " << (*res_it->second) )		
            }
            if (arg_it!=this->arguments_->end()) return (*arg_it->second);
            if (res_it!=y.results_->end()) return (*res_it->second);
            */
            return zeroDerivative_;
        }	
        //! Assignment operators
        //! Variable = Variable assignment
        inline VARIABLE_TYPE& operator = (const VARIABLE_TYPE& argument) {
            MSG( "Variable Active Assignment: " << this << " " << &argument )
            if (this != &argument) {
                this->eliminate();
                this->value_ = argument.value_;
                this->addArgument( argument, (DerivativeType)1);
            }
            // otherwise do nothing
            return *this;
        }
        //! Variable = ValueType assignment
        inline VARIABLE_TYPE& operator = (const ValueType& argument) {
            MSG( "Variable ValueType Assignment: " << this << " " << argument )
            this->eliminate();
            this->value_ = argument;
            return *this;
        }
        //! Arithmetic-Assignment operators
        //! Variable (OP)= Variable
        inline VARIABLE_TYPE& operator += (const VARIABLE_TYPE &x) { return *this = *this + x; }
        inline VARIABLE_TYPE& operator -= (const VARIABLE_TYPE &x) { return *this = *this - x; }
        inline VARIABLE_TYPE& operator *= (const VARIABLE_TYPE &x) { return *this = *this * x; }
        inline VARIABLE_TYPE& operator /= (const VARIABLE_TYPE &x) { return *this = *this / x; }
        //! Variable (OP)= ValueType
        inline VARIABLE_TYPE& operator += (const ValueType &x) { return *this = *this + x; }
        inline VARIABLE_TYPE& operator -= (const ValueType &x) { return *this = *this - x; }
        inline VARIABLE_TYPE& operator *= (const ValueType &x) { return *this = *this * x; }
        inline VARIABLE_TYPE& operator /= (const ValueType &x) { return *this = *this / x; }        
        //! Relational operators
        //! Variable x Variable
        inline bool operator >  (const VARIABLE_TYPE &x) { return this->value_ >  x.value_; }
        inline bool operator >= (const VARIABLE_TYPE &x) { return this->value_ >= x.value_; }
        inline bool operator <  (const VARIABLE_TYPE &x) { return this->value_ <  x.value_; }
        inline bool operator <= (const VARIABLE_TYPE &x) { return this->value_ <= x.value_; }
        inline bool operator != (const VARIABLE_TYPE &x) { return this->value_ != x.value_; }
        inline bool operator == (const VARIABLE_TYPE &x) { return this->value_ == x.value_; }
        //! Variable x Variable (for left const operands)
        inline friend bool operator >  (const VARIABLE_TYPE &x, const VARIABLE_TYPE &y) { return x.value_ >  y.value_; }
        inline friend bool operator >= (const VARIABLE_TYPE &x, const VARIABLE_TYPE &y) { return x.value_ >= y.value_; }
        inline friend bool operator <  (const VARIABLE_TYPE &x, const VARIABLE_TYPE &y) { return x.value_ <  y.value_; }
        inline friend bool operator <= (const VARIABLE_TYPE &x, const VARIABLE_TYPE &y) { return x.value_ <= y.value_; }
        inline friend bool operator != (const VARIABLE_TYPE &x, const VARIABLE_TYPE &y) { return x.value_ != y.value_; }
        inline friend bool operator == (const VARIABLE_TYPE &x, const VARIABLE_TYPE &y) { return x.value_ == y.value_; }
        //! Variable x ValueType
        inline bool operator >  (const ValueType &x) { return this->value_ >  x; }
        inline bool operator >= (const ValueType &x) { return this->value_ >= x; }
        inline bool operator <  (const ValueType &x) { return this->value_ <  x; }
        inline bool operator <= (const ValueType &x) { return this->value_ <= x; }
        inline bool operator != (const ValueType &x) { return this->value_ != x; }
        inline bool operator == (const ValueType &x) { return this->value_ == x; }
        //! ValueType x Variable
        inline friend bool operator >  (const ValueType &x, const VARIABLE_TYPE &y) { return x >  y.value_; }
        inline friend bool operator >= (const ValueType &x, const VARIABLE_TYPE &y) { return x >= y.value_; }
        inline friend bool operator <  (const ValueType &x, const VARIABLE_TYPE &y) { return x <  y.value_; }
        inline friend bool operator <= (const ValueType &x, const VARIABLE_TYPE &y) { return x <= y.value_; }
        inline friend bool operator != (const ValueType &x, const VARIABLE_TYPE &y) { return x != y.value_; }
        inline friend bool operator == (const ValueType &x, const VARIABLE_TYPE &y) { return x == y.value_; }
        //! IO steaming
        //! Apply IO only to the value component
        inline friend std::ostream& operator << (std::ostream &output, const VARIABLE_TYPE &x) {
            return output << x.value_;
        }
        //! Binary arithmetic operators
        //! Variable x Variable
        inline VARIABLE_TYPE operator + (const VARIABLE_TYPE &y ) {
            return VARIABLE_TYPE( this->value_ + y.value_, *this, y, (DerivativeType)1, (DerivativeType)1 );
        }
        inline VARIABLE_TYPE operator - (const VARIABLE_TYPE &y ) {
            return VARIABLE_TYPE( this->value_ - y.value_, *this, y, (DerivativeType)1, (DerivativeType)(-1) );
        }
        inline VARIABLE_TYPE operator * (const VARIABLE_TYPE &y ) {
            return VARIABLE_TYPE( this->value_ * y.value_, *this, y, y.value_, this->value_ );
        }
        inline VARIABLE_TYPE operator / (const VARIABLE_TYPE &y ) {
            ValueType tmp = this->value_ / y.value_;
            return VARIABLE_TYPE( tmp, *this, y, 1/y.value_, -tmp/y.value_ );
        }
        //! Variable x Variable (for left const operands)
        inline friend VARIABLE_TYPE operator + (const VARIABLE_TYPE &x, const VARIABLE_TYPE &y ) {
            return VARIABLE_TYPE( x.value_ + y.value_, x, y, (DerivativeType)1, (DerivativeType)1 );
        }
        inline friend VARIABLE_TYPE operator - (const VARIABLE_TYPE &x, const VARIABLE_TYPE &y ) {
            return VARIABLE_TYPE( x.value_ - y.value_, x, y, (DerivativeType)1, (DerivativeType)(-1) );
        }
        inline friend VARIABLE_TYPE operator * (const VARIABLE_TYPE &x, const VARIABLE_TYPE &y ) {
            return VARIABLE_TYPE( x.value_ * y.value_, x, y, y.value_, x.value_ );
        }
        inline friend VARIABLE_TYPE operator / (const VARIABLE_TYPE &x, const VARIABLE_TYPE &y ) {
            ValueType tmp = x.value_ / y.value_;
            return VARIABLE_TYPE( tmp, x, y, 1/y.value_, -tmp/y.value_ );
        }

        //! Variable x ValueType
        inline VARIABLE_TYPE operator + (const ValueType &y ) {
            return VARIABLE_TYPE( this->value_ + y, *this, (DerivativeType)1 );
        }
        inline VARIABLE_TYPE operator - (const ValueType &y ) {
            return VARIABLE_TYPE( this->value_ - y, *this, (DerivativeType)1 );
        }
        inline VARIABLE_TYPE operator * (const ValueType &y ) {
            return VARIABLE_TYPE( this->value_ * y, *this, y );
        }
        inline VARIABLE_TYPE operator / (const ValueType &y ) {
            return VARIABLE_TYPE( this->value_ / y, *this, 1/y );
        }
        //! ValueType x Variable
        inline friend VARIABLE_TYPE operator + (const ValueType &x, const VARIABLE_TYPE &y ) {
            return VARIABLE_TYPE( x + y.value_, y, (DerivativeType)1 );
        }
        inline friend VARIABLE_TYPE operator - (const ValueType &x, const VARIABLE_TYPE &y ) {
            return VARIABLE_TYPE( x - y.value_, y, (DerivativeType)(-1) );
        }
        inline friend VARIABLE_TYPE operator * (const ValueType &x, const VARIABLE_TYPE &y ) {
            return VARIABLE_TYPE( x * y.value_, y, x );
        }
        inline friend VARIABLE_TYPE operator / (const ValueType &x, const VARIABLE_TYPE &y ) {
            ValueType tmp = x / y.value_;
            return VARIABLE_TYPE( tmp, y, -tmp/y.value_ );
        }
        //! Unary arithmetic operators
        inline friend VARIABLE_TYPE operator + (const VARIABLE_TYPE &x) {
            return VARIABLE_TYPE( x.value_, x, (DerivativeType)1 );
        }
        inline friend VARIABLE_TYPE operator - (const VARIABLE_TYPE &x) {
            return VARIABLE_TYPE( -x.value_, x, (DerivativeType)(-1) );
        }
        //! Intrinsic functions
        inline friend VARIABLE_TYPE fabs(const VARIABLE_TYPE &x) {
            return (x.value_>=0.0) ? x : (-x);
        }
        inline friend VARIABLE_TYPE exp(const VARIABLE_TYPE &x) {
            ValueType tmp = exp(x.value_);
            return VARIABLE_TYPE( tmp, x, tmp );
        }
        inline friend VARIABLE_TYPE log(const VARIABLE_TYPE &x) {
            return VARIABLE_TYPE( log(x.value_), x, 1/x.value_ );
        }
        inline friend VARIABLE_TYPE sqrt(const VARIABLE_TYPE &x) {
            ValueType tmp = sqrt(x.value_);
            return VARIABLE_TYPE( tmp, x, 0.5/tmp );
        }
        inline friend VARIABLE_TYPE sin(const VARIABLE_TYPE &x) {
            return VARIABLE_TYPE( sin(x.value_), x, cos(x.value_) );
        }
        inline friend VARIABLE_TYPE cos(const VARIABLE_TYPE &x) {
            return VARIABLE_TYPE( cos(x.value_), x, -sin(x.value_) );
        }
        inline friend VARIABLE_TYPE tan(const VARIABLE_TYPE &x) {
            ValueType tmp = tan(x.value_);
            return VARIABLE_TYPE( tmp, x, 1+tmp*tmp );
        }
        inline friend VARIABLE_TYPE erf(const VARIABLE_TYPE &x) {
            return VARIABLE_TYPE( erf(x.value_), x, M_2_SQRTPI * exp(-x.value_*x.value_) );
        }
        inline friend VARIABLE_TYPE erfc(const VARIABLE_TYPE &x) {
            return VARIABLE_TYPE( erfc(x.value_), x, (-M_2_SQRTPI) * exp(-x.value_*x.value_) );
        }
        inline friend VARIABLE_TYPE atan2(const VARIABLE_TYPE &y, const VARIABLE_TYPE &x) {
            ValueType den = x.value_ * x.value_ + y.value_ * y.value_;
            return VARIABLE_TYPE( atan2(y.value_,x.value_), y, x, x.value_/den, -y.value_/den );
        }
    };   // class Variable 

    //! Definition of static class attributes
    TEMPLATE_FUNCT_VARIABLE
    DerivativeType VARIABLE_TYPE::zeroDerivative_(0);
    TEMPLATE_FUNCT_VARIABLE
    double VARIABLE_TYPE::eliminations = 0;

    //! Double cast
    template <class Type> inline double value        ( Type   x ) { return value(x.value()); }
    template <>           inline double value<double>( double x ) { return x;                }


    //! Additional double x Variable operators; required for more complex Variable instantiations
    //! Variable x double
    TEMPLATE_FUNCT_VARIABLE
    inline VARIABLE_TYPE operator + (const double &x, const VARIABLE_TYPE &y ) {
        return VARIABLE_TYPE( x + y.value(), y, (DerivativeType)1 );
    }
    TEMPLATE_FUNCT_VARIABLE
    inline VARIABLE_TYPE operator - (const double &x, const VARIABLE_TYPE &y ) {
        return VARIABLE_TYPE( x - y.value(), y, (DerivativeType)(-1) );
    }
    TEMPLATE_FUNCT_VARIABLE
    inline VARIABLE_TYPE operator * (const double &x, const VARIABLE_TYPE &y ) {
        return VARIABLE_TYPE( x * y.value(), y, (DerivativeType)(x) );
    }
    TEMPLATE_FUNCT_VARIABLE
    inline VARIABLE_TYPE operator / (const double &x, const VARIABLE_TYPE &y ) {
        ValueType tmp = x / y.value();
        return VARIABLE_TYPE( tmp, y, -tmp/y.value() );
    }    
    //! Double x Variable
    TEMPLATE_FUNCT_VARIABLE
    inline VARIABLE_TYPE operator + (const VARIABLE_TYPE &x, const double &y ) {
        return VARIABLE_TYPE( x.value() + y, x, (DerivativeType)1 );
    }
    TEMPLATE_FUNCT_VARIABLE
    inline VARIABLE_TYPE operator - (const VARIABLE_TYPE &x, const double &y ) {
        return VARIABLE_TYPE( x.value() - y, x, (DerivativeType)1 );
    } 
    TEMPLATE_FUNCT_VARIABLE
    inline VARIABLE_TYPE operator * (const VARIABLE_TYPE &x, const double &y ) {
        return VARIABLE_TYPE( x.value() * y, x, (DerivativeType)y );
    } 
    TEMPLATE_FUNCT_VARIABLE
    inline VARIABLE_TYPE operator / (const VARIABLE_TYPE &x, const double &y ) {
        return VARIABLE_TYPE( x.value() / y, x, (DerivativeType)1/y );
    }

#undef TEMPLATE_CLASS_VARIABLE
#undef TEMPLATE_FUNCT_VARIABLE
#undef VARIABLE_TYPE

}   // namespace MinimAD

#endif   // MinimAD_MinimADVariable_hpp

#ifdef EXECUTABLE

// Dummy userdefined double type
template <class Type>
class MyDouble {
public:
    Type value_;
    //! Standard constructor
    MyDouble() : value_(0) { }
    //! ValueType constructor
    MyDouble( const Type &argument) : value_(argument) { }
    //! Copy constructor
    MyDouble( const MyDouble<Type> &argument ) : value_(argument.value_) { }            
};


// 2D Rosenbrock function
// f(x,y)    = 100(x^2 - y)^2 + (x - 1)^2
// df/dx     = 400x(x^2 - y) + 2(x - 1)
// df/dy     = -200(x^2 - y)
// d^2f/dx^2 = 400(3x^2 - y) + 2
// d^2f/dxdy = -400x = d^2f/dydxs
// d^2f/dy^2 = 200
template <class Type>
Type rosenbrockFunction(Type& x, Type& y) {    
    return (x*x - y)*(x*x - y)*100 + (x - 1.0)*(x - 1.0);
}

#include <time.h>

template <class Type>
void poissonTest(std::vector<Type>& x, std::vector<Type>& y, size_t dim) {
    if (dim<2) return;
    Type a_(-1.0), b_(2.0), c_(-1.0);
    std::vector<Type> z(dim), a(dim), b(dim), c(dim);
    x.resize(dim);
    y.resize(dim);
    // initialisation
    for (size_t i=0; i<dim; ++i) {
        x[i] = 1.0;
        a[i] = a_;
        b[i] = b_;
        c[i] = c_;
    }
    // in place LU decomposition; no error handling if LU decomposition does not exist
    for (size_t i=1; i<dim; ++i) {
        a[i] /= b[i-1];
        b[i] -= c[i-1]*a[i];
    }
    // forward substitution
    z[0] = x[0];
    for (size_t i=1; i<dim; ++i) z[i] = x[i] - a[i]*z[i-1];
    // backward substitution
    y[dim-1] = z[dim-1]/b[dim-1];
    for (int i=dim-2; i>=0; --i) y[i] = (z[i] - c[i]*y[i+1])/b[i];
    // elimination
    //for (size_t i=0; i<dim; ++i) c[i] = MinimAD::value(c[i]);
    //for (size_t i=0; i<dim; ++i) b[i] = 0.0;
    //for (size_t i=0; i<dim; ++i) a[i] = 0.0;
    //for (size_t i=0; i<dim; ++i) x[i] = 0.0;
    //for (size_t i=dim; i>0; --i) z[i-1] = 0.0;
    // multiplication
    /*
    z[0] = 2.0*y[0] - y[1];
    for (size_t i=1; i<dim-1; ++i)  z[i] = -y[i-1] + 2.0*y[i] - y[i+1];
    z[dim-1] = -y[dim-2] + 2.0*y[dim-1];
    // copy result
    for (size_t i=0; i<dim; ++i) y[i] = z[i];
    */
    return;
}


void run_poissonTest(size_t base, size_t count) {
    std::vector< MinimAD::Variable<double> > x, y;
    clock_t start, end;
    size_t dim=base;
    double res;
    std::cout << "Poisson Test" << std::endl;
    for (size_t k=0; k<count; ++k) {
        MinimAD::Variable<double>::eliminations = 0.0;
        res = 0.0;
        start = clock();
        poissonTest(x,y,dim);
        for (size_t i=0; i<dim; ++i) {
            for (size_t j=0; j<dim; ++j) {
                res += y[i] % x[j];
            }
        }
        end = clock();
        std::cout << "  dim = "      << dim 
            << "  res = "      << res
            << "  elims = "    << MinimAD::Variable<double>::eliminations
            << "  CPU-sec. = " << 1.0*(end - start)/CLOCKS_PER_SEC
            << std::endl;
        dim *= base;
    }
}


template <typename Type> inline
    Type Phi( Type x) { return 0.5*(erf(x*M_SQRT1_2)+1.0); }

template <typename Type> inline
    Type PhiInv(const Type x) { return M_SQRT2*erf_inv(2.0*x-1.0); }


template <typename Type> inline
    Type normalExpectation( std::vector<Type>&  v,  //  payoff
    size_t              dim
    ) {
        std::vector<double> x(dim);  //  grid points of payoff
        double x0=0, s=3.0;
        for (size_t i=0; i<dim; ++i) x[i] = (2.0*i/(dim-1)-1.0)*s + x0;	
        v.resize(dim);
        for (size_t i=0; i<dim; ++i) v[i] = 1.0;
        Type  mu(0.0), var(1.0), res;
        std::vector<Type> sums(x.size());
        sums[0] = (Type)0.0;
        Type Q1, Q2 = Phi((x[0]-mu)/sqrt(var));
        for (size_t i=0; i<x.size()-1; ++i) {
            Q1 = Q2;
            Q2 = Phi((x[i+1]-mu)/sqrt(var));
            sums[i+1] = sums[i] + v[i+1]*Q2 - v[i]*Q1 - 0.5*(Q1 + Q2)*(v[i+1] - v[i]) ;
            //std::cout << " i="    << i << "  " << "elims=" << MinimAD::Variable<double>::eliminations << std::endl;
        }
        res = sums[dim-1];
        Q2 = 0;
        Q1 = 0;
        for (size_t i=dim; i>0; --i) sums[i-1] = 0.0;
        //std::cout << " i="    << i << "  " << "elims=" << MinimAD::Variable<double>::eliminations << std::endl;
        return res;
}

void run_NormalExpectation(size_t base, size_t count) {
    std::vector< MinimAD::Variable<double> > v;
    MinimAD::Variable<double> y;
    std::vector<double> a;
    double b;
    clock_t start, end;
    size_t dim=base;
    double res;
    std::cout << "Normal Expectation Test double" << std::endl;
    for (size_t k=0; k<count; ++k) {
        start = clock();
        b = normalExpectation(a,dim);
        end = clock();
        std::cout << "  dim = "      << dim 
            << "  CPU-sec. = " << 1.0*(end - start)/CLOCKS_PER_SEC
            << std::endl;
        dim *= base;
    }
    dim=base;
    std::cout << "Normal Expectation Test Variable<double>" << std::endl;
    for (size_t k=0; k<count; ++k) {
        res = 0.0;
        start = clock();
        y = normalExpectation(v,dim);
        for (size_t i=0; i<dim; ++i) res += y % v[i];
        end = clock();
        std::cout << "  dim = "      << dim 
            << "  res = "      << res
            << "  elims = "    << MinimAD::Variable<double>::eliminations
            << "  CPU-sec. = " << 1.0*(end - start)/CLOCKS_PER_SEC
            << std::endl;
        dim *= base;
    }
}

template <typename Type> inline
    Type fibonacci(Type x1, Type x2, size_t n) {
        switch(n) {
        case 0 : return 0;
        case 1 : return x1;
        case 2 : return x2;
        }
        for (size_t k=3; k<=n; ++k) {
            if (k%2) x1 += x2/x1;
            else     x2 += x1/x2;
        }
        if (n%2) return x1;
        else     return x2;
}

void run_Fibonacci(size_t base, size_t count) {
    clock_t start, end;
    size_t n;
    std::cout << "Fibonacci Test double" << std::endl;
    n=base;
    for (size_t k=0; k<count; ++k) {
        double x1=1.0, x2=2.0, z;
        start = clock();
        z = fibonacci(x1,x2,n);
        end = clock();
        std::cout << "  n = "      << n
            << "  z = "      << z	    
            << "  CPU-sec. = " << 1.0*(end - start)/CLOCKS_PER_SEC
            << std::endl;
        n *= base;
    }
    std::cout << "Fibonacci Test  Variable<double>" << std::endl;
    n=base;
    for (size_t k=0; k<count; ++k) {
        MinimAD::Variable<double>::eliminations = 0;
        MinimAD::Variable<double> x1=1.0, x2=2.0, z;
        start = clock();
        z = fibonacci(x1,x2,n);
        end = clock();
        std::cout << "  n = "      << n 
            << "  z = "      << z
            << "  elims = "    << MinimAD::Variable<double>::eliminations
            << "  CPU-sec. = " << 1.0*(end - start)/CLOCKS_PER_SEC
            << std::endl;
        n *= base;
    }

}

int main() {
    std::cout << "MinimAD (C) Sebastian Schlenktich (2011)"<< std::endl;
    std::cout << std::endl;

    // Test case 1, Rosenbrock function
    typedef MinimAD::Variable<double>  Active1;
    typedef MinimAD::Variable<Active1> Active2;
    typedef MinimAD::Variable<Active2> Active3;
    Active2 t;
    t = 1;
    Active3 x, y, z;
    std::cout << "+++ Text case 1 +++" <<std::endl;

    std::cout << "Rosenbrock function evaluation z = f(x,y) = 100(x^2 - y)^2 + (x - 1)^2 with x = "
        << (x=t) << " and y = " << (y=t) << " (global minimum)" << std::endl;
    z = rosenbrockFunction(x, y);
    std::cout << std::endl;

    std::cout << "objective function z = " << z << std::endl;
    std::cout << std::endl;

    std::cout << "gradient [ df/dx, df/dy ] = [ (z % x)  (z % y) ] = [ "
        << (z % x) << "  " << (z % y) << " ]" << std::endl;
    std::cout << std::endl;

    std::cout << "Hessian:" << std::endl;
    std::cout << "  [ d^2f/dx^2  d^2f/dydx ] = [ ((z % x) % x.value())  ((z % y) % x.value()) ] "
        << "= [ " << std::setw(4) << ((z % x) % x.value()) << "  " << std::setw(4) << ((z % y) % x.value()) << " ]" << std::endl;
    std::cout << "  [ d^2f/dxdy  d^2f/dy^2 ] = [ ((z % x) % y.value())  ((z % y) % y.value()) ] "
        << "= [ " << std::setw(4) << ((z % x) % y.value()) << "  " << std::setw(4) << ((z % y) % y.value()) << " ]" << std::endl;
    std::cout << std::endl;

    std::cout << "Third order derivative: d^3f/dx^3 =  (((z % x) % x.value()) % x.value().value())  " << (((z % x) % x.value()) % x.value().value()) << std::endl;
    std::cout << std::endl;

    std::cout << "+++ Text case 2 +++" <<std::endl;
    t = M_PI/6;    

    std::cout << "Trigonometric function evaluation y = f(x) = sin(x)/cos(x) with x = " << (x=t) << " = Pi/6" << std::endl;
    y = sin(x)/cos(x);
    std::cout << "y       = " << y                                           << std::endl;
    std::cout << "f'(x)   = " << (y % x)                                     << std::endl;
    std::cout << "f''(x)  = " << ((y % x) % x.value())                       << std::endl;
    std::cout << "f'''(x) = " << (((y % x) % x.value()) % x.value().value()) << std::endl;
    std::cout << "Check for consistency y = f(x) = tan(x) with x = " << (x=t) << " = Pi/6" << std::endl;
    y = tan(x);
    std::cout << "y       = " << y << std::endl;
    std::cout << "f'(x)   = " << (y % x)                                     << std::endl;
    std::cout << "f''(x)  = " << ((y % x) % x.value())                       << std::endl;
    std::cout << "f'''(x) = " << (((y % x) % x.value()) % x.value().value()) << std::endl;
    std::cout << std::endl;    

    std::cout << "+++ Text case 3 +++" <<std::endl;
    t = 1;    
    std::cout << "Exponential function evaluation y = f(x) = exp(x) with x = " << (x=t) << std::endl;
    y = exp(x);
    std::cout << "y       = " << y                                           << std::endl;
    std::cout << "f'(x)   = " << (y % x)                                     << std::endl;
    std::cout << "f''(x)  = " << ((y % x) % x.value())                       << std::endl;
    std::cout << "f'''(x) = " << (((y % x) % x.value()) % x.value().value()) << std::endl;
    std::cout << std::endl;    
    {
        std::cout << "+++ Text case 4 +++" <<std::endl;
        MinimAD::Variable<double> x, y, z;
        double h;
        std::cout << "g(y,x) = atan2(y,x) function with y = " << (y=1) << ", x = " << (x=1)
            << ", and h = " << (h=1.0e-8) << std::endl;
        std::cout << "double   [ g, g_y, g_x ] = [ "
            << (atan2(y.value(),x.value())) << ", "
            << ((atan2(y.value()+h,x.value())-atan2(y.value()-h,x.value()))/h/2) << ", "
            << ((atan2(y.value(),x.value()+h)-atan2(y.value(),x.value()-h))/h/2) << "] "
            << std::endl;
        z=atan2(y,x);
        std::cout << "Variable [ g, g_y, g_x ] = [ "
            << (z) << ", " << (z%y) << ", " << (z%x) << "] " << std::endl;
    }
    /*
    Active2 t1, t2, t3;
    t1 = 3; t2 = 2;
    t3 = t1 * t2;
    std::cout << "(t3 % t1) = " << (t3 % t1) << std::endl;
    std::cout << "(t3 % t2) = " << (t3 % t2) << std::endl;
    std::cout << "((t3 % t1) % t1.value()) = " << ((t3 % t1) % t1.value()) << std::endl;
    std::cout << "((t3 % t2) % t2.value()) = " << ((t3 % t2) % t2.value()) << std::endl;
    std::cout << "((t3 % t1) % t2.value()) = " << ((t3 % t1) % t2.value()) << std::endl;

    for (size_t k=0; k<10; ++k) 
    std::cout << "fibonacci(1,2,"<< k <<") = " << fibonacci(1,2,k) << std::endl;

    //run_NormalExpectation(2, 16);
    run_Fibonacci(2, 20);
    */
    return 0;
}

#endif   // EXECUTABLE

/*  Example output

d90329@t4nb117 ~/my/PROJECTS/MinimAD
$ g++ -o MinimAD MinimADVariable.cpp 

d90329@t4nb117 ~/my/PROJECTS/MinimAD
$ ./MinimAD
MinimAD (C) Sebastian Schlenktich (2011)

+++ Text case 1 +++
Rosenbrock function evaluation z = f(x,y) = 100(x^2 - y)^2 + (x - 1)^2 with x = 1 and y = 1 (global minimum)

objective function z = 0

gradient [ df/dx, df/dy ] = [ (z % x)  (z % y) ] = [ 0  -0 ]

Hessian:
[ d^2f/dx^2  d^2f/dydx ] = [ ((z % x) % x.value())  ((z % y) % x.value()) ] = [  802  -400 ]
[ d^2f/dxdy  d^2f/dy^2 ] = [ ((z % x) % y.value())  ((z % y) % y.value()) ] = [ -400   200 ]

Third order derivative: d^3f/dx^3 =  (((z % x) % x.value()) % x.value().value())  2400

+++ Text case 2 +++
Trigonometric function evaluation y = f(x) = sin(x)/cos(x) with x = 0.523599 = Pi/6
y       = 0.57735
f'(x)   = 1.33333
f''(x)  = 1.5396
f'''(x) = 5.33333
Check for consistency y = f(x) = tan(x) with x = 0.523599 = Pi/6
y       = 0.57735
f'(x)   = 1.33333
f''(x)  = 1.5396
f'''(x) = 5.33333

+++ Text case 3 +++
Exponential function evaluation y = f(x) = exp(x) with x = 1
y       = 2.71828
f'(x)   = 2.71828
f''(x)  = 2.71828
f'''(x) = 2.71828

*/
