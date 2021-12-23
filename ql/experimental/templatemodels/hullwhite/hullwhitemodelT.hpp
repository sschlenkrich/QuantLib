/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2010 Sebastian Schlenkrich

*/

/*! \file hullwhitemodelT.hpp
    \brief Analytic pricing formulas in the Hull White model
               dr(t) = [a(t) - b r(t)]dt + sigma(t) dW(t)
           The model assumes
             - piecewise left-constant volatility sigma given as a std:vector of
               (active) reals
             - continuous forward rate interpololation in yield curve
             - time given as day count independent year fractions
           The interaction with instruments is realised by BondOptionEngine pricing engine
           All methods are template based to allow incorporation of Automatic Differentiation
           tools
*/


#ifndef quantlib_templatehullwhitemodel_hpp
#define quantlib_templatehullwhitemodel_hpp

#include <ql/termstructures/yieldtermstructure.hpp>
#include <ql/option.hpp>
#include <ql/experimental/templatemodels/auxilliaries/auxilliariesT.hpp>
#include <ql/experimental/templatemodels/auxilliaries/integratorsT.hpp>


namespace QuantLib {

    class TemplateModel : public virtual Observable { };

    // Declaration of the Hull White model class
    template <class DateType, class PassiveType, class ActiveType>
    class HullWhiteModelT : public TemplateModel {
    protected:
        // attributes defining the model
        Handle<YieldTermStructure> termStructure_;  // the yield curve is assumed to be passive
        PassiveType                mean_;           // constant mean reversion parameter
        std::vector<DateType>      volaTimes_;      // time-grid of left-constant short rate volatility values
        std::vector<ActiveType>    volaValues_;     // volatility values
        size_t                     Nvol_;           // number of volatility values
        // calibration reference prices
        bool                       calibrateMean_;      // calibrate mean reversion parameter
        std::vector<PassiveType>   bermudanB76Prices_;  // reference Black'76 prices
        std::vector<PassiveType>   calibIters_;         // iterations required for convergence
        // attributes of bermudan option
        bool                       estimateAccuracy_;
        std::vector<PassiveType>   shortRateGrid_;
        std::vector<ActiveType>    europeansAnalytical_;
        std::vector<ActiveType>    europeansNumerical_;  // only evaluated if estimateAccuracy_ == true
        // auxilliary functions
        // set up the short rate grid
        void evaluateShortRateGrid( PassiveType r0,      // center of short rate grid
                                    PassiveType s,       // distance to boundaries
                                    size_t dim );        // number of grid points
        // helper class
        class func_F {
            PassiveType mean_;
            DateType T_;
        public:
            func_F( PassiveType mean, DateType T) : mean_(mean), T_(T) {}
            PassiveType operator() ( PassiveType u ) { return exp(-2.0*mean_*(T_ - u)) / 2.0 / mean_; }
        };
    public:
        // Construct a Hull White model by passing the attributes
        HullWhiteModelT( const Handle<YieldTermStructure>& termStructure,
                         const PassiveType                mean,
                         const std::vector<DateType>&     volaTimes,
                         const std::vector<ActiveType>&   volaValues );

        // Inspectors
        const Handle<YieldTermStructure>& termStructure() const { return termStructure_; }
        bool                              calibrateMean() const { return calibrateMean_; }
        PassiveType                       mean() const { return mean_; }
        const std::vector<DateType>&      volaTimes() const { return volaTimes_; }
        const std::vector<ActiveType>&    volaValues() const { return volaValues_; }
        const std::vector<PassiveType>&   bermudanB76Prices() const { return bermudanB76Prices_; }
        const std::vector<PassiveType>&   calibIters() const { return calibIters_; }
        bool                              estimateAccuracy() const { return estimateAccuracy_; }
        const std::vector<ActiveType>&    europeansAnalytical() const { return europeansAnalytical_; }
        const std::vector<ActiveType>&    europeansNumerical() const { return europeansNumerical_; }
        // Change settings
        void                              setCalibrateMean(bool calibrateMean) { calibrateMean_ = calibrateMean; }
        void                              setEstimateAccuracy(bool estimateAccuracy) { estimateAccuracy_ = estimateAccuracy; }

        virtual
        ActiveType ZeroBond( const DateType    settlement,
                             const DateType    maturity,
                             const ActiveType  shortRate,
                             ActiveType        *dBond_dr = 0 );
     
        virtual
        ActiveType CouponBond( const DateType                  settlement,   // bond settlement date in year fractions
                               const std::vector<DateType>&    payTimes,     // pay dates in year fraction of coupon period
                               const std::vector<PassiveType>& cashFlows,    // fixed coupon payments (absolut value)
                               const ActiveType                shortRate,
                               ActiveType                      *dBond_dr = 0,
                               const size_t                    idx_start = 0); // consider only coupons idx_start,...,*.size()

        virtual
        ActiveType ZeroBondOption( const DateType     excercise,     // option's excercise date (equal settlment)
                                   const ActiveType   strike,        // strike payed at excercise date
                                   const DateType     maturity,      // payment date of notional 1
                                   const Option::Type cop);          // call (1) or put (-1) option

        virtual
        ActiveType CouponBondOption( const DateType                  excercise,    // option's excercise date (equal settlment)
                                     const PassiveType               strike,       // strike payed at excercise date
                                     /* option's underlying */
                                     const std::vector<DateType>&     startTimes,   // start dates of coupon period
                                     const std::vector<DateType>&     payTimes,     // pay dates of coupon perid
                                     const std::vector<PassiveType>&  cashFlows,    // fixed coupon payments (absolut value)
                                     const Option::Type               cop);         // call (1) or put (-1) option

        virtual
        ActiveType BermudanBondOption( const std::vector<DateType>&     exercTimes,    // option's exercise dates (equal settlment)
                                       const std::vector<PassiveType>&  strikeValues,  // strike payed at excercise dates
                                       // option's underlying
                                       const std::vector<DateType>&     startTimes,   // start dates of coupon period
                                       const std::vector<DateType>&     payTimes,     // pay dates of coupon perid
                                       const std::vector<PassiveType>&  cashFlows,    // fixed coupon payments (absolut value)
                                       const Option::Type               cop,          // call (1) or put (-1) option
                                       // discretisation properties
                                       const size_t                     dim,          // number of short rate grid points
                                       const PassiveType                gridRadius,   // radius s of short rate grid [r0-s, r0+s]
                                       const PassiveType                tol );        // absolute tolerance for numerical integration

        virtual
        const std::vector<ActiveType>& CalibrateVolatility (
                                       const std::vector<DateType>&     exercTimes,   // option's exercise dates (equal settlment)
                                       const std::vector<PassiveType>&  strikeValues, // strike payed at excercise dates
                                       const std::vector<PassiveType>&  b76Prices,    // reference European prices
                                       // option's underlying
                                       const std::vector< std::vector<DateType> >&     startTimes,   // start dates of coupon period
                                       const std::vector< std::vector<DateType> >&     payTimes,     // pay dates of coupon perid
                                       const std::vector< std::vector<PassiveType> >&  cashFlows,    // fixed coupon payments (absolut value)
                                       const std::vector< Option::Type >&              cop,          // call (1) or put (-1) option
                                       // calibration parameters
                                       const PassiveType                tol_vola      // absolut tolerance in short rate volatility
                                );

        // calibrate short rate volatility and mean reversion s.t. minimal relative change
        // this method is uses method CalibrateVolatility(...)
        virtual
        const std::vector<ActiveType>& BermudanCalibration (  
                                       const std::vector<DateType>&     exercTimes,   // option's exercise dates (equal settlment)
                                       const std::vector<PassiveType>&  strikeValues, // strike payed at excercise dates
                                       const std::vector<PassiveType>&  b76Prices,    // reference European prices
                                       // option's underlying
                                       const std::vector< std::vector<DateType> >&     startTimes,   // start dates of coupon period
                                       const std::vector< std::vector<DateType> >&     payTimes,     // pay dates of coupon perid
                                       const std::vector< std::vector<PassiveType> >&  cashFlows,    // fixed coupon payments (absolut value)
                                       const std::vector< Option::Type >&              cop,          // call (1) or put (-1) option
                                       // calibration parameters
                                       const PassiveType                tol_vola      // absolut tolerance in short rate volatility
                                );
    };


    template <class DateType, class PassiveType, class ActiveType>
    HullWhiteModelT<DateType,PassiveType,ActiveType>::HullWhiteModelT(
                                const Handle<YieldTermStructure>& termStructure,
                                const PassiveType                mean,
                                const std::vector<DateType>&     volaTimes,
                                const std::vector<ActiveType>&   volaValues )
        : termStructure_(termStructure), mean_(mean), volaTimes_(volaTimes), volaValues_(volaValues) {
            calibrateMean_    = false; // default setting
            estimateAccuracy_ = true;  // default setting
            Nvol_ = std::min(volaTimes_.size(),volaValues_.size());
            europeansAnalytical_.resize(1);
            europeansNumerical_.resize(1);
            QL_REQUIRE (Nvol_ > 0, "no short rate volatilities provided.");
    }

    template <class DateType, class PassiveType, class ActiveType> ActiveType 
    HullWhiteModelT<DateType,PassiveType,ActiveType>::ZeroBond(
                            const DateType    settlement,
                            const DateType    maturity,
                            const ActiveType  shortRate,
                            ActiveType        *dBond_dr /*= 0*/ ) { 
        PassiveType fwRate = termStructure_->forwardRate(settlement,settlement,Continuous);	
        PassiveType B      = (1.0 - exp(-mean_ * (maturity-settlement))) / mean_;
        PassiveType DF1    = termStructure_->discount(settlement);
        PassiveType DF2    = termStructure_->discount(maturity);
        ActiveType  Integr = TemplateAuxilliaries::PieceWiseConstantIntegral<PassiveType,ActiveType,func_F>(volaTimes_, TemplateAuxilliaries::sqr(volaValues_), func_F(mean_,0.0))(0.0,settlement);
        ActiveType  A      = DF2/DF1*exp(B*fwRate - B*B/2.0*exp(-2.0*mean_*(settlement-0))*Integr);
        ActiveType  result = A * exp(-B*shortRate);
        if (dBond_dr) (*dBond_dr) = - B * result;
        return result;
    }

    template <class DateType, class PassiveType, class ActiveType> 
    ActiveType HullWhiteModelT<DateType,PassiveType,ActiveType>::CouponBond(
                            const DateType                  settlement,   // bond settlement date in year fractions
                            const std::vector<DateType>&    payTimes,     // pay dates in year fraction of coupon period
                            const std::vector<PassiveType>& cashFlows,    // fixed coupon payments (absolut value)
                            const ActiveType                shortRate,
                            ActiveType                      *dBond_dr,
                            const size_t                    idx_start) {
        size_t Ncfs = std::min(payTimes.size(),cashFlows.size());
        ActiveType PV=0, dPV=0;
        for (size_t i=idx_start; i<Ncfs; ++i) {
            PV += cashFlows[i]*ZeroBond(settlement,payTimes[i],shortRate,dBond_dr);
            if (dBond_dr) dPV += cashFlows[i]*(*dBond_dr);
        }
        if (dBond_dr) (*dBond_dr) = dPV;
        return PV;
    }

    template <class DateType, class PassiveType, class ActiveType> 
    ActiveType HullWhiteModelT<DateType,PassiveType,ActiveType>::ZeroBondOption(
                        const DateType     excercise,     // option's excercise date (equal settlment)
                        const ActiveType   strike,        // strike payed at excercise date
                        const DateType     maturity,      // payment date of notional 1
                        const Option::Type cop) {         // call (1) or put (-1) option
        ActiveType integral  = TemplateAuxilliaries::PieceWiseConstantIntegral<PassiveType,ActiveType,func_F>(volaTimes_, TemplateAuxilliaries::sqr(volaValues_), func_F(mean_,0.0))(0.0,excercise);
        ActiveType sigmaBond = (exp(-mean_*(excercise-0)) - exp(-mean_*(maturity-0)))/mean_ * sqrt(integral);
        return termStructure_->discount(excercise) * TemplateAuxilliaries::Black76<ActiveType>(
                    termStructure_->discount(maturity)/termStructure_->discount(excercise),
                    strike, sigmaBond, 1.0, cop);
    }

    template <class DateType, class PassiveType, class ActiveType> 
    ActiveType HullWhiteModelT<DateType,PassiveType,ActiveType>::CouponBondOption(
                        const DateType                  excercise,    // option's excercise date (equal settlment)
                        const PassiveType               strike,       // strike payed at excercise date
                        /* option's underlying, copy vectors */
                        const std::vector<DateType>&           startTimes,   // start dates of coupon period
                        const std::vector<DateType>&           payTimes,     // pay dates of coupon perid
                        const std::vector<PassiveType>&        cashFlows,    // fixed coupon payments (absolut value)
                        const Option::Type              cop){         // call (1) or put (-1) option
        // use Jamschdian's trick
        ActiveType rStar, cbStar, meritFunction, derivative, newtonStep, zeroStrike, res;
        double rTolerance = 1.0e-12; /* intervall Newton's method */
        size_t i, idx_start, max_iter = 50;
        size_t Ncfs = std::min(startTimes.size(),payTimes.size());
        Ncfs = std::min(Ncfs,cashFlows.size());

        if ((Ncfs<1)||(startTimes[Ncfs-1]<excercise)) { // no cashflows
            return 0; /* or discounted strike for put */
        }

        /* consider only coupons with start date after excercise date */
        idx_start = 0;
        while ((startTimes[idx_start]<excercise)&(idx_start<Ncfs-1)) ++idx_start;

        /* solve for rStar by Newton's method */
        rStar  = 0.0;
        for (i=0; i<max_iter; ++i) {
            cbStar = CouponBond( excercise, payTimes, cashFlows, rStar, &derivative, idx_start);
            meritFunction = cbStar - strike;
            newtonStep = - meritFunction / derivative;
            if (fabs(newtonStep)<rTolerance) break;
            rStar += newtonStep;
        }
        // if (fabs(newtonStep)>=rTolerance) report error
        QL_REQUIRE (fabs(newtonStep) < rTolerance, "can not solve for rStar in Jamschdian's decomposition.");
        // eliminate intermediates
        cbStar        = 0;
        meritFunction = 0;
        newtonStep    = 0;
        derivative    = 0;

        // evaluate the corresp. sum of ZCB options
        res = 0;
        for (i=idx_start; i<Ncfs; ++i) {
            zeroStrike = ZeroBond(excercise,payTimes[i],rStar);
            res += cashFlows[i] * ZeroBondOption(excercise,zeroStrike,payTimes[i],cop);
        }
        zeroStrike = 0;
        return res;
    }


    template <class DateType, class PassiveType, class ActiveType> 
    void HullWhiteModelT<DateType,PassiveType,ActiveType>::evaluateShortRateGrid( 
                           PassiveType r0,      // center of short rate grid
                           PassiveType s,       // distance to boundaries
                           size_t dim ) {       // number of grid points
        shortRateGrid_.resize(dim);
        PassiveType vol = s / TemplateAuxilliaries::PhiInv(dim /(dim+1.0));
        for (size_t i=0; i<dim; ++i) shortRateGrid_[i] = TemplateAuxilliaries::PhiInv((i+1.0)/(dim+1.0))*vol + r0;
        //shortRateGrid_ = TemplateAuxilliaries::Gau�TschebyschowIntegral().getGrid(r0-s, r0+s, dim);
    }

    template <class DateType, class PassiveType, class ActiveType> 
    ActiveType HullWhiteModelT<DateType,PassiveType,ActiveType>::BermudanBondOption(
                            const std::vector<DateType>&     exercTimes,    // option's exercise dates (equal settlment)
                            const std::vector<PassiveType>&  strikeValues,  // strike payed at excercise dates
                            // option's underlying
                            const std::vector<DateType>&     startTimes,   // start dates of coupon period
                            const std::vector<DateType>&     payTimes,     // pay dates of coupon perid
                            const std::vector<PassiveType>&  cashFlows,    // fixed coupon payments (absolut value)
                            const Option::Type               cop,          // call (1) or put (-1) option
                            // discretisation properties
                            const size_t                     dim,          // number of short rate grid points
                            const PassiveType                gridRadius,   // radius s of short rate grid [r0-s, r0+s]
                            const PassiveType                tol  ) {      // absolute tolerance for numerical integration
        DateType    startTime;		
        PassiveType r0, f0, f1;
        ActiveType  discretePV; //, variance, integral1, integral2, forwardDF, expectation;
        size_t i, j, k, idx_start, Nexc, Ncfs;
        
        std::vector< std::vector< std::vector< ActiveType > > > V, G, Z;
        std::vector< std::vector< ActiveType > > forwardDF, expectation;
        std::vector< ActiveType > variance, integral1, integral2;

        Ncfs = std::min(startTimes.size(),payTimes.size());
        Ncfs = std::min(Ncfs,cashFlows.size());
        QL_REQUIRE (Ncfs > 0, "no cashflows found.");

        r0 = termStructure_->forwardRate(0,0,Continuous);	
        evaluateShortRateGrid( r0, gridRadius, dim );

        Nexc = std::min(exercTimes.size(),strikeValues.size());
        V.resize(Nexc+1);      // V[0] Bermudan, V[1..Nexc] European numeric
        G.resize(Nexc+1);      // slope dV/dr for C2 interpolation
        Z.resize(Nexc+1);      // intermediates in linear solver
        V[0].resize(Nexc+1);   // Bermudan at each exercise + initiat time
        G[0].resize(Nexc+1);
        Z[0].resize(Nexc+1);
        for (i=0; i<Nexc+1; ++i) {
            V[0][i].resize(dim);
            G[0][i].resize(dim);
            Z[0][i].resize(dim);
        }
        for (k=1; k<Nexc+1; ++k) {
            V[k].resize(k+1);
            G[k].resize(k+1);
            Z[k].resize(k+1);
            for (i=0; i<k+1; ++i) {
                V[k][i].resize(dim);
                G[k][i].resize(dim);
                Z[k][i].resize(dim);
            }
        }
        forwardDF.resize(Nexc);
        expectation.resize(Nexc);
        variance.resize(Nexc);
        integral1.resize(Nexc);
        integral2.resize(Nexc);
        for (k=0; k<Nexc; ++k) {
            forwardDF[k].resize(dim);
            expectation[k].resize(dim);
        }

        for (i=0; i<dim; ++i) V[0][Nexc][i] = 0.0;  // ensure European equals Bermudan at last exercise

        for (long k=Nexc-1; k>=0; --k) {
            // consider only coupons with start date after excercise date
            idx_start = 0;
            while ((startTimes[idx_start]<exercTimes[k])&(idx_start<Ncfs-1)) ++idx_start;
            // evaluate pay-off at excercise
            for (i=0; i<dim; ++i) {
                V[k+1][k+1][i] = CouponBond( exercTimes[k], payTimes, cashFlows, shortRateGrid_[i],
                                             (ActiveType *) 0, idx_start );
                V[k+1][k+1][i] = (Integer(cop)*(V[k+1][k+1][i]-strikeValues[k])>0) ? Integer(cop)*(V[k+1][k+1][i]-strikeValues[k]) : (ActiveType)0.0;         
                /* evaluate Bermudan option */
                V[0][k+1][i] = (V[k+1][k+1][i]>V[0][k+1][i]) ? V[k+1][k+1][i] : V[0][k+1][i];
            }
            // evaluate ZCB() E^T [ V(r) ]
            // cubic spline interpolation of solution V only if required for integration (for non trivial tolerance)
            // otherwise use simple integration scheme
            if (1) {  // (tol>0)
                TemplateAuxilliaries::c2splineDerivatives(shortRateGrid_,V[0][k+1],G[0][k+1],Z[0][k+1]);
                if (estimateAccuracy_) {
                    for (j=k+1; j<Nexc+1; ++j) TemplateAuxilliaries::c2splineDerivatives(shortRateGrid_,V[j][k+1],G[j][k+1],Z[j][k+1]);
                }
            }
            {   // separate scope for temporary auxilliary objects
                // variance is independent of r
                startTime = (k>0) ? exercTimes[k-1] : 0;
                variance[k] = TemplateAuxilliaries::PieceWiseConstantIntegral<PassiveType,ActiveType,func_F>(volaTimes_, TemplateAuxilliaries::sqr(volaValues_), func_F(mean_,exercTimes[k]))(startTime,exercTimes[k]);
                // evaluate the nasty integrals
                integral1[k] = TemplateAuxilliaries::PieceWiseConstantIntegral<PassiveType,ActiveType,func_F>(volaTimes_, TemplateAuxilliaries::sqr(volaValues_), func_F(mean_,startTime    ))(0,startTime);
                integral1[k] *= exp(-mean_*(exercTimes[k]-startTime)) / mean_;
                integral2[k] = TemplateAuxilliaries::PieceWiseConstantIntegral<PassiveType,ActiveType,func_F>(volaTimes_, TemplateAuxilliaries::sqr(volaValues_), func_F(mean_,0            ))(0,startTime);
                integral2[k] *= exp(-2.0*mean_*(exercTimes[k]-0)) / mean_;
            }
            // forward correction in time-T neutral measure
            integral1[k] -= integral2[k];
            // forward rate f(0,t)	    
            f0 = termStructure_->forwardRate(startTime,startTime,Continuous);	
            // forward rate f(0,T)	    
            f1 = termStructure_->forwardRate(exercTimes[k],exercTimes[k],Continuous);	
            // now we can evaluate E[ r ] and integrate the payoff
            for (i=0; i<dim; ++i) {
                expectation[k][i] = f1 + exp(-mean_*(exercTimes[k]-startTime))*(shortRateGrid_[i] - f0);
                expectation[k][i] += integral1[k];
                // evaluate expectation for given payoff and density (given by variance and r-expectation)				
                V[0][k][i] = TemplateAuxilliaries::normalExpectation(shortRateGrid_, V[0][k+1], G[0][k+1], expectation[k][i], variance[k], tol);
                if (estimateAccuracy_) {
                    for (j=k+1; j<Nexc+1; ++j)
                        V[j][k][i] = TemplateAuxilliaries::normalExpectation(shortRateGrid_, V[j][k+1], G[j][k+1], expectation[k][i], variance[k], tol);
                }
            }
            for (i=0; i<dim; ++i) {
                forwardDF[k][i] = ZeroBond(startTime,exercTimes[k],shortRateGrid_[i],(ActiveType *) 0);
                V[0][k][i] *=  forwardDF[k][i];
                if (estimateAccuracy_) {
                    for (j=k+1; j<Nexc+1; ++j) V[j][k][i] *= forwardDF[k][i];
                }		
            }
        }
        // cubic spline interpolation of final solution V
        TemplateAuxilliaries::c2splineDerivatives(shortRateGrid_,V[0][0],G[0][0],Z[0][0]);
        if (estimateAccuracy_) {
            for (j=1; j<Nexc+1; ++j) TemplateAuxilliaries::c2splineDerivatives(shortRateGrid_,V[j][0],G[j][0],Z[j][0]);
        }
        // Bermudan solution
        discretePV = TemplateAuxilliaries::interpolCSpline( r0, shortRateGrid_, V[0][0], G[0][0] );
        // calculate analytic prices and compare numerical results
        europeansAnalytical_.resize(Nexc);
        europeansNumerical_.resize(Nexc);
        for (k=0; k<Nexc; ++k) {
            europeansNumerical_[k] = TemplateAuxilliaries::interpolCSpline( r0, shortRateGrid_, V[k+1][0], G[k+1][0] );
            europeansAnalytical_[k] = CouponBondOption( exercTimes[k], strikeValues[k], startTimes, payTimes, cashFlows, cop);
        }

        // Eliminate in reverse order
        // at t=0
        for (k=0; k<Nexc+1; ++k) {
            for (i=0; i<dim; ++i) V[k][0][i] = 0;
            for (i=0; i<dim; ++i) G[k][0][i] = 0;
            for (i=dim; i>0; --i) Z[k][0][i-1] = 0;
        }
        // at T_1 .. T_Nexc
        for (j=1; j<Nexc+1; ++j) {
            for (i=0; i<dim; ++i) V[0][j][i] = 0;
            for (k=j; k<Nexc+1; ++k) {
                for (i=0; i<dim; ++i) V[k][j][i] = 0;
            }
            for (i=0; i<dim; ++i) {
                forwardDF[j-1][i] = 0.0;
                expectation[j-1][i] = 0.0;
            }
            integral2[j-1] = 0.0;
            integral1[j-1] = 0.0;
            variance[j-1] = 0.0;
            for (i=0; i<dim; ++i) G[0][j][i] = 0;
            for (k=j; k<Nexc+1; ++k) {
                for (i=0; i<dim; ++i) G[k][j][i] = 0;
            }
            for (i=dim; i>0; --i) Z[0][j][i-1] = 0;
            for (k=j; k<Nexc+1; ++k) {
                for (i=dim; i>0; --i) Z[k][j][i-1] = 0;
            }
        }

        return discretePV;
    }

    template <class DateType, class PassiveType, class ActiveType> 
    const std::vector<ActiveType>& HullWhiteModelT<DateType,PassiveType,ActiveType>::CalibrateVolatility ( 
                               const std::vector<DateType>&     exercTimes,   // option's exercise dates (equal settlment)
                               const std::vector<PassiveType>&  strikeValues, // strike payed at excercise dates
                               const std::vector<PassiveType>&  b76Prices,    // reference European prices
                               // option's underlying
                               const std::vector< std::vector<DateType> >&     startTimes,   // start dates of coupon period
                               const std::vector< std::vector<DateType> >&     payTimes,     // pay dates of coupon perid
                               const std::vector< std::vector<PassiveType> >&  cashFlows,    // fixed coupon payments (absolut value)
                               const std::vector< Option::Type >&              cop,          // call (1) or put (-1) option
                               // calibration parameters
                               const PassiveType                tol_vola      // absolut tolerance in short rate volatility
                               ) {
        // solve successively for volaValues_ by secant method
        size_t k, Nexc;
        ActiveType sigma;      // iteration variable
        ActiveType sigma_min;  // lower bound for sigma
        ActiveType sigma_max;  // upper bound for sigma
        ActiveType sigma_last; // last iterate for secant evaluation
        ActiveType sigma_step; // secant step

        ActiveType price;      // corresponding Hull White prices
        ActiveType price_min;
        ActiveType price_max;
        ActiveType price_last;

        Nexc = std::min(exercTimes.size(),strikeValues.size());
        Nexc = std::min(Nexc,b76Prices.size());
        QL_REQUIRE (Nexc > 0, "no exercises found.");

        volaTimes_.resize(Nexc);
        volaValues_.resize(Nexc);
        calibIters_.resize(Nexc);
        bermudanB76Prices_.resize(Nexc);
        Nvol_ = Nexc;
        for (k=0; k<Nexc; ++k) volaTimes_[k]         = exercTimes[k];
        for (k=0; k<Nexc; ++k) volaValues_[k]        = 0.0; // initialisation
        for (k=0; k<Nexc; ++k) bermudanB76Prices_[k] = b76Prices[k]; 

        for (k=0; k<Nexc; ++k) {
            // initialisation
            calibIters_[k] = 0;
            sigma_min      = 1.0e-4;
            sigma_max      = 1.0e-2;
            sigma_last     = sigma_min;
            sigma          = sigma_max;
            // boundaries for prices
            volaValues_[k] = sigma_min;
            price_min      = CouponBondOption( exercTimes[k], strikeValues[k], startTimes[std::min(k,startTimes.size()-1)], payTimes[std::min(k,payTimes.size()-1)], cashFlows[std::min(k,cashFlows.size()-1)], cop[std::min(k,cop.size()-1)]);
            // QL_REQUIRE (price_min <= b76Prices[k], "Black'76 price too small.");
            if (price_min > b76Prices[k]) continue;  // we relax the calibration condition
            volaValues_[k] = sigma_max;
            price_max      = CouponBondOption( exercTimes[k], strikeValues[k], startTimes[std::min(k,startTimes.size()-1)], payTimes[std::min(k,payTimes.size()-1)], cashFlows[std::min(k,cashFlows.size()-1)], cop[std::min(k,cop.size()-1)]);
            while ((price_max < b76Prices[k])&&(sigma_max<1.0)) {
                sigma_max *= 2.0;
                volaValues_[k] = sigma_max;
                price_max      = CouponBondOption( exercTimes[k], strikeValues[k], startTimes[std::min(k,startTimes.size()-1)], payTimes[std::min(k,payTimes.size()-1)], cashFlows[std::min(k,cashFlows.size()-1)], cop[std::min(k,cop.size()-1)]);
            }
            QL_REQUIRE (price_max >= b76Prices[k], "Black'76 price too big.");
            // initialisations for iteration
            price_last     = price_min;
            price          = price_max;
            while (fabs(sigma-sigma_last)>tol_vola) {
                sigma_step  = - (price - b76Prices[k])*(sigma-sigma_last)/(price-price_last);
                sigma_last  = sigma;
                price_last  = price;
                sigma       = sigma_last + sigma_step;
                // ensure bounded iteration
                if ((sigma<sigma_min)||(sigma>sigma_max)) sigma = (sigma_min + sigma_max)/2;
                volaValues_[k] = sigma;
                price          = CouponBondOption( exercTimes[k], strikeValues[k], startTimes[std::min(k,startTimes.size()-1)], payTimes[std::min(k,payTimes.size()-1)], cashFlows[std::min(k,cashFlows.size()-1)], cop[std::min(k,cop.size()-1)]);
                // new boundaries
                if (price<=b76Prices[k]){
                    sigma_min = sigma;
                    price_min = price;
                }
                if (price>=b76Prices[k]){
                    sigma_max = sigma;
                    price_max = price;
                }
                calibIters_[k] += 1;
            }
        }
        return volaValues_;
    }

    template <class DateType, class PassiveType, class ActiveType> 
    const std::vector<ActiveType>& HullWhiteModelT<DateType,PassiveType,ActiveType>::BermudanCalibration ( 
                               const std::vector<DateType>&     exercTimes,   // option's exercise dates (equal settlment)
                               const std::vector<PassiveType>&  strikeValues, // strike payed at excercise dates
                               const std::vector<PassiveType>&  b76Prices,    // reference European prices
                               // option's underlying
                               const std::vector< std::vector<DateType> >&     startTimes,   // start dates of coupon period
                               const std::vector< std::vector<DateType> >&     payTimes,     // pay dates of coupon perid
                               const std::vector< std::vector<PassiveType> >&  cashFlows,    // fixed coupon payments (absolut value)
                               const std::vector< Option::Type >&              cop,          // call (1) or put (-1) option
                               // calibration parameters
                               const PassiveType                tol_vola      // absolut tolerance in short rate volatility
                               ) {
        // in case mean reversion should not be calibrated simply pass call to CalibrateVolatility(...)
        CalibrateVolatility(exercTimes,strikeValues,b76Prices,startTimes,payTimes,cashFlows,cop,tol_vola);
        if ((!calibrateMean_) | (volaValues_.size()<2)) return volaValues_;  // we require at least two vola values for optimization
        // solve min_a {0.5 ||F(sigma)||^2} s.t. CBO(sigma,a) = Swaption with F(sigma) = C log(sigma),
        // C = tridiag{ [0..0], [1..1, 0], [-1..-1] }, i.e. minimum relative change
        std::vector<ActiveType> vola0, vola1, dvola, F, dF;
        PassiveType mean0, mean1, mean2, lambda, mean_min, mean_max;
        PassiveType tol_mean;
        ActiveType  tmp1, tmp2;
        Size        iter=0, max_iter;
        // initial function evaluations
        tol_mean = sqrt(tol_vola);  // this is probably a bit sloppy
        max_iter = 10;
        mean_min = 1.0e-4;
        mean_max = 1.0;
        mean0    = mean_ + 1.0; // force recalculation
        mean1    = mean_;
        dvola.resize(volaValues_.size());
        F.resize(volaValues_.size());
        dF.resize(volaValues_.size());
        vola0 = volaValues_;
        while ((fabs(mean1 - mean0)>=tol_mean)&(iter<max_iter)) {
            ++iter;
            // for large steps the derivative approximation does not hold and we need an additional function evaluation
            if (fabs(mean1 - mean0)>1.0e-2) {
                mean0 = mean1 + 1.0e-2;
                mean_ = mean0;
                CalibrateVolatility(exercTimes,strikeValues,b76Prices,startTimes,payTimes,cashFlows,cop,tol_vola);
                vola0 = volaValues_;
            }
            mean_ = mean1;
            CalibrateVolatility(exercTimes,strikeValues,b76Prices,startTimes,payTimes,cashFlows,cop,tol_vola);
            vola1 = volaValues_;
            for (Size k=0; k<volaValues_.size(); ++k) dvola[k] = (vola1[k] - vola0[k])/(mean1 - mean0);
            for (Size k=0; k<volaValues_.size(); ++k) {
                F[k]  = log(vola1[k]);
                dF[k] = dvola[k] / vola1[k];
            }
            for (Size k=0; k<volaValues_.size()-1; ++k) {
                F[k]  -= F[k+1];
                dF[k] -= dF[k+1];
            }
            tmp1 = 0;
            tmp2 = 0;
            for (Size k=0; k<volaValues_.size()-1; ++k) {
                tmp1 += F[k]*dF[k];
                tmp2 += dF[k]*dF[k];
            }
            lambda = TemplateAuxilliaries::DBL( tmp1 / tmp2 );  // check for division by zero...
            // mean2 is the new iterate s.t. box constraints
            mean2 = mean1 - lambda;
            if (mean2 < mean_min) mean2 = mean_min;
            if (mean2 > mean_max) mean2 = mean_max;
            // iterate state variables
            mean0  = mean1;
            vola0  = vola1;
            mean1  = mean2;
        }
        return volaValues_;
    }

}

#endif  /* ifndef quantlib_templatehullwhitemodel_hpp */
