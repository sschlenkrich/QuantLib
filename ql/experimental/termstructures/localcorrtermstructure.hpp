/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2017 Cord Harms

 This file is part of QuantLib, a free-software/open-source library
 for financial quantitative analysts and developers - http://quantlib.org/

 QuantLib is free software: you can redistribute it and/or modify it
 under the terms of the QuantLib license.  You should have received a
 copy of the license along with this program; if not, please email
 <quantlib-dev@lists.sf.net>. The license is also available online at
 <http://quantlib.org/license.shtml>.

 This program is distributed in the hope that it will be useful, but WITHOUT
 ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 FOR A PARTICULAR PURPOSE.  See the license for more details.
*/

/*! \file localcorrtermstructure.hpp
    \brief Local correlation term structure base class
*/

#ifndef quantlib_local_corr_term_structures_hpp
#define quantlib_local_corr_term_structures_hpp

#include <ql/experimental/termstructures/corrtermstructureStrike.hpp>
#include <ql/patterns/visitor.hpp>
#include <ql/experimental/templatemodels/stochasticprocessT.hpp>
#include <ql/experimental/processes/hestonslvprocess.hpp>
#include <ql/processes/blackscholesprocess.hpp>

namespace QuantLib {

    /*! This abstract class defines the interface of concrete
        local-correlation term structures which will be derived from this one.

        Correlations are assumed to be expressed on an annual basis.
    */
    class LocalCorrTermStructure : public CorrelationTermStructureStrike {
      public:
        /*! \name Constructors
            See the TermStructure documentation for issues regarding
            constructors.
        */
        //@{
        //! default constructor
        /*! \warning term structures initialized by means of this
                     constructor must manage their own reference date
                     by overriding the referenceDate() method.
        */
        LocalCorrTermStructure(const std::vector<ext::shared_ptr<QuantLib::GeneralizedBlackScholesProcess>>& processes,
                               const ext::shared_ptr<QuantLib::GeneralizedBlackScholesProcess>&			   processToCal );
        LocalCorrTermStructure(const std::vector<ext::shared_ptr<QuantLib::HestonSLVProcess>>& processes,
            const ext::shared_ptr<QuantLib::GeneralizedBlackScholesProcess>&			   processToCal);

        //! initialize with a fixed reference date
        //@}
        virtual ~LocalCorrTermStructure() = default;
        //! \name Local correlation
        //@{
        void localCorr(RealStochasticProcess::MatA& correlationMatrix, 
                            const Date& d,
                            const RealStochasticProcess::VecA& X0,
                            bool extrapolate = false);
        void localCorr(RealStochasticProcess::MatA& correlationMatrix,
                            Time t,
                            const RealStochasticProcess::VecA& X0, 
                            bool extrapolate = false);
        //@}
        //! \name Visitability
        //@{
        virtual void accept(AcyclicVisitor&);
        //@}
        //! \name TermStructure interface
        //@{
        const Date& referenceDate() const;
        DayCounter dayCounter() const;
        Date maxDate() const;
        //@}
        //! \name CorrelationTermStructure interface
        //@{
        Real minStrike(Natural ulId) const;
        Real maxStrike(Natural ulId) const;


        std::vector<ext::shared_ptr<QuantLib::HestonSLVProcess>>& getProcesses() { return processes_; };
        ext::shared_ptr<QuantLib::GeneralizedBlackScholesProcess>& getProcessToCal() { return processToCal_; };

      protected:
        /*! \name Calculations

            These methods must be implemented in derived classes to perform
            the actual correlation calculations. When they are called,
            range check has already been performed; therefore, they must
            assume that extrapolation is required.
        */
        //@{
        //! local corr calculation
        virtual void localCorrImpl(RealStochasticProcess::MatA& corrMatrix, Time t, const RealStochasticProcess::VecA& X0,
            bool extrapolate = false) = 0;
        //@}
        std::vector<ext::shared_ptr<QuantLib::HestonSLVProcess>>				processes_;
        ext::shared_ptr<QuantLib::GeneralizedBlackScholesProcess>			    processToCal_;
    };

}

#endif
