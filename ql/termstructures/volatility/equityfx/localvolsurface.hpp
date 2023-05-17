/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2003 Ferdinando Ametrano

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

/*! \file localvolsurface.hpp
    \brief Local volatility surface derived from a Black vol surface
*/

#ifndef quantlib_localvolsurface_hpp
#define quantlib_localvolsurface_hpp

#include <ql/termstructures/volatility/equityfx/localvoltermstructure.hpp>
#include <ql/termstructures/volatility/equityfx/fixedlocalvolsurface.hpp>

namespace QuantLib {

    class BlackVolTermStructure;
    class YieldTermStructure;
    class Quote;

    //! Local volatility surface derived from a Black vol surface
    /*! For details about this implementation refer to
        "Stochastic Volatility and Local Volatility," in
        "Case Studies and Financial Modelling Course Notes," by
        Jim Gatheral, Fall Term, 2003

        see www.math.nyu.edu/fellows_fin_math/gatheral/Lecture1_Fall02.pdf

        \bug this class is untested, probably unreliable.
    */
    class LocalVolSurface : public LocalVolTermStructure {
      public:
        LocalVolSurface(const Handle<BlackVolTermStructure>& blackTS,
                        Handle<YieldTermStructure> riskFreeTS,
                        Handle<YieldTermStructure> dividendTS,
                        Handle<Quote> underlying);
        LocalVolSurface(const Handle<BlackVolTermStructure>& blackTS,
                        Handle<YieldTermStructure> riskFreeTS,
                        Handle<YieldTermStructure> dividendTS,
                        Real underlying);
        //! \name TermStructure interface
        //@{
        const Date& referenceDate() const override;
        DayCounter dayCounter() const override;
        Date maxDate() const override;
        //@}
        //! \name VolatilityTermStructure interface
        //@{
        Real minStrike() const override;
        Real maxStrike() const override;
        //@}
        //! \name Visitability
        //@{
        void accept(AcyclicVisitor&) override;
        //@}

        inline Handle<Quote> getUnderlying() { return underlying_; };
        inline Handle<YieldTermStructure>& getDividendTS() { return dividendTS_; };
        inline Handle<YieldTermStructure>& getInterestRateTS() { return riskFreeTS_; };
        inline Handle<BlackVolTermStructure>& getBlackSurface() { return blackTS_; };
      protected:
        Volatility localVolImpl(Time, Real) const override;

      private:
        Handle<BlackVolTermStructure> blackTS_;
        Handle<YieldTermStructure> riskFreeTS_, dividendTS_;
        Handle<Quote> underlying_;
    };


    class InterpolatedLocalVolSurface : public LocalVolSurface {
    public:
        InterpolatedLocalVolSurface(const Handle<BlackVolTermStructure>& blackTS,
            const Handle<YieldTermStructure>& riskFreeTS,
            const Handle<YieldTermStructure>& dividendTS,
            const Handle<Quote>& underlying, Size strikeGridAmt,
            Size timeStepsPerYear);

        inline Matrix getSurface() {
            Matrix mat(gridTimes_.size(), strikes_[0]->size());

            for (size_t i = 0; i < mat.rows(); i++)
            {
                for (size_t j = 0; j < mat.columns(); j++)
                {
                    mat[i][j] = surface_->localVol(gridTimes_[i], strikes_[i]->at(j), true);
                }
            }

            return mat;
        };

    protected:
          Volatility localVolImpl(Time, Real) const override;

        private:
        ext::shared_ptr<FixedLocalVolSurface> surface_;
        std::vector<Time> gridTimes_;
        std::vector<ext::shared_ptr<std::vector<Real> > > strikes_;
    };

}

#endif
