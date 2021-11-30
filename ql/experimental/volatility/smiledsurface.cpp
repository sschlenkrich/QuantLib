/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2016 Sebastian Schlenkrich

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

/*! \file smiledsurface.cpp
    \brief BlackVolTermStructure based on interpolated SmileSections
*/

#include <ql/experimental/volatility/smiledsurface.hpp>
#include <ql/math/interpolations/cubicinterpolation.hpp>


namespace QuantLib {

    Rate SmiledSurface::minStrike() const {
        QL_REQUIRE(smiles_.size()>0, "SmiledSurface error: no smiles provided");
        Rate mStrike = smiles_[0]->minStrike();
        for (Size k=1; k<smiles_.size(); ++k) if (smiles_[k]->minStrike()<mStrike) mStrike = smiles_[k]->minStrike();
        return mStrike;
    }

    Rate SmiledSurface::maxStrike() const {
        QL_REQUIRE(smiles_.size()>0, "SmiledSurface error: no smiles provided");
        Rate mStrike = smiles_[0]->maxStrike();
        for (Size k=1; k<smiles_.size(); ++k) if (smiles_[k]->maxStrike()>mStrike) mStrike = smiles_[k]->maxStrike();
        return mStrike;
    }

    Date SmiledSurface::maxDate() const {
        QL_REQUIRE(smiles_.size()>0, "SmiledSurface error: no smiles provided");
        Date mDate = smiles_[0]->exerciseDate();
        for (Size k=1; k<smiles_.size(); ++k) if (smiles_[k]->exerciseDate()>mDate) mDate = smiles_[k]->exerciseDate();
        return mDate;
    }

    Real SmiledSurface::blackVarianceImpl(Time t, Real strike) const {
        QL_REQUIRE(smiles_.size()>0, "SmiledSurface error: no smiles provided");
        std::vector<Time>       times(smiles_.size());
        std::vector<Real>       vars(smiles_.size());
        for (Size k=0; k<smiles_.size(); ++k) {
            times[k] = smiles_[k]->exerciseTime();       // assume consistent time calculation
            vars[k]  = smiles_[k]->variance(strike);     // assume no volatility type conversion
        }
        MonotonicCubicNaturalSpline interp(times.begin(),times.end(),vars.begin());  // assume ascending times
        Real unSafeVariance = interp(t, true);
        Real minVariance = 1.0e-8 * t;  // assume 1bp minimum volatility
        Real safeVariance = (unSafeVariance < minVariance) ? (minVariance) : (unSafeVariance);
        return safeVariance;
    }

    Volatility SmiledSurface::blackVolImpl(Time t, Real strike) const {
        QL_REQUIRE(t>0, "SmiledSurface error: positive time required");
        return sqrt(blackVarianceImpl(t,strike) / t);
    }
}

