/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2016 Sebastian Schlenkrich

*/

/*! \file mcscriptT.hpp
    \brief script payoffs for MC simulation
    
*/


#ifndef quantlib_templatemcscript_hpp
#define quantlib_templatemcscript_hpp



#include <ql/experimental/templatemodels/montecarlo/mcpayoffT.hpp>



#include <ql/experimental/templatemodels/montecarlo/scripting/flexbisondriver.hpp>
#include <ql/experimental/templatemodels/montecarlo/scripting/expression.hpp>

#include <ql/time/date.hpp>

#include <boost/regex.hpp>



namespace QuantLib {

    template <class DateType, class PassiveType, class ActiveType>
    class MCScriptT : public MCPayoffT<DateType, PassiveType, ActiveType> {
    protected:
        // easy use of templated types
        typedef MCSimulationT<DateType, PassiveType, ActiveType>        SimulationType;
        typedef MCPayoffT<DateType, PassiveType, ActiveType>            PayoffType;
        typedef BasePayoffT<DateType, PassiveType, ActiveType>          MCBase;
        typedef typename MCSimulationT<DateType, PassiveType, ActiveType>::Path  PathType;

    private:
        std::map<std::string, ext::shared_ptr<PayoffType>>   payoffs_;      // the actual payoffs which may be accessed
        std::vector<std::string>                             expressions_;  // resulting expressions after parsing the script but before syntactic analysis
        std::vector<std::string>                             scriptLog_;    // log messages when parsing the script
        ext::shared_ptr<PayoffType>                          result_;       // result payoff for MCPayoff interface implementation
    public:
        MCScriptT(const std::vector<std::string>&                    keys,
                  const std::vector <ext::shared_ptr<PayoffType>>&   payoffs,
                  const std::vector<std::string>&                    script,
                  const bool                                         overwrite=true) : MCPayoffT<DateType, PassiveType, ActiveType>(0.0) {
            QL_REQUIRE(keys.size()==payoffs.size(), "MCScript error: key vs. value size missmatch");
            for (Size k = 0; k < keys.size(); ++k) { // initialize map
                typename std::map< std::string, ext::shared_ptr<PayoffType> >::iterator it = payoffs_.find(keys[k]);
                if (it == payoffs_.end()) { // insert a new element
                    payoffs_.insert( std::make_pair(keys[k], payoffs[k]) );
                } else { // potentially overwrite existing element
                    QL_REQUIRE(overwrite, "MCScript error: overwrite not allowed");
                    it->second = payoffs[k];
                }
            }
            if ((script.size() > 0) && (script[0].compare("NonRecursive") == 0)) {
                //parseScript(script, overwrite);           // deprecated and for debugging purpose
            } else {
                parseFlexBisonScript(script, overwrite);  // for briefty we delegate parsing to separate method
            }
            QL_REQUIRE(payoffs_.size()>0, "MCScript error: no payoffs stored.");
            result_ = payoffs_.rbegin()->second; // pick the last element as fall back

            // we need to find a 'result' payoff
            typename std::map< std::string, ext::shared_ptr<PayoffType> >::iterator it = payoffs_.find("payoff");
            if (it != payoffs_.end()) result_ = it->second;
            else result_ = payoffs_.rbegin()->second; // pick the last element as fall back
            PayoffType::observationTime_ = result_->observationTime();
        }

        inline virtual ActiveType at(const ext::shared_ptr<PathType>& p) {
            return result_->at(p);
        }

        // inspector
        inline const std::map<std::string, ext::shared_ptr<PayoffType>>&  payoffs()     { return payoffs_;   }
        inline const std::vector<std::string>&                            expressions() { return expressions_; }
        inline const std::vector<std::string>&                            scriptLog()   { return scriptLog_; }

        // return all keys
        std::vector<std::string> payoffsKeys() {
            std::vector<std::string> keyVector;
            for (typename std::map<std::string, ext::shared_ptr<PayoffType>>::iterator it = payoffs_.begin(); it != payoffs_.end(); ++it) {
                keyVector.push_back(it->first);
            }
            return keyVector;
        }

        // return all payoffs (values in map)
        std::vector < ext::shared_ptr<PayoffType> > payoffValues() {
            std::vector < ext::shared_ptr<PayoffType> > payoffVector;
            for (typename std::map<std::string, ext::shared_ptr<PayoffType> >::iterator it = payoffs_.begin(); it != payoffs_.end(); ++it) {
                payoffVector.push_back(it->second);
            }
            return payoffVector;
        }

        std::vector<DateType> observationTimes(const std::vector<std::string>& keys) {
            std::vector<ext::shared_ptr<PayoffType>> payoffs = findPayoffs(keys);
            std::set<DateType> s;
            for (size_t k = 0; k < payoffs.size(); ++k) s = PayoffType::unionTimes(s, payoffs[k]->observationTimes());
            return std::vector<DateType>(s.begin(), s.end());
        }

        // MC valuation
        inline std::vector<ActiveType> NPV(const ext::shared_ptr<SimulationType>&    simulation,
                                           const std::vector<std::string>&           keys) {
            std::vector<ext::shared_ptr<PayoffType>> payoffs = findPayoffs(keys);
            std::vector<ActiveType> npv(payoffs.size(), (ActiveType)0.0);
            for (Size n = 0; n < simulation->nPaths(); ++n) {
                const ext::shared_ptr<PathType> p(simulation->path(n));
                for (Size k = 0; k < payoffs.size(); ++k) {
                    npv[k] += payoffs[k]->discountedAt(p);
                }
            }
            for (Size k = 0; k < payoffs.size(); ++k) npv[k] /= simulation->nPaths();
            return npv;
        }

        // some helper functions to simplify Asset payoff handling
        inline static bool add_FixingTimes_to_Asset(
            ext::shared_ptr<PayoffType>        payoff,
            const std::vector<DateType>&       fixingTimes,
            const std::vector<PassiveType>&    fixingValues) {
            ext::shared_ptr<typename MCBase::Asset>  assetPayoff =
                ext::dynamic_pointer_cast<typename MCBase::Asset>(payoff);
            QL_REQUIRE(assetPayoff, "Payoff is no Asset");
            QL_REQUIRE(fixingTimes.size() == fixingValues.size(), "fixingTimes.size()==fixingValues.size() required");
            std::vector< std::pair<DateType, PassiveType> > history;
            for (size_t k = 0; k < fixingTimes.size(); ++k)
                history.push_back(std::make_pair(fixingTimes[k], fixingValues[k]));
            assetPayoff->addFixings(history);
            return true;
        }

        inline static bool add_FixingDates_to_Asset(
            ext::shared_ptr<PayoffType>        payoff,
            const std::vector<Date>&           fixingDates,
            const std::vector<PassiveType>&    fixingValues) {
            std::vector<DateType> fixingTimes(fixingDates.size());
            Date today = Settings::instance().evaluationDate();
            for (size_t k = 0; k < fixingDates.size(); ++k)
                fixingTimes[k] = (DateType)(((fixingDates[k].serialNumber() - today.serialNumber()) / 365.0));			
            return add_FixingTimes_to_Asset(payoff, fixingTimes, fixingValues);
        }

    private:

        // convert string to number
        inline static bool to_Number(const std::string str, ActiveType& number) {
            double res;
            std::string::size_type sz;
            try {
                res = std::stod(str, &sz);
            }
            catch (std::exception e) {
                return false;
            }
            number = res;
            return true;
        }

        // convert a Date to number
        inline static DateType date_to_Number(const Date& d) {
            Date today = Settings::instance().evaluationDate();
            DateType number = (DateType)(((d.serialNumber() - today.serialNumber()) / 365.0));
            return number;
        }

        // convert date string with format ddmmmyyyy to number
        inline static bool date_to_Number(const std::string str, ActiveType& number) {
            if (str.length() != 9) return false;
            std::string::size_type sz;
            try {
                Day  day  = std::stol(str.substr(0,2), &sz);
                Year year = std::stol(str.substr(5,4), &sz);
                Month month;
                std::string s = str.substr(2, 3);
                if (s.compare("Jan") == 0) month = Month::Jan;
                else if (s.compare("Feb") == 0) month = Month::Feb;
                else if (s.compare("Mar") == 0) month = Month::Mar;
                else if (s.compare("Apr") == 0) month = Month::Apr;
                else if (s.compare("May") == 0) month = Month::May;
                else if (s.compare("Jun") == 0) month = Month::Jun;
                else if (s.compare("Jul") == 0) month = Month::Jul;
                else if (s.compare("Aug") == 0) month = Month::Aug;
                else if (s.compare("Sep") == 0) month = Month::Sep;
                else if (s.compare("Oct") == 0) month = Month::Oct;
                else if (s.compare("Nov") == 0) month = Month::Nov;
                else if (s.compare("Dec") == 0) month = Month::Dec;
                else return false;
                Date d(day, month, year);
                Date today = Settings::instance().evaluationDate();
                number = (ActiveType)(((d.serialNumber() - today.serialNumber()) / 365.0));
            }
            catch (std::exception e) {
                return false;
            }
            return true;
        }



        // we define that helper function to simplify code in forthcoming expression parsing
        inline bool hasChilds(const ext::shared_ptr<Scripting::Expression> tree, Size nArgs, Size lineNr) {
            // make sure we can actually do something with the tree
            if (!tree) {
                scriptLog_.push_back(std::string("Error line " + std::to_string(lineNr) + ": Empty expression tree."));
                return false;
            }
            if (tree->childs().size() != nArgs) {
                scriptLog_.push_back(std::string("Error line " + std::to_string(lineNr) + ": " + std::to_string(nArgs) + " child expressions expected, but " + std::to_string(tree->childs().size()) + " found." ));
                return false;
            }
            return true;
        }

        // we define that helper function to simplify code in forthcoming expression parsing
        inline bool hasChildsInRange(const ext::shared_ptr<Scripting::Expression> tree, Size nArgsMin, Size nArgsMax, Size lineNr) {
            // make sure we can actually do something with the tree
            if (!tree) {
                scriptLog_.push_back(std::string("Error line " + std::to_string(lineNr) + ": Empty expression tree."));
                return false;
            }
            if ((tree->childs().size() < nArgsMin)||(tree->childs().size() > nArgsMax)) {
                scriptLog_.push_back(std::string("Error line " + std::to_string(lineNr) + ": [" + std::to_string(nArgsMin) + ", " + std::to_string(nArgsMax) +
                    "] child expressions expected, but " + std::to_string(tree->childs().size()) + " found."));
                return false;
            }
            return true;
        }

        // we define that helper function to simplify code in forthcoming expression parsing
        inline bool hasLeafs(const ext::shared_ptr<Scripting::Expression> tree, Size nArgs, Size lineNr) {
            // make sure we can actually do something with the tree
            if (!tree) {
                scriptLog_.push_back(std::string("Error line " + std::to_string(lineNr) + ": Empty expression tree."));
                return false;
            }
            if (tree->leafs().size() != nArgs) {
                scriptLog_.push_back(std::string("Error line " + std::to_string(lineNr) + ": " + std::to_string(nArgs) + " leafs expected, but " + std::to_string(tree->leafs().size()) + " found."));
                return false;
            }
            return true;
        }

        // we define that helper function to simplify code in forthcoming expression parsing
        inline bool hasLeafsInRange(const ext::shared_ptr<Scripting::Expression> tree, Size nArgsMin, Size nArgsMax, Size lineNr) {
            // make sure we can actually do something with the tree
            if (!tree) {
                scriptLog_.push_back(std::string("Error line " + std::to_string(lineNr) + ": Empty expression tree."));
                return false;
            }
            if ((tree->leafs().size()<nArgsMin)|| (tree->leafs().size()>nArgsMax)) {
                scriptLog_.push_back(std::string("Error line " + std::to_string(lineNr) + ": [" + std::to_string(nArgsMin) + ", " + std::to_string(nArgsMax) +
                    "] leafs expected, but " + std::to_string(tree->leafs().size()) + " found."));
                return false;
            }
            return true;
        }

        // check if a list of payoffs exists before doing some computationally expensive stuff with them
        std::vector< ext::shared_ptr<PayoffType> > findPayoffs(const std::vector<std::string>& keys, bool throwException = true) {
            std::vector<ext::shared_ptr<PayoffType>> payoffs;
            for (Size k = 0; k < keys.size(); ++k) {
                typename std::map<std::string, ext::shared_ptr<PayoffType> >::iterator it = payoffs_.find(keys[k]);
                if (it != payoffs_.end()) {
                    payoffs.push_back(it->second);
                    continue; // all done for this key
                }
                if (throwException) QL_FAIL("MCScript error: payoff '" + keys[k] + "' not found");
            }
            return payoffs;
        }

        // convert an abstract expression tree into a payoff
        // this function does the actual work...
        ext::shared_ptr<PayoffType> payoff(const ext::shared_ptr<Scripting::Expression> tree, Size k) {
            // make sure we can actually do something with the tree
            if (!tree) {
                scriptLog_.push_back(std::string("Error line " + std::to_string(k) + ": Empty expression tree."));
                QL_FAIL("Cannot interprete payoff");
            }
            // check any possible expression
            switch (tree->type()) {
            // expressions based on tokens
            case Scripting::Expression::NUMBER: {
                if (!hasChilds(tree, 0, k)) QL_FAIL("Cannot interprete payoff");
                if (!hasLeafs(tree, 1, k)) 	QL_FAIL("Cannot interprete payoff");
                ActiveType number;
                if (to_Number(tree->leafs()[0], number)) {
                    return ext::shared_ptr<PayoffType>(new typename MCBase::FixedAmount(number));
                }
                scriptLog_.push_back(std::string("Error line " + std::to_string(k) + ": cannot convert " + tree->leafs()[0] + " to number."));
                QL_FAIL("Cannot interprete payoff");
            }
            case Scripting::Expression::IDENTIFIER: {
                if (!hasChilds(tree, 0, k)) QL_FAIL("Cannot interprete payoff");
                if (!hasLeafs(tree, 1, k))  QL_FAIL("Cannot interprete payoff");
                // check for existing payoff in map
                typename std::map<std::string, ext::shared_ptr<PayoffType>>::iterator it = payoffs_.find(tree->leafs()[0]);
                if (it != payoffs_.end()) {
                    scriptLog_.push_back(std::string("Payoff line " + std::to_string(k) + ": '" + tree->leafs()[0] + "' is in map"));
                    return it->second;
                }
                // if we end up here no conversion was successfull
                scriptLog_.push_back(std::string("Error line " + std::to_string(k) + ": '" + tree->leafs()[0] + "' is no payoff"));
                QL_FAIL("Cannot interprete payoff");
            }
            // expressions basen on unary operators
            case Scripting::Expression::UNARYPLUS: {
                if (!hasChilds(tree, 1, k)) QL_FAIL("Cannot interprete payoff");
                if (!hasLeafs(tree, 0, k))  QL_FAIL("Cannot interprete payoff");
                return this->payoff(tree->childs()[0], k);
            }
            case Scripting::Expression::UNARYMINUS: {
                if (!hasChilds(tree, 1, k)) QL_FAIL("Cannot interprete payoff");
                if (!hasLeafs(tree, 0, k))  QL_FAIL("Cannot interprete payoff");
                return ext::shared_ptr<PayoffType>(new typename MCBase::Axpy(-1.0, payoff(tree->childs()[0], k), 0));
            }
            case Scripting::Expression::PLUS: {
                if (!hasChilds(tree, 2, k)) QL_FAIL("Cannot interprete payoff");
                if (!hasLeafs(tree, 0, k))  QL_FAIL("Cannot interprete payoff");
                return ext::shared_ptr<PayoffType>(new typename MCBase::Axpy(1.0, payoff(tree->childs()[0], k), payoff(tree->childs()[1], k)));
            }
            case Scripting::Expression::MINUS: {
                if (!hasChilds(tree, 2, k)) QL_FAIL("Cannot interprete payoff");
                if (!hasLeafs(tree, 0, k))  QL_FAIL("Cannot interprete payoff");
                return ext::shared_ptr<PayoffType>(new typename MCBase::Axpy(-1.0, payoff(tree->childs()[1], k), payoff(tree->childs()[0], k)));
            }
            case Scripting::Expression::MULT: {
                if (!hasChilds(tree, 2, k)) QL_FAIL("Cannot interprete payoff");
                if (!hasLeafs(tree, 0, k))  QL_FAIL("Cannot interprete payoff");
                return ext::shared_ptr<PayoffType>(new typename MCBase::Mult(payoff(tree->childs()[0], k), payoff(tree->childs()[1], k)));
            }
            case Scripting::Expression::DIVISION: {
                if (!hasChilds(tree, 2, k)) QL_FAIL("Cannot interprete payoff");
                if (!hasLeafs(tree, 0, k))  QL_FAIL("Cannot interprete payoff");
                return ext::shared_ptr<PayoffType>(new typename MCBase::Division(payoff(tree->childs()[0], k), payoff(tree->childs()[1], k)));
            }
            case Scripting::Expression::IFTHENELSE: {
                if (!hasChilds(tree, 3, k)) QL_FAIL("Cannot interprete payoff");
                if (!hasLeafs(tree, 0, k))  QL_FAIL("Cannot interprete payoff");
                return ext::shared_ptr<PayoffType>(new typename MCBase::IfThenElse(payoff(tree->childs()[0], k), payoff(tree->childs()[1], k), payoff(tree->childs()[2], k)));
            }
            case Scripting::Expression::MIN: {
                if (!hasChilds(tree, 2, k)) QL_FAIL("Cannot interprete payoff");
                if (!hasLeafs(tree, 0, k))  QL_FAIL("Cannot interprete payoff");
                return ext::shared_ptr<PayoffType>(new typename MCBase::Min(payoff(tree->childs()[0], k), payoff(tree->childs()[1], k)));
            }
            case Scripting::Expression::MAX: {
                if (!hasChilds(tree, 2, k)) QL_FAIL("Cannot interprete payoff");
                if (!hasLeafs(tree, 0, k))  QL_FAIL("Cannot interprete payoff");
                return ext::shared_ptr<PayoffType>(new typename MCBase::Max(payoff(tree->childs()[0], k), payoff(tree->childs()[1], k)));
            }
            case Scripting::Expression::EXPONENTIAL: {
                if (!hasChilds(tree, 1, k)) QL_FAIL("Cannot interprete payoff");
                if (!hasLeafs(tree, 0, k))  QL_FAIL("Cannot interprete payoff");
                return ext::shared_ptr<PayoffType>(new typename MCBase::Exponential(payoff(tree->childs()[0], k)));
            }
            case Scripting::Expression::LOGARITHM: {
                if (!hasChilds(tree, 1, k)) QL_FAIL("Cannot interprete payoff");
                if (!hasLeafs(tree, 0, k))  QL_FAIL("Cannot interprete payoff");
                return ext::shared_ptr<PayoffType>(new typename MCBase::Logarithm(payoff(tree->childs()[0], k)));
            }
            case Scripting::Expression::SQUAREROOT: {
                if (!hasChilds(tree, 1, k)) QL_FAIL("Cannot interprete payoff");
                if (!hasLeafs(tree, 0, k))  QL_FAIL("Cannot interprete payoff");
                return ext::shared_ptr<PayoffType>(new typename MCBase::Squareroot(payoff(tree->childs()[0], k)));
            }
            case Scripting::Expression::LOGICAL: {
                if (!hasChilds(tree, 2, k)) QL_FAIL("Cannot interprete payoff");
                if (!hasLeafs(tree, 1, k))  QL_FAIL("Cannot interprete payoff");
                return ext::shared_ptr<PayoffType>(new typename MCBase::Logical(payoff(tree->childs()[0], k), payoff(tree->childs()[1], k), tree->leafs()[0]));
            }
            case Scripting::Expression::PAY: {
                if (!hasChilds(tree, 1, k)) QL_FAIL("Cannot interprete payoff");
                if (!hasLeafs(tree, 1, k))  QL_FAIL("Cannot interprete payoff");
                ActiveType number;
                if (!to_Number(tree->leafs()[0], number)) {
                    scriptLog_.push_back(std::string("Error line " + std::to_string(k) + ": cannot convert " + tree->leafs()[0] + " to number."));
                    QL_FAIL("Cannot interprete payoff");
                }
                return ext::shared_ptr<PayoffType>(new typename MCBase::Pay(payoff(tree->childs()[0], k), number));
            }
            case Scripting::Expression::PAY_WITHDATE: {
                if (!hasChilds(tree, 1, k)) QL_FAIL("Cannot interprete payoff");
                if (!hasLeafs(tree, 1, k))  QL_FAIL("Cannot interprete payoff");
                ActiveType number;
                if (!date_to_Number(tree->leafs()[0], number)) {
                    scriptLog_.push_back(std::string("Error line " + std::to_string(k) + ": cannot convert " + tree->leafs()[0] + " to number."));
                    QL_FAIL("Cannot interprete payoff");
                }
                return ext::shared_ptr<PayoffType>(new typename MCBase::Pay(payoff(tree->childs()[0], k), number));
            }
            case Scripting::Expression::CACHE: {
                if (!hasChilds(tree, 1, k)) QL_FAIL("Cannot interprete payoff");
                if (!hasLeafs(tree, 0, k))  QL_FAIL("Cannot interprete payoff");
                return ext::shared_ptr<PayoffType>(new typename MCBase::Cache(payoff(tree->childs()[0], k)));
            }
            case Scripting::Expression::PAYOFFAT: {
                if (!hasChilds(tree, 1, k)) QL_FAIL("Cannot interprete payoff");
                if (!hasLeafsInRange(tree, 1, 2, k)) QL_FAIL("Cannot interprete payoff");
                ActiveType number;
                if (!to_Number(tree->leafs()[0], number)) {
                    scriptLog_.push_back(std::string("Error line " + std::to_string(k) + ": cannot convert " + tree->leafs()[0] + " to number."));
                    QL_FAIL("Cannot interprete payoff");
                }
                ActiveType sign = 1.0;
                if (tree->leafs().size()>1) { // we expect a negative number
                    if (tree->leafs()[1].compare("-") != 0) {
                        scriptLog_.push_back(std::string("Error line " + std::to_string(k) + ": cannot convert " + tree->leafs()[1] + tree->leafs()[0] + " to negative number."));
                        QL_FAIL("Cannot interprete payoff");
                    }
                    sign = -1.0;
                }
                ext::shared_ptr<PayoffType> p = payoff(tree->childs()[0], k);
                if (p) return p->at(sign*number);
                QL_FAIL("Cannot interprete payoff");
            }
            case Scripting::Expression::PAYOFFAT_WITHDATE: {
                if (!hasChilds(tree, 1, k)) QL_FAIL("Cannot interprete payoff");
                if (!hasLeafs(tree, 1, k))  QL_FAIL("Cannot interprete payoff");
                ActiveType number;
                if (!date_to_Number(tree->leafs()[0], number)) {
                    scriptLog_.push_back(std::string("Error line " + std::to_string(k) + ": cannot convert " + tree->leafs()[0] + " to number."));
                    QL_FAIL("Cannot interprete payoff");
                }
                ext::shared_ptr<PayoffType> p = payoff(tree->childs()[0], k);
                if (p) return p->at(number);
                QL_FAIL("Cannot interprete payoff");
            }
            case Scripting::Expression::NEXT: {
                QL_FAIL("Cannot interprete payoff");
            }
            case Scripting::Expression::ASSIGNMENT: {
                QL_FAIL("Cannot interprete payoff");
            }
            case Scripting::Expression::UNKNOWNTYPE: {
                QL_FAIL("Cannot interprete payoff");
            }
            // we don't need a default because we returned in each of the previous cases
            } // finished all switch types
            // if we end up here there is an expression which we didn't interprete 
            scriptLog_.push_back(std::string("Error line " + std::to_string(k) + ": unknown expression type."));
            QL_FAIL("Cannot interprete payoff");
            return 0; // this should never be reached
        }

        // parse the script and set up payoffs
        inline void parseFlexBisonScript(const std::vector<std::string>&  script,
                                         const bool                       overwrite = true) {
            for (Size k = 0; k < script.size(); ++k) {  // first line should equal 'FlexBison' and is skipped anyway
                Scripting::FlexBisonDriver driver(script[k], false, false);
                // in any case we want to know the parsing result
                if (driver.expressionTree()) expressions_.push_back("L" + std::to_string(k) + ":" + driver.expressionTree()->toString());
                if (driver.returnValue() == 0) {
                    if (!driver.expressionTree()) {
                        scriptLog_.push_back(std::string("Error line " + std::to_string(k) + ": Empty expression tree."));
                        continue;
                    }
                    if (driver.expressionTree()->type() != Scripting::Expression::ASSIGNMENT ) {
                        scriptLog_.push_back(std::string("Error line " + std::to_string(k) + ": Assignment expected."));
                        continue;
                    }
                    if (!hasChilds(driver.expressionTree(), 1, k)) continue;
                    if (!hasLeafs(driver.expressionTree(), 1, k)) continue;
                    // interprete right side of assignment
                    ext::shared_ptr<PayoffType> p;
                    try {
                        p = payoff(driver.expressionTree()->childs()[0], k);
                    }
                    catch (std::exception e) {  // something went wrong, for details check scriptLog_
                        scriptLog_.push_back(std::string("Error line " + std::to_string(k) + ": Exception caught: " + e.what()));
                        continue;
                    }
                    if (!p) {
                        scriptLog_.push_back(std::string("Error line " + std::to_string(k) + ": No payoff found."));
                        continue;
                    }
                    // just define an abbreviation
                    std::string var = driver.expressionTree()->leafs()[0];
                    if (var.compare("") == 0) {
                        scriptLog_.push_back(std::string("Error line " + std::to_string(k) + ": Non-empty identifier expected."));
                        continue;
                    }
                    // now we have a payoff which we may store in the map
                    typename std::map<std::string, ext::shared_ptr<PayoffType> >::iterator it = payoffs_.find(var);
                    if (it == payoffs_.end()) { // insert a new element
                        payoffs_.insert(std::make_pair(var, p));
                        scriptLog_.push_back(std::string("Insert line " + std::to_string(k) + ": '" + script[k] + "'"));
                        continue;
                    }
                    if (overwrite) {
                        it->second = p;
                        scriptLog_.push_back(std::string("Replace line " + std::to_string(k) + ": '" + script[k] + "'"));
                        continue;
                    } 
                    else scriptLog_.push_back(std::string("Error line " + std::to_string(k) + ": Cannot replace line '" + script[k] + "'"));
                }
                else {
                    scriptLog_.push_back(std::string("Error line " + std::to_string(k) + ": " + driver.errorMsg()));
                    continue;
                }
            }
        }

        /*

        we implement the following non-recursive grammar

        line  =  var '=' expr
        var   =  [a-zA-Z][a-zA-Z0-9]*           { RegEx }
        expr  =  operator | function | payoff   { apply from left to right }

        operator   =  operator1 | operator2
        operator1  =  ['+' | '-'] payoff
        operator2  =  payoff ['+' | '-' | '*' | == | != | < | <= | > | >= | && | || ] payoff

        function   =  function3 | function2 | function1
        function3  =  fname3 '(' payoff ',' payoff ',' payoff ')'
        function2  =  fname2 '(' payoff ',' payoff ')'
        function1  =  fname1 '(' payoff ')'

        fname3     =  'IfThenElse'
        fname2     =  'Min' | 'Max | Pay'
        fname1     =  'Cache'

        payoff  =  number | string              { try double conversion and lookup in map }

        */

        // parse the script and set up payoffs
        inline void parseScript(const std::vector<std::string>&  script,
            const bool                       overwrite = true) {
            for (Size k = 0; k < script.size(); ++k) {  // parse lines
                std::string line = boost::regex_replace(script[k], boost::regex(" "), ""); // remove whitespaces
                boost::smatch what;

                bool isAssignment = boost::regex_match(line, what, boost::regex("([a-zA-Z][a-zA-Z0-9]*)(=)(.+)"));
                if (!isAssignment) {
                    scriptLog_.push_back(std::string("Error line " + std::to_string(k) + ": '" + line + "' is no valid assignment"));
                    continue; // move to next item in script
                }
                std::string var(what[1]), expr(what[3]);

                ext::shared_ptr<PayoffType> p;

                if (boost::regex_match(expr, what, boost::regex("(\\+|-)(.+)"))) {
                    p = operator1(std::string(what[1]), std::string(what[2]), k);
                }
                else if (boost::regex_match(expr, what, boost::regex("(.+)(\\+|-|\\*|==|!=|<=|<|>=|>|&&|\\|\\|)(.+)"))) {
                    p = operator2(std::string(what[2]), std::string(what[1]), std::string(what[3]), k);
                }
                else if (boost::regex_match(expr, what, boost::regex("([a-zA-Z]+)\\((.+),(.+),(.+)\\)"))) {
                    p = function3(std::string(what[1]), std::string(what[2]), std::string(what[3]), std::string(what[4]), k);
                }
                else if (boost::regex_match(expr, what, boost::regex("([a-zA-Z]+)\\((.+),(.+)\\)"))) {
                    p = function2(std::string(what[1]), std::string(what[2]), std::string(what[3]), k);
                }
                else if (boost::regex_match(expr, what, boost::regex("([a-zA-Z]+)\\((.+)\\)"))) {
                    p = function1(std::string(what[1]), std::string(what[2]), k);
                }
                else p = payoff(expr, k); // action of last resort

                if (!p) {
                    scriptLog_.push_back(std::string("Error line " + std::to_string(k) + ": '" + expr + "' is no valid expression"));
                    continue; // move to next item in script
                }

                // now we have a payoff which we may store in the map
                typename std::map<std::string, ext::shared_ptr<PayoffType> >::iterator it = payoffs_.find(var);
                if (it == payoffs_.end()) { // insert a new element
                    payoffs_.insert(std::make_pair(var, p));
                    scriptLog_.push_back(std::string("Insert line " + std::to_string(k) + ": '" + line + "'"));
                    continue;
                }
                if (overwrite) {
                    it->second = p;
                    scriptLog_.push_back(std::string("Replace line " + std::to_string(k) + ": '" + line + "'"));
                    continue;
                }
                // if we end up here we have a valid payoff but are not allowed to overwrite existing map entry
                scriptLog_.push_back(std::string("Error line " + std::to_string(k) + ": '" + var + "' can not be replaced"));
            }
            if (script.size() == 0) { // in this case the previous loop was not executed and we just print some help details
                scriptLog_.push_back(std::string("we implement the following non-recursive grammar                                     "));
                scriptLog_.push_back(std::string("                                                                                     "));
                scriptLog_.push_back(std::string("line  =  var '=' expr                                                                "));
                scriptLog_.push_back(std::string("var   =  [a-zA-Z][a-zA-Z0-9]*           { RegEx }                                    "));
                scriptLog_.push_back(std::string("expr  =  operator | function | payoff   { apply from left to right }                 "));
                scriptLog_.push_back(std::string("                                                                                     "));
                scriptLog_.push_back(std::string("operator   =  operator1 | operator2                                                  "));
                scriptLog_.push_back(std::string("operator1  =  ['+' | '-'] payoff                                                     "));
                scriptLog_.push_back(std::string("operator2  =  payoff ['+' | '-' | '*' |                                              "));
                scriptLog_.push_back(std::string("                      '==' | '!=' | '<=' |'<' | '>=' | '>' | '&&' | '||' ] payoff    "));
                scriptLog_.push_back(std::string("                                                                                     "));
                scriptLog_.push_back(std::string("function   =  function3 | function2 | function1                                      "));
                scriptLog_.push_back(std::string("function3  =  fname3 '(' payoff ',' payoff ',' payoff ')'                            "));
                scriptLog_.push_back(std::string("function2  =  fname2 '(' payoff ',' payoff ')'                                       "));
                scriptLog_.push_back(std::string("function1  =  fname1 '(' payoff ')'                                                  "));
                scriptLog_.push_back(std::string("                                                                                     "));
                scriptLog_.push_back(std::string("fname3     =  'IfThenElse'                                                           "));
                scriptLog_.push_back(std::string("fname2     =  'Min' | 'Max' | 'Pay'                                                  "));
                scriptLog_.push_back(std::string("fname1     =  'Cache'                                                                "));
                scriptLog_.push_back(std::string("                                                                                     "));
                scriptLog_.push_back(std::string("payoff  =  number | string              { try double conversion and lookup in map }  "));
            }
        }


        // compile fixed cash flow or lookup in map
        inline ext::shared_ptr<PayoffType> payoff(const std::string expr, const Size lineNr) {
            ActiveType amount;
            bool isFixed = to_Number(expr, amount);
            if (isFixed) {
                scriptLog_.push_back(std::string("Payoff line " + std::to_string(lineNr) + ": '" + boost::lexical_cast<std::string>(amount) + "' is fixed amount"));
                return ext::shared_ptr<PayoffType>(new typename MCBase::FixedAmount(amount));
            }
            typename std::map<std::string, ext::shared_ptr<PayoffType> >::iterator it = payoffs_.find(expr);
            if (it != payoffs_.end()) {
                scriptLog_.push_back(std::string("Payoff line " + std::to_string(lineNr) + ": '" + expr + "' is in map"));
                return it->second;
            }
            // if we end up here no conversion was successfull
            scriptLog_.push_back(std::string("Error line " + std::to_string(lineNr) + ": '" + expr + "' is no payoff"));
            return 0;
        }

        // compile single operand function
        inline ext::shared_ptr<PayoffType> function1(const std::string fname, const std::string operand, const Size lineNr) {
            ext::shared_ptr<PayoffType> p = payoff(operand, lineNr);
            if (!p) {
                scriptLog_.push_back(std::string("Error line " + std::to_string(lineNr) + ": '" + operand + "' is no valid operand"));
                return 0;
            }
            boost::smatch what;
            if (boost::regex_match(fname, what, boost::regex("Cache")))
                return ext::shared_ptr<PayoffType>(new typename MCBase::Cache(p));
            // if we end up here the function name is not valid
            scriptLog_.push_back(std::string("Error line " + std::to_string(lineNr) + ": '" + fname + "' is no valid unary function name"));
            return 0;
        }

        // compile dual operand function
        inline ext::shared_ptr<PayoffType> function2(const std::string fname, const std::string oper1, const std::string oper2, const Size lineNr) {
            ext::shared_ptr<PayoffType> p1 = payoff(oper1, lineNr);
            if (!p1) {
                scriptLog_.push_back(std::string("Error line " + std::to_string(lineNr) + ": '" + oper1 + "' is no valid operand"));
                return 0;
            }
            ext::shared_ptr<PayoffType> p2 = payoff(oper2, lineNr);
            if (!p2) {
                scriptLog_.push_back(std::string("Error line " + std::to_string(lineNr) + ": '" + oper2 + "' is no valid operand"));
                return 0;
            }
            boost::smatch what;
            if (boost::regex_match(fname, what, boost::regex("Min")))
                return ext::shared_ptr<PayoffType>(new typename MCBase::Min(p1, p2));
            if (boost::regex_match(fname, what, boost::regex("Max")))
                return ext::shared_ptr<PayoffType>(new typename MCBase::Max(p1, p2));
            if (boost::regex_match(fname, what, boost::regex("Pay"))) {
                DateType t;
                if (to_Number(oper2, t)) return ext::shared_ptr<PayoffType>(new typename MCBase::Pay(p1, t)); // usual application
                return ext::shared_ptr<PayoffType>(new typename MCBase::Pay(p1, p2->observationTime())); // fall back
            }
            // if we end up here the function name is not valid
            scriptLog_.push_back(std::string("Error line " + std::to_string(lineNr) + ": '" + fname + "' is no valid binary function name"));
            return 0;
        }

        // compile three operand function
        inline ext::shared_ptr<PayoffType> function3(const std::string fname, const std::string oper1, const std::string oper2, const std::string oper3, const Size lineNr) {
            ext::shared_ptr<PayoffType> p1 = payoff(oper1, lineNr);
            if (!p1) {
                scriptLog_.push_back(std::string("Error line " + std::to_string(lineNr) + ": '" + oper1 + "' is no valid operand"));
                return 0;
            }
            ext::shared_ptr<PayoffType> p2 = payoff(oper2, lineNr);
            if (!p2) {
                scriptLog_.push_back(std::string("Error line " + std::to_string(lineNr) + ": '" + oper2 + "' is no valid operand"));
                return 0;
            }
            ext::shared_ptr<PayoffType> p3 = payoff(oper3, lineNr);
            if (!p3) {
                scriptLog_.push_back(std::string("Error line " + std::to_string(lineNr) + ": '" + oper3 + "' is no valid operand"));
                return 0;
            }
            boost::smatch what;
            if (boost::regex_match(fname, what, boost::regex("IfThenElse")))
                return ext::shared_ptr<PayoffType>(new typename MCBase::IfThenElse(p1, p2, p3));
            // if we end up here the function name is not valid
            scriptLog_.push_back(std::string("Error line " + std::to_string(lineNr) + ": '" + fname + "' is no valid function name"));
            return 0;
        }

        // compile unary operators
        inline ext::shared_ptr<PayoffType> operator1(const std::string opname, const std::string operand, const Size lineNr) {
            ext::shared_ptr<PayoffType> p = payoff(operand, lineNr);
            if (!p) {
                scriptLog_.push_back(std::string("Error line " + std::to_string(lineNr) + ": '" + operand + "' is no valid operand"));
                return 0;
            }
            boost::smatch what;
            if (boost::regex_match(opname, what, boost::regex("\\+"))) return p;
            if (boost::regex_match(opname, what, boost::regex("-"))) return ext::shared_ptr<PayoffType>(new typename MCBase::Axpy(-1.0, p, 0));
            // if we end up here the function name is not valid
            scriptLog_.push_back(std::string("Error line " + std::to_string(lineNr) + ": '" + opname + "' is no valid unary operator name"));
            return 0;
        }

        // compile binary operators
        inline ext::shared_ptr<PayoffType> operator2(const std::string opname, const std::string oper1, const std::string oper2, const Size lineNr) {
            ext::shared_ptr<PayoffType> p1 = payoff(oper1, lineNr);
            if (!p1) {
                scriptLog_.push_back(std::string("Error line " + std::to_string(lineNr) + ": '" + oper1 + "' is no valid operand"));
                return 0;
            }
            ext::shared_ptr<PayoffType> p2 = payoff(oper2, lineNr);
            if (!p2) {
                scriptLog_.push_back(std::string("Error line " + std::to_string(lineNr) + ": '" + oper2 + "' is no valid operand"));
                return 0;
            }
            boost::smatch what;
            if (boost::regex_match(opname, what, boost::regex("\\+")))
                return ext::shared_ptr<PayoffType>(new typename MCBase::Axpy(1.0, p1, p2));
            if (boost::regex_match(opname, what, boost::regex("-")))
                return ext::shared_ptr<PayoffType>(new typename MCBase::Axpy(-1.0, p2, p1));
            if (boost::regex_match(opname, what, boost::regex("\\*")))
                return ext::shared_ptr<PayoffType>(new typename MCBase::Mult(p1, p2));
            if (boost::regex_match(opname, what, boost::regex("==|!=|<|<=|>|>=|&&|\\|\\|")))
                return ext::shared_ptr<PayoffType>(new typename MCBase::Logical(p1, p2, opname));
            // if we end up here the function name is not valid
            scriptLog_.push_back(std::string("Error line " + std::to_string(lineNr) + ": '" + opname + "' is no valid binary operator name"));
            return 0;
        }

    };

}

#endif  /* ifndef quantlib_templatemcscript_hpp */ 
