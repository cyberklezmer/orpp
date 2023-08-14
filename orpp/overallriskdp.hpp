#ifndef OVERALLRISKDP_HPP
#define OVERALLRISKDP_HPP

#include "orpp/dp.hpp"
#include "orpp/random.hpp"

namespace orpp
{


using finitedpstatespace = integerspace;
using finitedpactionspace =  constrainedspace<unsigned int,unsigned int>;
using finitedptransition = finitetransition<finitedpactionspace, finitedpactionspace>;
using finitedpreward = dpreward<finitedpactionspace, finitedpactionspace>;

template <typename Criterion>
class finitedpproblem : public dpproblem
        <Criterion,
         integerspace,
         finitedpactionspace,
         finitedptransition,
         finitedpreward
        >
{
public:
    using Statespace_t = Statespace;
    using ConstrainedActionSpace_t = ConstrainedActionSpace;
public:
    finiteoverallriskdp(const Criterion& crit,
              const Statespace& state,
              const ConstrainedActionSpace& constraint,
              const Transition& transition,
               const Reward& reward,
              double gamma ) :
         fcrit(crit), fstatespace(state),
         fconstraint(constraint),
         ftransition(transition), freward(reward), fgamma(gamma)
    {}
    double gamma() const { return fgamma; }
    const ConstrainedActionSpace& constraint() const { return fconstraint; }
protected:
    Criterion fcrit;
    Statespace fstatespace;
    ConstrainedActionSpace fconstraint;
    Transition ftransition;
    Reward freward;
    double fgamma;
};

class overallriskproblem : public finitedpproblem<testcrit, teststatespace,
        testactionspace, testtransition, testreward>
{
public:
    overallriskproblem(probability alpha, probability pincrease, double gamma) :
        finitedpproblem<testcrit, teststatespace,
                    testactionspace, testtransition, testreward>
          (testcrit(alpha),teststatespace(), testactionspace(),
           testtransition(pincrease), testreward(), gamma) {}
    void setriskaversion(probability alpha)
    {
        fcrit = testcrit(alpha);
    }
};



} // namespace

#endif // OVERALLRISKDP_HPP
