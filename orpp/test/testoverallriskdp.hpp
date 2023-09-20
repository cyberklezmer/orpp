#ifndef TESTOVERALLRISKDP_HPP
#define TESTOVERALLRISKDP_HPP

#include "orpp/overallriskdp.hpp"

namespace orpp
{

template <typename Problem>
inline bool testlipshitzproperty(
                     const Problem& problem,
                     index s0index,
                     const finitepolicy p,
                     double accuracy,
                     const typename Problem::computationparams& params,
                     unsigned numpoints
                     )
{
    sys::logline() << "testlipshitzproperty" << std::endl;
    double ra = problem.crit().getparam();
    if(fabs(ra) <= std::numeric_limits<double>::epsilon())
        throw exception("Zero param of crit!");
    auto res0 = problem.evaluatecrit(s0index,p,accuracy, params);
    double crit = res0.x;
    double l=problem.lipschitzconstant();
    sys::logline() << "Reference ra=" << ra << " crit=" << crit
                   << " assumed lipschitz=" << l <<std::endl;
    bool lfailed = false;
    bool mfailed = false;
    for(unsigned i=0; i<numpoints; i++)
    {
        Problem problemcopy(problem);
        double newra = ra * static_cast<double>(i)/numpoints;
        problemcopy.setriskaversion(newra);
        auto res = problemcopy.evaluatecrit(s0index,p,accuracy, params);
        double newcrit = res.x;
        double realizedlipschitz = (newcrit-crit) / (ra-newra);
        if(realizedlipschitz < - 2*accuracy)
            mfailed = true;
        if(realizedlipschitz - 2*accuracy / (ra-newra) > l)
            lfailed = true;
        sys::logline() << "newra=" << newra << " crit=" << newcrit
                       << " realizedlipschitz=" << realizedlipschitz;
        if(mfailed)
            sys::log() << " positivity test failed";
        if(lfailed)
            sys::log() << " lipschitz test failed";
        sys::log() << std::endl;
    }
    if(lfailed || mfailed)
        return false;
    else
        return true;
}

} // namespace

#endif // TESTOVERALLRISKDP_HPP
