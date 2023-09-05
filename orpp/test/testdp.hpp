#ifndef TESTDP_HPP
#define TESTDP_HPP

#include "orpp/dp.hpp"
#include <thread>

namespace orpp
{

struct accuracytestresult 
{
    double requiredaccuracy;
    statcounter details;
    statcounter errors;
    bool maxiterswarning = false;
    bool testresult() const
    {
        return details.stdev() < requiredaccuracy + sqrt(details.varstdev());
    }
};

template <bool raw, typename Problem>
inline accuracytestresult testevaluate(
                     const Problem& problem,
                     index s0index,
                     const std::vector<finitepolicy>& p,
                     double accuracy,
                     const typename Problem::computationparams& params,
                     unsigned testmaxiters
                     )
{
    accuracytestresult res;
    sys::log() << "testevaluate accuracy="
               << accuracy << std::endl;
    res.maxiterswarning = false;
    for(unsigned i=0; i<testmaxiters; i++)
    {
        sys::log() << "Iteration " << i << " ";
        double r;
        double e;
        unsigned num;
        if constexpr(raw)
        {
            auto resr = problem.evaluateraw(s0index,p,accuracy,params);
            r = resr.average();
            e = resr.averagestdev();
            num = resr.num;
            if(num==params.computationparams)
                res.maxiterswarning = true;
        }
        else
        {
           auto res = problem.evaluatecrit(s0index,p,accuracy,params);
           r = res.x;
           e = res.sd;
           num = 1;
        }

        res.details.add(r);
        res.errors.add(e);
        double acc = sqrt(res.details.varstdev());
        sys::log() << " value=" << r << "(" << e << ")"
                      << " n=" << num <<
                      " err=" << res.details.stdev()
                   << "(" << acc << ")" << std::endl;
        if(i > 4 && res.details.stdev() + acc < accuracy)
            break;
    }
    res.requiredaccuracy=accuracy;
    sys::log() << "Result: required=" << res.requiredaccuracy
               << " reported=" << res.errors.average()
               << " acquired=" << res.details.stdev()
               << "(" << sqrt(res.details.varstdev())
               << ") -> "
               <<  (res.testresult() ? "Passed" : "Failed") << std::endl;
    return res;
}

template <typename Problem>
inline accuracytestresult testevaluatehomo(
                     const Problem& problem,
                     index s0index,
                     finitepolicy& p,
                     double accuracy,
                     const typename Problem::computationparams& params,
                     unsigned testmaxiters
                     )
{
    accuracytestresult res;
    sys::log() << "testevaluatehomo accuracy="
               << accuracy << std::endl;

    finitevaluefunction initialV(problem);

    for(unsigned i=0; i<testmaxiters; i++)
    {
        for(unsigned i=0; i<initialV.size(); i++)
            initialV[i] = problem.maxreward() / (1-problem.gamma()) * sys::uniform();
        sys::log() << "Iteration " << i << " ";
        auto rr = problem.evaluate(initialV,p,accuracy,params);
        double r = rr.x[s0index];
        res.details.add(r);
        res.errors.add(rr.sd);
        double acc = sqrt(res.details.varstdev());
        sys::log() << " value=" << r << " err=" << res.details.stdev()
                   << "(" << acc << ")" << std::endl;
        if(i > 4 && res.details.stdev() + acc < accuracy)
            break;
    }
    res.requiredaccuracy=accuracy;
    sys::log() << "Result: required=" << res.requiredaccuracy
               << " reported=" << res.errors.average()
               << " acquired=" << res.details.stdev()
               << "(" << sqrt(res.details.varstdev())
               << ") -> "
               <<  (res.testresult() ? "Passed" : "Failed") << std::endl;
    return res;
}


} // namespace

#endif // TESTDP_HPP
