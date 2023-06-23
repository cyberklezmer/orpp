#ifndef SIMULATIONS_HPP
#define SIMULATIONS_HPP

#include <cmath>
#include "orpp/random.hpp"

namespace orpp
{

class empiricaldistribution : public equipdistribution<double, true>
{
public:
    empiricaldistribution()
        : equipdistribution<double, true> (std::vector<double>(),true)
    {
    }

    empiricaldistribution(const std::vector<double>& values)
        : equipdistribution<double, true> (values,true)
    {
    }
    void add(double v)
    {
        auto insertionPos = std::lower_bound(values().begin(), values().end(), v);
        values().insert(insertionPos, v);
    }
    double value(unsigned i) const { return (*this)(i).x; }
    double p() const { return 1.0 / (*this).natoms(); }
};

struct statcounter
{
    empiricaldistribution dist;
    double sum = 0.0;
    double sumsq = 0.0;
    unsigned num = 0;
    void add(double x) { sum += x; sumsq += x*x; num++; dist.add(x); }
    double var() const { return sumsq/num - average()*average(); }
    double averagevar() const  { return var() / num; }
    double averagestdev() const { return sqrt(averagevar()); }
    double average() const { return sum / num; }
    realestimate estaverage() const { return {average(),averagestdev()};  }    
};


class empiricalMeanCVaR
{
public:
    empiricalMeanCVaR(probability alpha, double lambda) : falpha(alpha),
        flambda(lambda) {}
    virtual realestimate operator () (const empiricaldistribution& d) const
    {
        statcounter sc;
        unsigned n = d.natoms();
        assert(n);


        double p = 1-falpha;
        double s = 0;
        for(int i = d.natoms()-1;i>=0 ;i--)
        {
            double x = d.value(i);
            double delta = std::min(1.0 / n, p-s);
            double inc = x * (1-flambda);
            if(delta >= 0)
            {
               inc += d.value(i) / p * delta * n * flambda;
               s += delta;
            }
            sc.add(inc);
        }
        return sc.estaverage();
    }

private:
    probability falpha;
    double flambda;

};


} // namespace

#endif // SIMULATIONS_HPP
