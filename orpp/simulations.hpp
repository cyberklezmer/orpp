#ifndef SIMULATIONS_HPP
#define SIMULATIONS_HPP

namespace orpp
{

struct statcounter
{
    double sum = 0.0;
    double sumsq = 0.0;
    unsigned num = 0;
    void add(double x) { sum += x; sumsq += x*x; num++;}
    double average() const { return sum / num; }
    double var() const { return sumsq/num - average()*average(); }
    double averagevar() const  { return var() / num; }
    double averagestdev() const { return sqrt(averagevar()); }
};


} // namespace

#endif // SIMULATIONS_HPP
