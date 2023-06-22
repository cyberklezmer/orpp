#ifndef BOOSTDIST_H
#define BOOSTDIST_H

#include "orpp/random.hpp"
#include <boost/math/distributions.hpp>

namespace orpp
{


/// \addtogroup realdists Real distributions
/// \ingroup Distributions
/// @{

/// tbd own draw
template <typename B>
class boostdistribution : virtual public qdistribution<nothing>
{
public:
    boostdistribution() {}
    boostdistribution(const B& d) : fd(d) {}
    const B& d() const { return fd; }
private:
    virtual probability cdf_is(double x, const nothing& ) const
    {
        return boost::math::cdf(fd,x);
    }
    virtual double quantile_is(probability p, const nothing& ) const
    {
        return boost::math::quantile(fd,p);
    }
    double pdf(double x) const
    {
        return boost::math::pdf(fd,x);
    }
    B fd;
// tbd
//      virtual double do_draw(const nothing&) const
//    {
//        generator(,fd)
//        return fd();
//    }
};

using stdnormaldistribution
  =boostdistribution<boost::math::normal_distribution<double>>;

using normaldistribution=scaleddistribution<stdnormaldistribution>;

class lognormaldistribution :
        public boostdistribution<boost::math::lognormal_distribution<double>>
{
public:
    lognormaldistribution(double mu, double sigma) :
        boostdistribution(boost::math::lognormal_distribution<double>(mu, sigma))
    {}
    double mu() const { return this->d().location(); }
    double sigma() const { return this->d().scale(); }
};


/// @} // realdists



} // namespace

#endif // BOOSTDIST_H

