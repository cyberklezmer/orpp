#ifndef BOOSTDIST_H
#define BOOSTDIST_H

#include "orpp/random.hpp"
#include <boost/math/distributions.hpp>
#include <boost/math/distributions/binomial.hpp>

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

template <typename B>
class boostfddistribution : virtual public fdistribution<double,nothing>
{
public:
    boostfddistribution(const B& d, unsigned n) : fd(d), fn(n) {}
    const B& d() const { return fd; }
private:
    virtual bool is_sorted() const
    {
        return true;
    }

    virtual unsigned int natoms_is(const nothing&) const
    {
        return fn;
    }

    virtual atom<double> atom_is(unsigned int i, const nothing&) const
    {
        return {static_cast<double>(i),pdf(fd,static_cast<double>(i))};
    }

    B fd;
    unsigned fn;
};

class binomialdistribution:
   public boostfddistribution<boost::math::binomial_distribution<double>>
{
public:
    binomialdistribution(unsigned n, probability p) :
        boostfddistribution<boost::math::binomial_distribution<double>>
        (boost::math::binomial_distribution<double>(n,p),n+1)
    {}
};

using stdnormaldistribution
  =boostdistribution<boost::math::normal_distribution<double>>;

using normaldistribution=scaleddistribution<stdnormaldistribution>;

class   lognormaldistribution :
        public boostdistribution<boost::math::lognormal_distribution<double>>
{
public:
    lognormaldistribution(double mu, double sigma) :
        boostdistribution(boost::math::lognormal_distribution<double>(mu, sigma))
    {}
    double mu() const { return this->d().location(); }
    double sigma() const { return this->d().scale(); }
};

//class binomialdistribution :
//        public boostdistribution<boost::math::lognormal_distribution<double>>


/// @} // realdists



} // namespace

#endif // BOOSTDIST_H

