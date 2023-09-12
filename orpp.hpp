#ifndef ORPP_HPP
#define ORPP_HPP

#include <stdexcept>
#include <iostream>
#include <vector>
#include <chrono>
#include <sstream>
#include <random>
#include <memory>
#include <assert.h>
#include <math.h>
#include <time.h>
#include <memory>


namespace orpp
{


/** \addtogroup sysclasses System Classes
 *  @{
 */

// TBD inherit all fron this
/// Base class
class object
{
public:
    virtual ~object() {}
};

/// \ingroup General Random
using probability = double;


class exception: public std::runtime_error
{
public:
    exception(const std::string& s) : std::runtime_error(s) {}
    exception(const char* s) : std::runtime_error(s) {}
    exception(const std::ostringstream& s) : std::runtime_error(s.str()) {}
};

class sys
{
    static sys& self()
    {
       static sys s;
       return s;
    }
public:
    sys() : fout(0), flog(0), floglevel(0), ferr(0), funiform(0.0,1.0),
        fstarttime(abstimems()), flastlinetime(0)
    {
    }
    ~sys()
    {
    }

    /// Sets an existing stream as the text output of the library.
    static void setout(std::ostream& o)
    {
        self().fout = &o;
    }

    /// Resets the standard ouptut of the library,
    /// redirecting it back to \p std::cout
    static void resetout()
    {
        self().fout = 0;
    }

    static void setlog(std::ostream& l)
    {
        self().flog = &l;
    }
    static void resetlog()
    {
        self().flog=0;
    }


    static void setloglevel(unsigned l)
    {
        self().floglevel = l;
    }

    static void seterr(std::ostream& e)
    {
        self().ferr = &e;
    }
    static void reseterr()
    {
        self().ferr = 0;
    }

    static void setoutputfolder(const std::string& afolder)
    {
        self().foutputfolder = afolder;
    }

    static void settmpfolder(const std::string& afolder)
    {
        self().ftmpfolder = afolder;
    }

    static std::ostream& out()
    {
        return self().fout ? *(self().fout) : std::cout;
    }
    static std::ostream& log()
    {
        return self().flog ? *(self().flog) : std::clog;
    }
    static std::ostream& logline(unsigned level=0)
    {
        auto tm = timems();
        log() << timems() << "," << tm - self().flastlinetime << ",";
        self().flastlinetime = tm;
        for(unsigned i=0; i<level; i++)
        {
           log() << ",";
        }
        return log();
    }
    static unsigned loglevel()
    {
        return self().floglevel;
    }

    static std::ostream& err()
    {
        return self().ferr ? *(self().ferr) : std::cerr;
    }

    static const std::string& outputfolder()
    {
        return self().foutputfolder;
    }

    static const std::string& tmpfolder()
    {
        return self().ftmpfolder;
    }

    static unsigned timems()
    {
        return abstimems() - self().fstarttime;
    }

    static unsigned abstimems()
    {
        auto now_ms = std::chrono::time_point_cast<std::chrono::milliseconds>(std::chrono::system_clock::now());
        auto epoch = now_ms.time_since_epoch();
        auto value = std::chrono::duration_cast<std::chrono::milliseconds>(epoch);
        return value.count();
    }


    /// Resets the seed of the random generator according to computer time.
    static void seed()
    {
        self().fengine.seed(time(0));
    }

    /// Resets the seed of the random generator.
    static void seed(unsigned int aseed)
    {
        self().fengine.seed(aseed);
    }

    /// Returns a pesudo-random observaton from the uniform distribution on [0,1]
    static double uniform()
    {
        return self().funiform(self().fengine);
    }

    static std::default_random_engine& engine()
    {
        return self().fengine;
    };

    static void set(const std::string& var, const std::string& val)
    {
        for(unsigned i=0; i<self().fenv.size(); i++)
        {
            if(self().fenv[i].first == var)
            {
                self().fenv[i].second = val;
                return;
            }
        }
        self().fenv.push_back({var,val});
    }
    static std::pair<bool,std::string> get(const std::string& var)
    {
        for(unsigned i=0; i<self().fenv.size(); i++)
        {
            if(self().fenv[i].first == var)
                return { true, self().fenv[i].second };
        }
        return { false, "" };
    }
private:
    std::ostream* fout;
    std::ostream* flog;
    unsigned floglevel;
    std::ostream* ferr;
    std::default_random_engine fengine;
    std::uniform_real_distribution<double> funiform;
    std::string foutputfolder;
    std::string ftmpfolder;
    std::vector<std::pair<std::string,std::string>> fenv;
    unsigned fstarttime;
    unsigned flastlinetime;
};

/// Heap pointer to \p T
template <typename T>
using ptr=std::shared_ptr<T>;

/// @}



/// \addtogroup general General Definitions
///  @{

template <typename T=double>
constexpr T nan=std::numeric_limits<T>::quiet_NaN();

template <typename T=double>
constexpr T infinity=std::numeric_limits<T>::infinity();

template <typename T=double>
constexpr T minfinity = -infinity<T>;

template <typename T=double>
constexpr T max = std::numeric_limits<T>::max();

template <typename T>
constexpr T min = std::numeric_limits<T>::min();


struct nothing {};

using index = unsigned int;


/// \addtogroup fns Functions
/// @{

/// A class bearing no information (alternative to \p void)

constexpr nothing na;

/// \brief Base class for mappings
/// \tparam D domain
/// \tparam R range
template <typename D,typename R>
class mapping: public object
{
public:
    using D_t=D;
    using R_t=R;
    virtual R operator() (const D& ) const = 0;
};

/// Function (mapping from \p D to real numbers)
template <typename D>
using function = mapping<D,double>;

/// Identity mapping
template <typename I>
class idmapping: public mapping<I,I>
{
public:
    virtual I operator() (const I& i) const { return i; }
};

/// \brief Mapping having no value
///
/// Used as a condition of unconditional distributions, for instance.
template <typename I>
class nomapping: public mapping<I,nothing>
{
public:
    virtual nothing operator() (const I& i) const { return na; }
};

/// Function from Eucleidean space
class efunction: public function<std::vector<double>>
{
public:
    efunction(unsigned int xdim) : fxdim(xdim) {}
    unsigned int xdim() const { return fxdim; }
    double operator() (const std::vector<double>& x) const
    {
        assert(x.size() == fxdim);
        return value_is(x);
    }
private:
    unsigned int fxdim;
    virtual double value_is(const std::vector<double>& x) const = 0;
};

/// Convex function from an Eucleidean space
class convexefunction: public efunction
{
public:
    convexefunction(unsigned int xdim) : efunction(xdim) {}
};

/// Function with defined subdifferential
class subdifefunction: public convexefunction
{
public:
    subdifefunction(unsigned int dim) : convexefunction(dim) {}

    std::vector<double> sg(const std::vector<double>& x) const
    {
        return sg_is(x);
    }
private:
   virtual std::vector<double> sg_is(const std::vector<double>& x) const = 0;
};

/// Linear function
class linearfunction: public subdifefunction
{
public:
    linearfunction(unsigned int dim):
        subdifefunction(dim), fc(dim) {}

    linearfunction(std::vector<double> c):
        subdifefunction(c.size()), fc(c) {}

    void setc(const std::vector<double>& c)
    {
        assert(c.size()==xdim());
        fc = c;
    }
    void setc(unsigned int i, double c)
    {
        assert(i<fc.size());
        fc[i] = c;
    }
    double c(unsigned int i) const
    {
        assert(i < fc.size());
        return fc[i];
    }
    const std::vector<double>& c() const
    {
        return fc;
    }

private:
    virtual std::vector<double> sg_is(const std::vector<double>& x) const
    {
        assert(x.size()==xdim());
        return fc;
    }
    virtual double value_is(const std::vector<double>& x) const
    {
        assert(x.size()==xdim());
        double s=0;
        for(unsigned int i=0; i<fc.size(); i++)
            s += fc[i] * x[i];
        return s;
    }
private:
    std::vector<double> fc;
};


///@}



} // namespace

#endif // ORPP_HPP

