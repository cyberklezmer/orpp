#ifndef ORPP_HPP
#define ORPP_HPP

#include <stdexcept>
#include <fstream>
#include <iostream>
#include <vector>
#include <sstream>
#include <random>
#include <assert.h>
#include <math.h>
#include <time.h>
#include <nlopt.hpp>


namespace orpp
{


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
    sys() : fout(0), flog(0), ferr(0), funiform(0.0,1.0)
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
    std::ostream* ferr;
    std::default_random_engine fengine;
    std::uniform_real_distribution<double> funiform;
    std::string foutputfolder;
    std::string ftmpfolder;
    std::vector<std::pair<std::string,std::string>> fenv;
};

constexpr double na=std::numeric_limits<double>::quiet_NaN();
constexpr double infinity=std::numeric_limits<double>::infinity();


/*

inline double pi() { return std::atan(1)*4; }

template <char sep>
inline void csvout(std::ostream& os, const std::string& s)
{
    os << '"';
    for(unsigned i=0; i<s.length(); i++)
        if(s[i] == '"')
            os << '"' << '"';
        else
            os << s[i];
    os << '"' << sep;
}

template <char sep = ','>
class csv : private std::vector<std::vector<std::string>>
{
    static std::string getstr(std::istream& is)
    {
        char c;
        is >> c;
        bool quotes = false;
        if(c=='"')
        {
            quotes = true;
            is >> c;
        }

        std::string r;
        for(;;)
        {
            if((signed)c > 128)
                throw "Czech characters cannot be used";
            if((c==sep && !quotes) || is.eof())
                break;
            if(c=='"' && quotes)
            {
               is >> c;
               if(c!='"')
                  break;
            }
            r.push_back(c);
            is >> c;
        }
        return r;
    }

    static std::vector<std::vector<std::string>> getvecs(const std::string& fn)
    {
        std::vector<vector<std::string>> res;
        std::ifstream file(fn);
        if(!file.is_open())
        {
            std::ostringstream es;
            es << "cannot open csv file: " << fn ;
            throw es.str();
        }
        std::string line="";
        while (getline(file, line))
        {
            std::vector<std::string> row;
            std::istringstream is(line);
            for(;;)
            {
                std::string s = getstr(is);
                row.push_back(s);
                if(is.eof())
                    break;
            }
            res.push_back(row);
        }

        return res;
    }

public:
    unsigned r() const { return (*this).size();}
    unsigned c(unsigned r) const
    {
        assert(r < this->size());
        return (*this)[r].size();
    }
    const std::string& operator()(unsigned i, unsigned j) const
    {
        std::ostringstream e;
        if(i>(*this).size())
        {
            e << "csv does not have" << i+1 << " lines";
            throw e.str();
        }
        if(j>(*this)[i].size())
        {
            e << "Line " << i+1 << " does not have" << j+1 << " columns";
            throw e.str();
        }
        return (*this)[i][j];
    }

    unsigned getunsigned(unsigned i, unsigned j) const
    {
        unsigned res;
        std::string c = (*this)(i,j);
        try {
              res = stoul(c);
          } catch (...) {
              std::ostringstream e;
              e << "Cannot convert '" << c << "' to unsigned (row="
                << i << ", col=" << j << ")";
              throw e.str();
          }
        return res;
    }

    double getdouble(unsigned i, unsigned j) const
    {
        double res;
        std::string c = (*this)(i,j);
        try {
              res = stod(c);
          } catch (...) {
              std::ostringstream e;
              e << "Cannot convert '" << c << "' to double (row="
                << i << ", col=" << j << ")";
              throw e.str();
          }
        return res;
    }




    csv(const std::string& fn) : std::vector<std::vector<std::string>>(getvecs(fn)) {}
};


template<char sep>
inline void csvout(std::ostream& os, const std::vector<std::string>& s)
{
    for(unsigned i=0; i<s.size(); i++)
        csvout<sep>(os, s[i]);
    os << std::endl;
}


template<typename T, char sep>
inline void csvout(std::ostream& os, const std::vector<T>& s)
{
    for(unsigned i=0; i<s.size(); i++)
        os << s[i] << sep;
    os << std::endl;
}

*/

/*using dvector=std::vector<double>;

inline dvector stackv(const dvector& x, const dvector& y)
{
    dvector ret = x;
    for(unsigned i=0; i<y.size(); i++)
        ret.push_back(y[i]);
    return ret;
}
*/

/*
std::ostream& operator << (std::ostream& o, const dvector& v)
{
    for(unsigned i=0; ; i++)
    {
        o << v[i];
        if(++i=v.size())
             break;
        o << ',';
    }
    o << std::endl;
    return o;
}

inline double scalarproduct(const dvector& a,
                            const dvector& b)
{
    assert(a.size()==b.size());
    double s = 0;
    for(unsigned i=0; i<a.size(); i++)
        s+=a[i]*b[i];
    return s;
}

*/


} // namespace

#endif // ORPP_HPP
