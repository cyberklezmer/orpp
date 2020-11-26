#ifndef FUTUREORPP_HPP
#define FUTUREORPP_HPP

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
#define BOOST_UBLAS_TYPE_CHECK 0
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>



namespace orpp
{

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

inline unsigned countfalses(const std::vector<bool>& m)
{
    unsigned n=0;
    for(unsigned i=0; i<m.size(); i++)
        if(!m[i])
            n++;
    return n;
}

template <char sep>
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

inline double scalarproduct(const std::vector<double>& a,
                            const std::vector<double>& b)
{
    assert(a.size()==b.size());
    double s = 0;
    for(unsigned i=0; i<a.size(); i++)
        s+=a[i]*b[i];
    return s;
}

class matrix;
inline std::ostream& operator<<(std::ostream& str,const matrix& r);

class matrix
{
public:
    matrix(const std::vector<double>& v) : fx(v.size(),std::vector<double>(1))
    {
        for(unsigned i=0; i<v.size(); i++)
            fx[i][0] = v[i];
    }
    matrix(unsigned r=0, unsigned c=0, double iv=0.0) :
      fx(r,std::vector<double>(c,iv)) {}

    matrix(const std::vector<std::vector<double>>& x): fx(x)
    {
        assert(x.size()>0);
        unsigned m=x[0].size();
        assert(m>0);
        for(unsigned i=1; i<x.size(); i++)
            assert(x[i].size() == m);
    }
    unsigned c() const { return fx[0].size(); }
    unsigned r() const { return fx.size(); }
    double& operator () ( unsigned i, unsigned j )
    { assert(i<r()); assert(j<c()); return fx[i][j]; }
    double operator () ( unsigned i, unsigned j ) const
    { assert(i<r()); assert(j<c()); return fx[i][j]; }

    template <char sep>
    void csvout(std::ostream& os, bool trans=false, unsigned emptycols=0) const
    {
        if(!trans)
            for(unsigned int i=0; i<r(); i++)
            {
                for(unsigned  j=0; j<emptycols; j++)
                    os << sep;
                for(unsigned  j=0; j<c(); j++)
                    os << fx[i][j] << sep;
                os << std::endl;
            }
        else
            for(unsigned int j=0; j<c(); j++)
            {
                for(unsigned  j=0; j<emptycols; j++)
                    os << sep;
                for(unsigned int i=0; i<r(); i++)
                    os << fx[i][j] << sep;
                os << std::endl;
            }
    }
    const std::vector<double>& operator[] (unsigned i) const
    {
        assert(i<r());
        return fx[i];
    }
/*    operator std::vector<double>() const
    {
        assert(c()==1);
        std::vector<double> ret(r());
        for(unsigned i=0; i<ret.size(); i++)
            ret[i] = (*this)(i,0);
        return ret;
    }*/

    matrix submatrix(unsigned r, unsigned c, unsigned nr, unsigned nc) const
    {
        assert(c+nc <= this->c());
        assert(r+nr <= this->r());
        matrix ret(nr,nc);
        for(unsigned i=0; i<nr; i++)
            for(unsigned j=0; j<nc; j++)
                ret(i,j)=(*this)(r+i,c+j);
        return ret;
    }

    matrix transpose() const
    {
        matrix ret(c(),r());
        for(unsigned i=0; i<r(); i++)
            for(unsigned j=0; j<c(); j++)
                ret(j,i)=(*this)(i,j);
        return ret;
    }

    matrix inverse()
    {
        assert(c()==r());
        unsigned n = c();
        matrix ret(n,n);
        boost::numeric::ublas::matrix<double> A(n,n);
        for(unsigned i=0; i<n; i++)
            for(unsigned j=0; j<n; j++)
                A(i,j)=(*this)(i,j);
        using namespace boost::numeric::ublas;
        typedef permutation_matrix<size_t> pmatrix;
        // create a permutation matrix for the LU-factorization
        pmatrix pm(A.size1());

        // perform LU-factorization
        int res = lu_factorize(A,pm);
           if( res != 0 )
           {
               std::cerr << "cannot invert matrix" << std::endl;
               std::cerr << *this << std::endl;
               throw "cannot invert matrix";
           }

        boost::numeric::ublas::matrix<double> inverse(n,n);
        // create identity matrix of "inverse"
        inverse.assign(boost::numeric::ublas::identity_matrix<double>(A.size1()));

        // backsubstitute to get the inverse
        lu_substitute(A, pm, inverse);

        for(unsigned i=0; i<n; i++)
            for(unsigned j=0; j<n; j++)
                ret(i,j)=inverse(i,j);
        return ret;
    }


    bool iszero()
    {
        for(unsigned i=0; i<r(); i++)
            for(unsigned j=0; j<c(); j++)
                if((*this)(i,j) != 0)
                    return false;
        return true;
    }
protected:
    std::vector<std::vector<double>> fx;
};

inline std::vector<double> vec(const matrix& m)
{
    assert(m.c()==1);
    std::vector<double> ret(m.r());
    for(unsigned i=0; i<ret.size(); i++)
        ret[i] = m(i,0);
    return ret;
}

inline std::vector<double> stackv(const std::vector<double>& x, const std::vector<double>& y)
{
    std::vector<double> ret = x;
    for(unsigned i=0; i<y.size(); i++)
        ret.push_back(y[i]);
    return ret;
}

inline matrix stackv(const matrix& x, const matrix& y)
{
    assert(x.c()==y.c());

    matrix ret(x.r()+y.r(),x.c());
    unsigned i=0;
    for(; i<x.r(); i++)
        for(unsigned j=0; j<x.c(); j++)
            ret(i,j) = x(i,j);
    for(unsigned k=0; k<y.r(); k++,i++)
        for(unsigned j=0; j<y.c(); j++)
            ret(i,j) = y(k,j);
    return ret;
}

inline matrix block(const matrix& x, const matrix& y,
                    const matrix& z, const matrix& zz)
{
    assert(x.r()==y.r());
    assert(z.r()==zz.r());
    assert(x.c()==z.c());
    assert(y.c()==zz.c());

    matrix ret(x.r()+z.r(),x.c()+y.c());
    unsigned i=0;
    for(; i<x.r(); i++)
    {
        unsigned j=0;
        for(j=0; j<x.c(); j++)
            ret(i,j) = x(i,j);
        for(unsigned k=0; k<y.c(); k++, j++)
            ret(i,j) = y(i,k);
    }
    for(unsigned k=0; k<z.r(); k++,i++)
    {
        unsigned j=0;
        for(; j<z.c(); j++)
            ret(i,j) = z(k,j);
        for(unsigned kk=0; kk<zz.c(); kk++,j++)
            ret(i,j) = zz(k,kk);
    }
    return ret;
}


inline matrix diag(const std::vector<double>& v)
{
    matrix ret(v.size(), v.size(), 0);
    for(unsigned i=0; i<v.size(); i++)
        ret(i,i) = v[i];
    return ret;
}


inline matrix operator * (double a, const matrix& m)
{
    matrix r(m.r(),m.c());
    for(unsigned i=0; i<m.r();i++)
        for(unsigned j=0; j<m.c();j++)
            r(i,j) = a * m(i,j);
    return r;
}



inline matrix operator + (const matrix& a, const matrix& b)
{
    assert(a.r()==b.r());
    assert(a.c()==b.c());
    matrix r(a.r(),a.c());
    for(unsigned i=0; i<a.r();i++)
        for(unsigned j=0; j<a.c();j++)
            r(i,j) = a(i,j) + b(i,j);
    return r;
}

inline matrix operator - (const matrix& a, const matrix& b)
{
    assert(a.r()==b.r());
    assert(a.c()==b.c());
    matrix r(a.r(),a.c());
    for(unsigned i=0; i<a.r();i++)
        for(unsigned j=0; j<a.c();j++)
            r(i,j) = a(i,j) - b(i,j);
    return r;
}


inline matrix operator * (const matrix& a, const matrix& b)
{
    assert(a.c() == b.r());
    matrix ret(a.r(),b.c(),0);
    for(unsigned i=0; i<ret.r(); i++)
        for(unsigned j=0; j<ret.c(); j++)
        {
            for(unsigned k=0; k<a.c(); k++)
                ret(i,j) += a(i,k)*b(k,j);
        }
    return ret;
}

inline std::ostream& operator<<(std::ostream& str,const matrix& r)
{
    for(unsigned i=0; i<r.r(); i++)
    {
        for(unsigned j=0; ; j++)
        {
            str << r(i,j);
            if(j==r.c()-1)
                break;
            str << ",";
        }
        str << std::endl;
    }
    return str;
}



inline double cdist(const matrix& a, const matrix& b)
{
    assert(a.r()==b.r());
    assert(a.c()==b.c());
    double sum = 0;
    for(unsigned i=0;i<a.r();i++)
        for(unsigned j=0;j<a.c();j++)
            sum += fabs(a(i,j)-b(i,j));
    return sum;
}

inline matrix transpose(std::vector<double> v)
{
    matrix res(1,v.size());
    for(unsigned i=0; i<v.size(); i++)
        res(0,i) = v[i];
    return res;
}

inline double minoffdiag(const matrix& s)
{
    assert(s.r()==s.c());
    double m=HUGE_VAL;
    for(unsigned i=0; i<s.r(); i++)
        for(unsigned j=0; j<s.r(); j++)
            if(i!=j && s(i,j) < m)
                m = s(i,j);
    return m;
}

inline double maxoffdiag(const matrix& s)
{
    assert(s.r()==s.c());
    double m=-HUGE_VAL;
    for(unsigned i=0; i<s.r(); i++)
        for(unsigned j=0; j<s.r(); j++)
            if(i!=j && s(i,j) > m)
                m = s(i,j);
    return m;
}

inline matrix normoffdiag(const matrix& aw)
{
    double sw = 0.0;
    for(unsigned int i=0; i<aw.r(); i++)
        for(unsigned int j=0; j<aw.c(); j++)
            if(i!=j)
                sw += aw(i,j);

    double avew = sw / (aw.r()*(aw.c()-1));

    return 1.0 / avew * aw;
}

inline double aveoffdiag(const matrix& s,
           const std::vector<bool>& muted = std::vector<bool>())
{
    assert(s.r()==s.c());
    assert(muted.size()==0 || muted.size()==s.r());
    unsigned n = muted.size() ?  countfalses(muted) : s.r();
    assert(n > 0);
    if(n=1)
        return 0;
    double sum = 0.0;
    for(unsigned i=0; i<n; i++)
        for(unsigned j=0; j<n; j++)
            if(i!=j && (muted.size()==0 || (!muted[i] && !muted[j])) )
                sum += s(i,j);
    return sum / ( n * (n-1) );
}

inline matrix relchangeoffdiag(const matrix& nm, const matrix& om)
{
    assert(nm.r()==nm.c());
    assert(om.r()==om.c());
    assert(om.r()==nm.c());
    unsigned n=nm.r();
    matrix res(n,n,0);
    for(unsigned i=0; i<n; i++)
        for(unsigned j=0; j<n; j++)
            if(i != j)
                res(i,j) = nm(i,j) / om(i,j);
    return res;
}

template <typename D>
inline D sum(const std::vector<std::vector<D>> c)
{
    D sum = 0;
    for(unsigned i=0; i<c.size(); i++)
        for(unsigned j=0; j<c[i].size(); j++)
            sum += c[i][j];
    return sum;
}

template <typename D>
inline D sum(const std::vector<D> c)
{
    D sum = 0;
    for(unsigned i=0; i<c.size(); i++)
            sum += c[i];
    return sum;
}


template<typename T>
inline void mcsvout(std::ostream& os, const std::vector<std::vector<T>>& s)
{
    for(unsigned i=0; i<s.size(); i++)
    {
        for(unsigned j=0; j<s[i].size(); j++)
            os << s[i][j] << ",";
        os << std::endl;
    }
    os << std::endl;
}


class gnuplot
{
    std::string name;
    std::ofstream sc;
    std::ofstream dat;
    static constexpr char* klabel = "gnuplot";
    static std::string scriptfn(const std::string& an)
       { return sys::tmpfolder() + an + ".plt"; }
    static std::string datfn(const std::string& an)
       { return sys::tmpfolder() + an + ".dat"; }
public:
    const std::string datfn() { return datfn(name); }
    std::ofstream& script() { return sc; }
    std::ofstream& datfile() { return dat; }
    gnuplot(const std::string& aname) :
      name(aname), sc(scriptfn(aname)), dat(datfn(aname))
      {
         sc << "set terminal postscript eps color " << std::endl;
         sc << "set output '"
            << sys::outputfolder() + aname + ".eps" << "'" << std::endl;
      }
    void process()
    {
        sc.flush();
        if(!sc)
            throw "Error was writing to script " + scriptfn(name);
        dat.flush();
        if(!dat)
            throw "Error was writing to script " + datfn(name);
        std::string cmd = sys::get(klabel).second + klabel + " " + scriptfn(name);
        if(system(cmd.c_str()))
        {
            std::cerr << "Error executing " << cmd << std::endl;
            throw 1;
        }
    }

};





#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>


/// \version 0.7
/// \author Martin Smid
/** \copyright Unless otherwise specified, the code and the documentation was exclusively written by Martin Smid and may be used only under his permission.
*/

/** \mainpage Description

\section Introduction

EPP is a C++ library implementing various econometric computations. In the current stage of development it contains only a single file - \p EPP.HPP

\section Installation

NLopt
boost/numerics/ublas


*/


using namespace std;

const double na=-HUGE_VAL;
const double infinity=HUGE_VAL;


namespace tools
{

/// Routine computing matrix inversion

/**
 *  Taken from http://www.crystalclearsoftware.com/cgi-bin/boost_wiki/wiki.pl?LU_Matrix_Inversion
 *  Uses lu_factorize and lu_substitute in uBLAS to invert a matrix
 */

template<class T>
bool invermatrix (const boost::numeric::ublas::matrix<T>& input, boost::numeric::ublas::matrix<T>& inverse) {
    using namespace boost::numeric::ublas;
    typedef boost::numeric::ublas::permutation_matrix<size_t> pmatrix;
    // create a working copy of the input
    boost::numeric::ublas::matrix<T> A(input);
    // create a permutation matrix for the LU-factorization
    pmatrix pm(A.size1());

    // perform LU-factorization
    int res = lu_factorize(A,pm);
       if( res != 0 ) return false;

    // create identity matrix of "inverse"
    inverse.assign(boost::numeric::ublas::identity_matrix<T>(A.size1()));

    // backsubstitute to get the inverse
    lu_substitute(A, pm, inverse);

    return true;
}

} // namespace




struct paraminfo
{
    string name;
    double lower;
    double upper;
    double initial;
    paraminfo() :
      lower(-infinity), upper(infinity), initial(na)
      {}
    paraminfo(string aname, double ainitial) :
      name(aname), lower(-infinity), upper(infinity), initial(ainitial)
      {}
    paraminfo(string aname, double ainitial, double alower, double aupper) :
      name(aname), lower(alower), upper(aupper), initial(ainitial)
      {}
};

/// \p HUGE_VAL codes indefinedness

struct paramresult
{
    paraminfo info;
    double value;
    double std;
    double grad;
    double z() const
    {
        return (std==na || value==na || std==0)
                  ? na : value / std;
    }
    paramresult(const paraminfo& ainfo) :
       info(ainfo), value(na), std(na), grad(na)
       {}
    paramresult(const paramresult& r) :
       info(r.info), value(r.value), std(r.std), grad(r.grad)
       {}
    paramresult() {}
    static string stars(double z)
    {
        if(z==na)
            return "";
        double q = fabs(z);
        if(q>3.090)
            return "***";
        if(q>2.326)
            return "**";
        if(q>1.645)
            return "*";
        return "";
    }
    string stars() const
    {
        return stars(z());
    }
    void output(ostream& str, bool latex) const
    {
        static const char* natxt="n/a";
        if(latex)
            str << "$";
        str << info.name << "=";
        if(value==na)
            str << natxt;
        else
        {
            str <<  value  << "(";
            if(std==na)
                str << natxt << ")";
            else
            {
                str << std << ")";
                if(latex)
                    str << "^{" << stars() << "}";
                else
                    str << stars();
            }
            if(grad!=na)
                str << " g=" << grad;
        }
        if(latex)
            str << "$";
    }
};

class olatexstream : public ofstream
{
public:
   olatexstream(const string& afn) : ofstream(afn.c_str())
        {}
};


inline olatexstream& operator<<(olatexstream& str,const paramresult& r)
{
    r.output(str,true);
    return str;
};

inline ostream& operator<<(ostream& str,const paramresult& r)
{
    r.output(str,false);
    return str;
};


class object
{
public:
    enum elogginglevel { nologging, basiclogging, extendedlogging };
    elogginglevel logging;
protected:
    ofstream nulstream;
    ostream& lout;
protected:
    object() : logging(nologging), lout(cout)
    {}
    virtual ~object() {}
public:
    void setlogging(elogginglevel al = basiclogging)
    {
        logging = al;
    }
    bool islogging() { return logging == nologging ? false : true; }
 };

/// work in progress, neodszkousene

class summary: public object
{
    vector<unsigned int> hist;
protected:
    virtual unsigned int getN() const = 0;
    virtual double X(unsigned int i) const = 0;
public:
    const vector<unsigned int>& histogram(double binw,
                    double lower, double upper)
    {
        unsigned int nbins = (upper-lower) / binw + 0.5 + 2;
        hist.resize(nbins);
        for(unsigned int i=0; i<nbins; i++)
            hist[i]=0;
        unsigned int n=getN();
        for(unsigned int i=0; i<n; i++)
        {
            int j=(X(i)-lower)/binw+1;
            if(j<0)
                j = 0;
            else if((unsigned)j>= nbins)
                j=nbins - 1;
            hist[j]++;
        }
        return hist;
    }
};

class estimator: public virtual object
{
    vector<paraminfo> params;
protected:
    const vector<paraminfo>& getparams()
    {
        return params;
    }
public:
    estimator(const vector<paraminfo>& aparams) :
        params(aparams)
        {}

    ///  \p initialpars.size()==0 means no ninitial parameter
    virtual double estimate(vector<paramresult>& result,
                   bool stds = false) = 0;

};


class nloptuser: public virtual object
{
    nlopt::algorithm alg;
    bool ismax;
    double maxtime;
    double xtolrel;
    double ftolrel;
    double xtolabs;
    double ftolabs;
    bool gradient;

    string what;

    static double eval(const vector<double> &x,
                     vector<double> &g, void* instance)
    {
        nloptuser* self = (nloptuser*) instance;

        double res = na;
        self->what = "";
        try
        {
            res = self->objf(x,g);
        }
        catch(exception& e)
        {
            self->what = e.what();
            throw;
        }
        return res;
    }
protected:
    nloptuser(nlopt::algorithm aalg, bool amax):
       alg(aalg), ismax(amax), maxtime(HUGE_VAL),
       xtolrel(na), ftolrel(na), xtolabs(na), ftolabs(na),
       gradient(true)
       {}
    virtual double objf(const vector<double> &x,
                                  vector<double> &g) = 0;

    double opt(const vector<paraminfo>& params,
                             vector<paramresult>& result  )
    {
        unsigned int q = params.size();

        vector<double> x(q);
        vector<double> lower(q);
        vector<double> upper(q);

        for(unsigned int i=0; i<q; i++)
        {
            double in = params[i].initial;
            x[i] = in==na ? 0 : in;
            lower[i] = params[i].lower;
            upper[i] = params[i].upper;
        }
        nlopt::opt op( // tbd use alg
          gradient ? nlopt::LN_BOBYQA : nlopt::LN_BOBYQA , q); //LD_LBFGS
        if(xtolrel != na)
            op.set_xtol_rel(xtolrel);
        if(ftolrel != na)
            op.set_ftol_rel(ftolrel);
        if(xtolabs != na)
            op.set_xtol_abs(xtolabs);
        if(ftolabs != na)
            op.set_ftol_abs(ftolabs);

        op.set_maxtime(maxtime);
        op.set_lower_bounds(lower);
        op.set_upper_bounds(upper);
        if(ismax)
            op.set_max_objective(eval, this);
        else
            op.set_min_objective(eval, this);

        double r;
        nlopt::result optr;

        try
        {
           optr = op.optimize(x, r);
clog << "nlopt res: " << optr << endl;
           result.clear();
           for(unsigned int i=0; i<q; i++)
           {
               result.push_back(paramresult(params[i]));
               result[i].value = x[i];
           }
        }
        catch(...) // added in order to bridge "swallowing" oo
                    // exceptions by nlopt
        {
            if(what.length() > 0)
            {
                cerr << "Exception occured: " << what << endl;
                throw runtime_error(what);
            }
            else
                throw;
        }
        const char *msg;
        switch(optr)
        {
            case NLOPT_SUCCESS:
                msg = "Generic success return value.";
                break;
            case NLOPT_STOPVAL_REACHED:
                msg = "Optimization stopped because stopval was reached.";
                break;
            case NLOPT_FTOL_REACHED:
                msg = "Optimization stopped because ftol_rel or ftol_abs was reached.";
                break;
            case NLOPT_XTOL_REACHED:
                msg = "Optimization stopped because xtol_rel or xtol_abs was reached.";
                break;
            case NLOPT_MAXEVAL_REACHED:
                msg = "Optimization stopped because maxeval (above) was reached.";
                break;
            case nlopt::MAXTIME_REACHED:
                 throw runtime_error("nlopt: axtime reached");
            case nlopt::FAILURE:
                 throw runtime_error("nlopt: failure");
            case nlopt::INVALID_ARGS:
                 throw runtime_error("nlopt: invalid args");
            case nlopt::OUT_OF_MEMORY:
                 throw runtime_error("nlopt: out of memory");
            case nlopt::ROUNDOFF_LIMITED:
                 throw runtime_error("nlopt: roundoff error");
            case nlopt::FORCED_STOP:
                 throw runtime_error("nlopt: forced stop");
            default:
                msg = "Unknown nlopt return value";
        }
        if(islogging())
                lout << msg << endl;

        return r;
    }
public:
    void setmaxtime(double m)
    {
        maxtime = m;
    }
    void setftolrel(double f)
    {
        ftolrel = f;
    }
    void setxtolrel(double f)
    {
        xtolrel = f;
    }
    void setftolabs(double f)
    {
        ftolabs = f;
    }
    void setxtolabs(double f)
    {
        xtolabs = f;
    }
    void setgradient(bool g)
    {
        gradient = g;
    }
};

//class sample
//{
//protected:
//    virtual unsigned int getn() = 0;
//    virtual unsigned int getdim() = 0;
//
//};


/// Empirical distribution on {0,1,...}

class empdistn
{
    vector<unsigned int> dist;
    unsigned int n;
public:
    empdistn() : n(0) {}
    void add(unsigned int obs)
    {
        if(obs > dist.size())
        {
            unsigned int olds = dist.size();
            dist.resize(olds+10);
            for(unsigned int i=olds; i<dist.size(); i++)
                dist[i] = 0;
            dist[obs]++;
            n++;
        }
    }
    void getdist(vector<double> d)
    {
        if(n==0)
        {
            d.resize(0);
        }
        else
        {
            d.resize(dist.size());
            for(unsigned int i=0; i<dist.size(); i++)
                d[i] = (double) dist[i] / (double) n;
        }
    }
};

/// Class computing MLE estinate


class mle : public estimator, public nloptuser
{
    // state variables
    boost::numeric::ublas::matrix<double> J;
    bool computeJ;

    double objf(const vector<double> &x, vector<double> &g)
    {

        unsigned int R = x.size();
        beforeloglikeval(x);
        for(unsigned int j=0; j<g.size(); j++)
            g[j] = 0;
        boost::numeric::ublas::matrix<double>
          I(boost::numeric::ublas::zero_matrix<double>(R,R));
        double s=0;

        unsigned int M = getN();
        for(unsigned int i=0; i<M; i++)
        {
            vector<double> grad(R);
            s += evallogdensity(i,x,grad);
            for(unsigned int j=0; j<g.size(); j++)
                g[j] += grad[j];
            if(computeJ)
                for(unsigned int j=0; j<R; j++)
                {
                    for(unsigned int k=0; k<R; k++)
                        I(j,k)+=grad[j]*grad[k] / (double) M;
                }
        }
        afterloglikeval(x, g, s);
        if(computeJ)
        {
            try
            {
                tools::invermatrix(I, J);
            }
            catch(...)
            {
                if(islogging())
                {
                    lout << "Cannot invert matrix" << endl;
                    for(unsigned int i=0; i<R; i++)
                    {
                        for(unsigned int j=0; j<R; j++)
                            lout << I(i,j) << ",";
                        lout << endl;
                    }
                }
                for(unsigned int i=0; i<R; i++)
                    for(unsigned int j=0; j<R; j++)
                        J(i,j) = na;
            }
        }
        if(logging >= extendedlogging)
        {
            lout << "#" << " LL=" << s << " args=(";
            for(unsigned int j=0;j<x.size();j++)
                lout << x[j] << " ";
            lout << ")";
            if(g.size())
            {
               lout << " grad=(";
               for(unsigned int j=0;j<g.size();j++)
                  lout << g[j] << " ";
               lout << ")";
            }
            lout << endl;
        }
//clog << "l(" << xx[0]<<","<< xx[1]<<"," << xx[2] << ") "
//        << "-l(" << x[0]<<","<< x[1]<<"," << x[2] << ") "
//        << ss << "-" << s << "..." << (ss-s) / 0.0001 << "=mlr=" << g[2] ;
//clog << " returning "<< s << endl;
        return s;
    }

protected:
    /// called once before repeated evaluation of evallogdensity (with the same params but different i)
    virtual int getN() = 0;
    virtual void beforeloglikeval(const vector<double> aparams) {}
    virtual double evallogdensity(int i, const vector<double>& aArgs,
                                                  vector<double>& g) = 0;
    virtual void afterloglikeval(const vector<double>& aparams,
                         const vector<double>& grad, double loglik) {}
public:

    mle(const vector<paraminfo>& aparams):
        estimator(aparams), nloptuser(nlopt::LD_LBFGS, true),
        J(aparams.size(),aparams.size()),
        computeJ(false)
      {}

    ///  \p initialpars.size()==0 means no ninitial parameter
    double estimate(vector<paramresult>& result,
                    bool stds = false)
    {
        unsigned int R=getparams().size();
        if(islogging())
        {
            lout << "MLE: Estimating " << R << " parameters using "
                << getN() << " observations." << endl;
        }

        computeJ = false;

        double loglik = opt(getparams(),result);

        if(stds)
        {
            computeJ = true;
            vector<double> x(R);
            vector<double> g(R);

            for(unsigned int i=0; i<R; i++)
                x[i] = result[i].value;

            loglik = objf(x,g);

            for(unsigned int i=0; i<R; i++)
            {
                double K = J(i,i);
                if(K!=na)
                    result[i].std=sqrt(K / (double) getN());
                result[i].grad = g[i];
            }
        }
        if(islogging())
        {
            lout << "Optimal params:" << endl;
            for(unsigned int j=0; j<result.size(); j++)
                lout << result[j] << endl;
            lout << endl;
        }

        return loglik;
    }

    vector<vector<double>> getJ()
    {
        unsigned int R=getparams().size();
        vector<vector<double>> res(R,vector<double>(R));
        for(unsigned i=0; i<R; i++)
            for(unsigned j=0; j<R; j++)
                res[i][j]=J(i,j);
        return res;
    }

/*	virtual string description()
    {
        string s = "MLE";
        return s;
    }*/
/*	void result(ostream& ofs, mle* refm = 0, bool tex = false)
    {
        using namespace std;
        if(!estimated)
            throw 1;

        ostream& t = ofs;

        if(!tex)
            ofs << endl << description() << " n=" << getN() << " r=" << getR() << endl << endl;

        const boost::numeric::ublas::matrix<double> Jm = J();

        for(unsigned int i=0; i<getR(); i++)
        {
            string s;
            double p = params[i];
            double sd = sqrt(Jm(i,i) / (double) getN());
            if(!tex)
                ofs << getlabel(i) << " " << p << "(" << sd << ") ";

            double z = fabs(p)/sd;

            if(!tex)
                ofs << z;

            string stars;
            if(z > 1.96)
                stars+= "*";
            if(z > 2.576)
                stars+= "*";
            if(z > 3.29)
                stars+= "*";

            if(!tex)
                ofs << stars << endl;
//			*mlog << m.getlabel(i) << "," << p << "," << sd
//					<< "," << z << "," << stars << endl;
            if(tex)
            {
                t << "$" << getlabel(i) << resetiosflags( ios::floatfield )
                  << "$ = $" << p << "$ ($" << sd << "$)" << stars << "\\\\" << endl;
            }
        }

        if(!tex)
            ofs << endl <<  "loglik=" << loglik() << endl;
    //    *mlog << endl <<  "loglik," << loglik << endl;


    //	*mlog << "rho," << rho << endl;

        double refll;
        int refdf;
        string rd;
        double rho;
        if(refm)
        {
            refll = refm->loglik();
            refdf = getR()-refm->getR();
            rd = refm->description();
            rho = 1.0 - loglik() / refm->loglik();
        }
        else
        {
            refll = loglik();
            refdf = getR();
            rd = "zero parameters";
            rho = 0;
        }
        if(!tex)
            ofs << "rho = " << rho << endl;

        if(!tex)
            ofs << "ratio = " << -2*(refll - loglik()) << " df=" << refdf << endl << endl;
        if(0 && tex)
        {
            t << "\\hline" << endl;
                        t << "\\end{tabular}" << endl << endl;
                        t << "\\vspace{5mm}" << endl << endl;

            t << "\\begin{tabular}{lrlr}" << endl;
            t << "$\\rho$ & " << rho << " & observations &" << getN() << "\\\\" << endl;
            t << "likelihood ratio & " << -2*(refll - loglik()) << " & d.f & "
                    << refdf << "\\\\" << endl;
            t << "\\end{tabular}" << endl;
        //    t << "%(" << m.description() << ") vs (" << rd << ")" << endl;
            t << "\\end{center}" << endl<< endl<< endl<< endl;
        }
    }*/

    double loglik(const vector<double> &x, vector<double> &g)
    {
        return objf(x,g);
    }

};

/// A class encapsulting calls of GRETL

class gretl
{
    string dir;
    string tmpfn()
    {
        return dir+"/_gretl_tmp";
    }
    string scriptresult;
public:
    gretl(const string& adir): dir(adir) {}
    string scdir()
    {
        return dir;
    }
    void runscript(const string& scfn)
    {
        using namespace std;
        string cmd = "gretlcli -b \"" + scfn + "\" > " + tmpfn();
        cout << cmd << endl;
        if(system(cmd.c_str()))
        {
            cerr << "Failed to run " << cmd << endl;
            throw 1;
        }
        ifstream p(tmpfn().c_str());
        if(!p)
        {
            cerr << "Error opening " << tmpfn() << endl;
            throw 1;
        }

        getline(p,scriptresult,'\0');
    }
    string getscriptoutput()
    {
        ifstream s(tmpfn().c_str());
        if(!s)
        {
                        cerr << "Error opening " << tmpfn().c_str() << endl;
            throw 1;
        }
        string str;
        getline(s,str);
        return str;
    }
    double findpar(const string& parname) // assumes that "? parname" appears in the output
    {
        string label = "? " + parname;
        unsigned int pos;
        if((pos = scriptresult.find(label))==string::npos)
        {
            cerr << "Cannot find parameter " << parname << " in gretl script output " << endl;
            throw 1;
        }
        while(scriptresult[pos++] != '\n')
            ;
        return atof(scriptresult.c_str()+pos);
    }
};

/// A class handling CSV input

class csvRow
{
    public:
        string const& operator[](size_t index) const
        {
            return m_data[index];
        }
        size_t size() const
        {
            return m_data.size();
        }
        void readNextRow(istream& str)
        {
            string         line;
            getline(str,line);

            stringstream   lineStream(line);
            string         cell;

            m_data.clear();
            while(getline(lineStream,cell,','))
            {
                m_data.push_back(cell);
            }
        }
    private:
        vector<string> m_data;
};

/// Operator reading a single csv row
/** TODO: handle end of line
*/

inline istream& operator>>(istream& str,csvRow& data)
{
    data.readNextRow(str);
    return str;
}


inline vector<double> operator+ (const vector<double>& x, const vector<double>& y)
{
    unsigned int q=x.size();
    if(y.size() != q)
        throw runtime_error("Different sizes in vector + oparation");
    vector<double> r(q);
    for(unsigned int s=0; s<q; s++)
        r[s] = x[s] + y[s];
    return r;
}



inline vector<double> operator- (const vector<double>& x, const vector<double>& y)
{
    unsigned int q=x.size();
    vector<double> r(q);
    for(unsigned int s=0; s<q; s++)
        r[s] = x[s] - y[s];
    return r;
}


inline vector<double> operator* (double c, const vector<double>& x)
{
    unsigned int q=x.size();
    vector<double> r(q);
    for(unsigned int s=0; s<q; s++)
        r[s] = c *x[s];
    return r;
}

inline vector<double> operator* (const vector<double>& y, double c)
{
    return c*y;
}

inline vector<double> operator/ (const vector<double>& y, double c)
{
    return y*(1.0/c);
}


inline vector<double> const_vector(double c, size_t q)
{
    vector<double> r(q);
    for(unsigned int s=0; s<q; s++)
        r[s] = c;
    return r;
}



inline vector<double> zero_vector(size_t q)
{
    return const_vector(0,q);
}


/*
class ols
{
    vector<double> **X;
    int r;
    const vector<double>& Y;
    bool c;
    bool estimated;
    boost::numeric::ublas::vector<double> bresult;
    vector<double> mhatY;
    vector<double> mhate;
    double mSSE;
public:
    ols(vector<double>** aX, int ar, vector<double>& aY, bool aConst = true):
        X(aX), r(ar), Y(aY), c(aConst), estimated(false),
        bresult(ar + (aConst?1:0)),
        mhatY(aY.size()), mhate(aY.size()) {}
    double b(int i) { if(estimated) return bresult[i]; else throw 1;}
    const vector<double>& hatY() { if(estimated) return mhatY; else throw 1;}
    const vector<double>& hate() { if(estimated) return mhate; else throw 1;}
    double SSE() { if(estimated) return mSSE; else throw 1;};
    void estimate()
    {
        int cshift = c ? 1 : 0;
        using namespace boost::numeric::ublas;
        int k=bresult.size();
        int n=Y.size();
        matrix<double> XtX(k,k);
        boost::numeric::ublas::vector<double> XtY(k);
        for(int i=0; i<k; i++)
        {
            double t=0;
            for(int h=0; h<n; h++)
            {
                int ii = i-cshift;
                double x = ii < 0 ? 1 : (*X[ii])[h];
                t +=x*Y[h];
            }
            XtY[i] = t;

            for(int j=0; j<k; j++)
            {
                double s=0;
                for(int h=0; h<n; h++)
                {
                    int ii = i-cshift;
                    double x1 = ii < 0 ? 1 : (*X[ii])[h];

                    int jj = j-cshift;
                    double x2 = jj < 0 ? 1 : (*X[jj])[h];

                    s += x1*x2;
                }
                XtX(i,j)=s;
            }
        }

        matrix<double> XtXinv(k,k);
        try
        {
            tools::InvertMatrix(XtX,XtXinv);
        }
        catch(...)
        {
            cerr << "unable to invert: " << endl;
            cerr << XtX << endl;
            throw;
        }
        bresult = prod(XtXinv, XtY);

        mSSE = 0;
        for(int i=0; i<n; i++)
        {
            double y = c ? bresult[0] : 0;
            for(int j=0; j<r; j++)
                y += bresult[j+cshift] * (*X[j])[i];
            mhatY[i] = y;
            double e = Y[i] - y;
            mhate[i] = e;
            mSSE += e * e;
        }
        estimated = true;
    }
};
*/

class uncertain
{
public:
    uncertain(const vector<double>& x, const matrix& var) : fx(x), fvar(var)
        {}
    uncertain(const vector<double>& x) : fx(x), fvar(matrix(x.size(),x.size(),0))
        {}
    uncertain() : fx(0), fvar(matrix(0,0,0))
        {}
    vector<double> x() const { return fx;}
    unsigned dim() const { return fx.size(); }
    matrix var() const { return fvar; }
private:
    vector<double> fx;
    matrix fvar;
};

} // namespace

#endif // FUTUREORPP_HPP
