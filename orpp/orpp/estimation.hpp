#ifndef ESTIMATION_HPP
#define ESTIMATION_HPP

#include <string>
#include "orpp/matrices.hpp"

namespace orpp
{

class uncertain
{
public:
    uncertain(const dvector& x, const dmatrix& var) : fx(x), fvar(var)
        {}
    uncertain(const dvector& x) : fx(x), fvar(dmatrix(x.size(),x.size()))
        {  fvar.setZero(); }
    uncertain(const std::vector<double>& x) : fx(dv(x)), fvar(dmatrix(x.size(),x.size()))
        {  fvar.setZero(); }
    uncertain() : fx(0), fvar(dmatrix(0,0))
        {}
    dvector x() const { return fx;}
    dvector& x() { return fx;}
    unsigned dim() const { return fx.size(); }
    dmatrix var() const { return fvar; }
    dmatrix& var()  { return fvar; }
private:
    dvector fx;
    dmatrix fvar;
};


struct paraminfo
{
    std::string name;
    double lower;
    double upper;
    double initial;
    paraminfo() :
      lower(-infinity), upper(infinity), initial(na)
      {}
    paraminfo(std::string aname, double ainitial) :
      name(aname), lower(-infinity), upper(infinity), initial(ainitial)
      {}
    paraminfo(std::string aname, double ainitial, double alower, double aupper) :
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
    static std::string stars(double z)
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
    std::string stars() const
    {
        return stars(z());
    }
    void output(std::ostream& str, bool latex) const
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

class olatexstream : public std::ofstream
{
public:
   olatexstream(const std::string& afn) : std::ofstream(afn.c_str())
        {}
};


inline olatexstream& operator<<(olatexstream& str,const paramresult& r)
{
    r.output(str,true);
    return str;
};

inline std::ostream& operator<<(std::ostream& str,const paramresult& r)
{
    r.output(str,false);
    return str;
};


} // namespace

#endif // ESTIMATION_HPP
