#ifndef OPTIMIZATION_HPP
#define OPTIMIZATION_HPP

#include "orpp.hpp"

namespace orpp
{

/// \addtogroup pvar Variables
/// @{

/// Real decision variable
using realvar=double;
/// Integer decision variable
using intvar=long;
/// Binary decision variable
using binvar=bool;
/// Mixed decision variable
union mixedvar
{
    realvar x;
    intvar n;
    binvar b;
    mixedvar(const realvar& ax) : x(ax) {}
    mixedvar(const intvar& an) : n(an) {}
    mixedvar(const binvar& ab) : b(ab) {}
    operator realvar() const { return x; }
    operator intvar() const { return n; }
    operator binvar() const { return b; }
//    mixedvar& operator =(realvar& v) { x=v; return *this; }
//    mixedvar& operator =(intvar& v) { n=v; return *this; }
};


/// Vector of variables of type \p V
template <typename V>
using variables=std::vector<V>;


/// Range (interval of possible values) of variables \p V
template <typename V>
class range : public object
{
public:
    using V_t = V;

    enum type { realt, intt, bint };

    range()
    {
        if constexpr( std::is_same<V,realvar>::value)
           ft = realt;
        else if constexpr( std::is_same<V,intvar>::value)
           ft = intt;
        else if constexpr( std::is_same<V,binvar>::value)
           ft = bint;
        else
        {
           assert(0);
        }
        setlimits();
    }


    range(type t, V l=minfinity<V>(), V h=infinity<V>())
    {
        if constexpr( std::is_same<V,realvar>::value)
        {
           assert(t==realt);
           ft = realt;
        }
        else if constexpr( std::is_same<V,intvar>::value)
        {
            assert(t==intt);
            ft = intt;
        }
        else if constexpr( std::is_same<V,binvar>::value)
        {
            assert(t==bint);
            ft = bint;
        }
        else if constexpr( std::is_same<V,mixedvar>::value)
        {
           assert(t==realt || t==intt || t==bint);
           ft = t;
        }
        else
            assert(0);
        setlimits(l,h);
    }

/*    void setreal(realvar l=minf<realvar>(), realvar h=inf<realvar>())
    {
        static_assert(std::is_same<V,mixedvar>::value);
        ft = realt;
        setlimits(l,h);
    }

    void setint(intvar l=minf<intvar>(), intvar h=inf<intvar>())
    {
        static_assert(std::is_same<V,mixedvar>::value);
        ft = intt;
        setlimits(l,h);
    }
*/
    void setlimits(V l=minfinity<V>(), V h=infinity<V>())
    {
        assert(ft==realt || ft==intt);
        fl=l; fh=h;
    }

    void setpositive()
    {
       setlimits(0);
    }
    V l() const
    {
        assert(ft==realt || ft==intt);
        return fl;
    }
    V h() const
    {
        assert(ft==realt || ft==intt);
        return fh;
    }
    type t() const
    {
        assert(ft>= realt && ft <= bint);
        return ft;
    }
    bool ishinf() const
    {
        assert(ft==realt || ft==intt);
        if(ft==realt)
        {
            realvar h = fh;
            return h==infinity<realvar>();
        }
        else
        {
            intvar h = fh;
            return h==infinity<int>();
        }
    }
    bool islinf() const
    {
        assert(ft==realt || ft==intt);
        if(ft==realt)
        {
            realvar l = fl;
            return l==-infinity<realvar>();
        }
        else
        {
            intvar l = fl;
            return l==minfinity<int>();
        }
    }
private:
    V fl;
    V fh;
    type ft;
};

// Container of decision variables ranges
template <typename V>
class ranges
{
public:
    ranges(unsigned int n) : frs(n)
    {
    }
    range<V>& operator [] (unsigned int i)
    {
        assert(i<frs.size());
        return frs[i];
    }
    operator std::vector<range<V>> () const
    {
        return frs;
    }
    unsigned int size() const { return frs.size(); }
private:
    std::vector<range<V>> frs;
};


///@}


/// \addtogroup problems Problems
/// @{

/// Base class of decision problems' constraints
class constraint: public object
{
public:
    enum type {eq, geq, leq};
    constraint(unsigned int xdim, type t=eq)
        : ft(t), fxdim(xdim) {}

    type t() const { return ft; }
    void settype(type t) { ft = t; }
    unsigned int xdim() const { return fxdim; }
private:
    type ft;
    unsigned int fxdim;
};

const constraint::type eq = constraint::eq;
const constraint::type leq = constraint::leq;
const constraint::type geq = constraint::geq;

/// \brief Vector-like container of constraints
/// \tparam G constraint type
template <typename G>
class msconstraints
{
public:
    msconstraints(unsigned int n) : fn(n) {}
    G& add(const G& g)
    {
        assert(g.xdim()==fn);
        fgs.push_back(g);
        return *fgs.rbegin();
    }
    G& add()
    {
        fgs.push_back(G(fn));
        return *fgs.rbegin();
    }
    G& operator [] (unsigned int i)
    {
        assert(i<fgs.size());
        return fgs[i];
    }

    operator std::vector<G> () const
    {
        return fgs;
    }

private:
    unsigned int fn;
    std::vector<G> fgs;
};


/// Base class for decision criterions
class criterion: public object
{
};

///@}

/// \addtogroup constraints Constraints
/// @{


class linearconstraint : virtual public constraint
{
public:
    linearconstraint(unsigned int xdim) :
        constraint(xdim), flhs(xdim), frhs(0) {}
    linearconstraint(const std::vector<double>& lhs, double rhs=0,
            constraint::type t=constraint::eq) :
        constraint(lhs.size(),t), flhs(lhs), frhs(rhs) {}
    double lhs(unsigned int i) const
    {
        assert(i<flhs.size());
        return flhs[i];
    }
    void setlhs(unsigned int i, double v)
    {
        assert(i<flhs.size());
        flhs[i]=v;
    }
    void setlhs(const std::vector<double>& lhs)
    {
        assert(lhs.size()==xdim());
        flhs = lhs;
    }
    void setrhs(double rhs ) { frhs = rhs; }
    double rhs() const { return frhs; }
private:
    std::vector<double> flhs;
    double frhs;
};

///@}


} // namespace

#endif // OPTIMIZATION_HPP
