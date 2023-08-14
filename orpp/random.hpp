#ifndef RANDOM_H
#define RANDOM_H

#include "orpp.hpp"
#include <cmath>
#include <algorithm>

namespace orpp
{

/// \addtogroup General Random
/// @{

enum esignificance { enotsignificant, e10percentsignificant, e5percensignificant,
                   e1percentsignificant, enumsignificances};

template <typename E>
struct valuewitherror
{
    E x;
    double sd;
};

struct realestimate : public valuewitherror<double>
{
    esignificance twosidedsignificance() const
    {
        assert(sd > 0);
        double s = fabs(x / sd);
        if(s > 2.5758)
            return e1percentsignificant;
        else if(s > 1.96)
            return e5percensignificant;
        else if(s > 1.6449)
            return e10percentsignificant;
        else
            return enotsignificant;
    }
};


inline std::string stars(esignificance s)
{
    switch(s)
    {
    case enotsignificant:
        return "";
    case e10percentsignificant:
        return "*";
    case e5percensignificant:
        return "**";
    case e1percentsignificant:
        return "***";
    default:
        throw exception("unknonw significance");
    }
}

inline realestimate operator+(const realestimate& a, realestimate& b)
{
    return { a.x + b.x, sqrt( a.sd * a.sd + b.sd * b.sd)};
}

inline realestimate operator-(const realestimate& a, realestimate& b)
{
    return { a.x - b.x, sqrt( a.sd * a.sd + b.sd * b.sd)};
}


/// @} - General Random


/// \addtogroup Distributions
/// @{

constexpr double probabilitytolerance = 1e-6;

/// \brief Base class.
///
///  \tparam I type of the values
///  \tparam C type of the condition
template <typename I=double, typename C=nothing>
class distribution: public object
{
public:
    distribution() {}
    using I_t = I;
    using C_t = C;
};

template <typename I=double>
class markovdistribution: public distribution <I,I>
{
public:
    markovdistribution() {}
    using I_t = I;
    using C_t = I;
};


/// \brief Joint distribution
/// \tparam F type of the first coordinate
/// \tparam S type of the second coordinate
/// \tparam C type of the condition
template <typename F, typename S, typename C>
class jdistribution: virtual public distribution
        <std::pair<F, S>, C>
{
public:
    using F_t = F;
    using S_t = S;
};


/// \brief std::vector distribution
/// \tparam X type of the std::vector coordinates
/// \tparam C type of the condition
template <typename X, typename C>
class vdistribution : virtual public distribution<std::vector<X>,C>
{
public:
    using X_t = X;
    unsigned int dim() const { return this->dim_is(); }
private:
    virtual unsigned int dim_is() const = 0;
};


/// \addtogroup mdisc MC distributions
/// \ingroup Distributions
/// @{

/// \brief Monte Carlo distribution
/// \tparam I type of the values
/// \tparam C type of the condition
///
/// Random numbers may be generated from the distribution.

template <typename I, typename C>
class mdistribution: virtual public distribution<I, C>
{
public:
    I draw(const C& c) const
    {
        return do_draw(c);
    }
    I draw() const
    {
        static_assert(std::is_same<C,nothing>::value);
        return do_draw(na);
    }
private:
    virtual I do_draw(const C& c) const = 0;
};

/// @}

/// \addtogroup realdists Real distributions
/// \ingroup Distributions
/// @{

template <typename C>
class realdistribution: virtual public distribution<double, C>
{
public:
    double mean(const C& c, double prec=0.0) const
    {
        return mean_is(c,prec);
    }
    double mean(double prec=0.0) const
    {
        return mean(na, prec);
    }
protected:
    virtual double mean_is(const C&, double  ) const
    {
        throw orpp::exception("Mean undefined!");
    }
};

template <typename C>
class cdfdistribution: virtual public realdistribution<C>
{
public:

    probability cdf(double x) const
    {
        return this->cdf(x,na);
    }
    probability cdf(double x, const C& c) const
    {
        if(x <= orpp::min<double>())
            return 0;
        else if(x >= orpp::max<double>())
            return 1;
        else
            cdf_is(x,c);
    }
private:
    virtual probability cdf_is(double, const C& c) const = 0;
};


template <typename C>
class qdistribution: virtual public cdfdistribution<C>,
        virtual public mdistribution<double,C>
{
public:
    double quantile(probability p) const
    {
        assert(p<=1);
        assert(p>=0);
        return this->quantile(p,na);
    }
    double quantile(probability p, const C& c) const
    {        
        return quantile_is(p,c);
    }
private:
    virtual double quantile_is(probability, const C&) const = 0;
       // a numeric method may be implemnted here
    virtual double do_draw(const C& c) const
    {
        probability u = sys::uniform();
        return quantile(u,c);
    }
protected:
    virtual double mean_is(const C& c, double prec) const
    {
        assert(prec>0);
        double s=0;
        double s2=0;
        unsigned int i=1;
        for(;;i++)
        {
            double x = this->draw(c);
            s += x;
            s2 += x*x;
            if(i%100 == 0)
            {
                double var = s2/i - (s/i) *(s/i);
                if(sqrt(var/i)<prec)
                    break;
            }
        }
        return s / i;
    }
};

template <typename D>
class scaleddistribution : public qdistribution<nothing>
{
    void check()
    {
        assert(fsd > 0);
        static_assert(std::is_base_of<qdistribution<nothing>,D>::value);
    }

public:
    scaleddistribution(double m, double sd) : fm(m), fsd(sd)
    {
        check();
    }

    scaleddistribution(const D& d, double m, double sd)
        : fm(m), fsd(sd), fd(d)
    {
        check();
    }
    /// returns the distribution which whas scalled
    const D& srcd() const
    {
        return fd;
    }
    double m() const { return fm; }
    double sd() const { return fsd; }
private:
    virtual double quantile_is(probability p, const nothing&) const
    {
        return fm + fsd*fd.quantile(p);
    }
    virtual probability cdf_is(double x, const nothing&) const
    {
        return fd.cdf((x-fm) / fsd);
    }

    virtual double mean_is(const nothing&, double prec) const
    {
        return (fd.mean(prec / fsd) + fm) * fsd;
    }

    double fm;
    double fsd;
    D fd;
};


/// not tested yet
template <typename D>
class truncateddistribution : public qdistribution<nothing>
{
    void check()
    {
        assert(fh > fl);
        static_assert(std::is_base_of<qdistribution<nothing>,D>::value);
    }

public:
    truncateddistribution(double l, double h) : fl(l), fh(h)
    {
        check();
    }
    truncateddistribution(double l, double h, const D& d)
        : fl(l), fh(h), fd(d)
    {
        check();
    }

    /// returns the distribution which whas scalled
    const D& srcd() const { return fd; }
    double l() const { return fl; }
    double h() const { return fh; }
private:
    virtual double quantile_is(probability p, const nothing&) const
    {
        double arg = p * fd.cdf(fh) + (1-p) * fd.cdf(fl);
        assert(arg);
        return fd.quantile(arg);
    }
    virtual probability cdf_is(double x, const nothing&) const
    {
        double cl = fd.cdf(fl);
        double denom=fd.cdf(fh)-cl;
        assert(denom>0);
        return (fd.cdf(x)-cl)/denom;
    }

    double fh;
    double fl;
    D fd;
};



/// tbd own draw
template <typename D>
class arqdistribution : public qdistribution<double>
{
    void check()
    {
        static_assert(std::is_base_of<qdistribution<nothing>,D>::value);
    }
public:
    arqdistribution(double a=1) : fa(a)
    {
        check();
    }

    arqdistribution(const D& d, double a=1) : fa(a), fd(d)
    {
        check();
    }
    double a() const { return fa; }
    const D& srcd() const { return fd; }
private:
    virtual double quantile_is(probability p, const double& c) const
    {
        return fa * c + fd.quantile(p);
    }
    virtual probability cdf_is(double x, const double& c) const
    {
        return fd.cdf(x-fa * c);
    }

    double fa;
    D fd;
};

/// multiplied distribution
/// tbd own draw
template <typename D>
class multiplieddistribution : public qdistribution<double>
{
    void check()
    {
        assert(fb>0);
        static_assert(std::is_base_of<qdistribution<nothing>,D>::value);
    }
public:
    multiplieddistribution(double b=1) : fb(b)
    {
        check();
    }

    multiplieddistribution(const D& d, double b=1) : fb(b), fd(d)
    {
        check();
    }
    double b() const { return fb; }
    const D& srcd() const { return fd; }
private:
    virtual double quantile_is(probability p, const double& c) const
    {
        assert(c>0);
        return fb * c * fd.quantile(p);
    }
    virtual probability cdf_is(double x, const double& c) const
    {
        assert(c>0);
        return fd.cdf(x/(fb * c));
    }

    double fb;
    D fd;
};


/// @} // real distributions


/// \addtogroup discretedisc Discrete distributions
/// \ingroup Distributions
/// @{



/// \brief Atom of #ddistribution
/// \tparam I type of the values
template <typename I>
struct atom
{
    static bool comparator(const atom<I>& l, const atom<I>& r)
    {
        return  (l.x < r.x);
    }

    using I_t = I;
    I x;
    probability p;
};

/// \brief Discrete distribution
/// \tparam I type of the values
/// \tparam C type of the condition
template <typename I, typename C>
class ddistribution: virtual public mdistribution<I,C>
{
public:
    atom<I> operator () (unsigned int i, const C& c) const
    {
        return atom_is(i,c);
    }
    atom<I> operator () (unsigned int i) const
    {
        static_assert(std::is_same<C,nothing>::value);
        return atom_is(i,na);
    }
    bool sorted() const { return is_sorted(); }
private:
    /// atoms must be sorted TBD
    virtual atom<I> atom_is(unsigned int i, const C& c) const  = 0;
    virtual I do_draw(const C& c) const
    {
        return this->drawuniv(c);
    }
protected:
    virtual bool is_sorted() const = 0;
    /// generates an uniform variable \p u and goes through the probabilities
    /// until their sum exceeds \p u. May be very slow
    I drawuniv(const C& c) const
    {
        if(!is_sorted())
            throw exception("distribution must be sorted to use drawuniv");
        double u = sys::uniform();
        double sum = 0;
#ifndef NDEBUG
        I last;
#endif
        for(unsigned int i=0; i<std::numeric_limits<unsigned int>::max() ; i++)
        {
            atom<I> a = (*this)(i,c) ;
            sum += a.p;
#ifndef NDEBUG
            assert(i == 0 || a.x >= last);
            last=a.x;
#endif
            if(sum >= u)
                return a.x;
        }
        assert(0);
    }
};

/// \brief Finite discrete distribution
/// \tparam I type of the values
/// \tparam C type of the condition
/// \tparam listdef vay of atoms definition:
///
/// \p listdef=true indicates
/// that the atoms are primarily determined by an overloading of
/// \p atoms_are while \p atom_is and \p natom_is need not be overoaded
/// as they use \p atoms_are. If
/// \p listdef=false then the atoms are
/// determided by overloading of \p natom_is and \p atom_is while \p atoms_are needs not
/// be overloaded.
template <typename I, typename C /*, bool listdef = false */ >
class fdistribution: virtual public ddistribution<I,C>
{
public:
//    static bool constexpr flistdef=listdef;
//    template <bool dosort=false>
    void atoms(std::vector<atom<I>>& a, const C& c) const
    {
        a.resize(0); // better safe than sorry
        index n = natoms_is(c);
        for(index i=0; i<n; i++)
        {
            a.push_back({atom_is(i,c)});
        }

        probability p=0;
        #ifndef NDEBUG
            for(unsigned int i=0; i< a.size(); i++)
                p += a[i].p;
            assert(fabs(p-1.0) < probabilitytolerance);
        #endif
    }

    void atoms(std::vector<atom<I>>& a) const
    {
        static_assert(std::is_same<C,nothing>::value);
        return atoms(a,na);
    }

    unsigned int natoms(const C& c) const
    {
        return natoms_is(c);
    }

    unsigned int natoms() const
    {
        static_assert(std::is_same<C,nothing>::value);
        return natoms_is(na);
    }
    double mean(const C& c, double prec=0) const
    {
        static_assert(std::is_convertible<I,double>::value);
        std::vector<atom<I>> a;
        this->atoms(a,c);
        double sum = 0;
        for(unsigned int i=0; i<a.size(); i++)
            sum += a[i].p * a[i].x;
        return sum;
    }
private:
    virtual unsigned int natoms_is(const C& c) const = 0;
    virtual atom<I> atom_is(unsigned int i, const C& c) const = 0;
    /// returning \p false means that slow \ref drawuniv will be used in MC.

    virtual bool is_equiprobable(const C& c) const { return false; }
    virtual I do_draw(const C& c) const
    {
        if(is_equiprobable(c))
        {
            double u = sys::uniform();
            unsigned int n = natoms(c);
            unsigned int i = static_cast<unsigned int>(floor(n*u));
            return (*this)(i,c).x;
        }
        else
            return this->drawuniv(c);
    }
};


/// \brief List defined distribution
/// \tparam I type of the values
///
///
template <typename I, bool sortable=false>
class ldistribution: public fdistribution<I,nothing>
{   

public:
    ldistribution(const std::vector<atom<I>>& atoms, bool equiprobable, bool sorted) :
        fatoms(atoms), fequiprobable(equiprobable), fsorted(sorted)
    {
        assert(fatoms.size());
#ifndef NDEBUG
        probability tp = fatoms[0].p;
        double sum = tp;
        auto last = fatoms[0].x;
        for(unsigned int i=1; i<fatoms.size(); i++)
        {
            auto p = fatoms[i].p;
            sum += p;
            if(tp != p && equiprobable)
                throw exception("ldistribution initialized by non-equiprobable data");
            if constexpr(sortable)
            {
                if(sorted && fatoms[i].x < last)
                   throw exception("ldistribution initialized by non-sorted data");
                last = fatoms[0].x;
            }
        }
        if(fabs(sum - 1.0) > probabilitytolerance)
            throw exception("ldistribution does not sum to one");
#endif
    }

    ldistribution(const std::vector<I>& values, bool sorted) :
        fatoms(values.size()), fequiprobable(true), fsorted(sorted)
    {
        assert(values.size());
        probability p=1.0/fatoms.size();
        for(unsigned int i=0; i<fatoms.size(); i++)
        {
            fatoms[i].x = values[i];
            fatoms[i].p = p;
        }
#ifndef NDEBUG
        auto last = fatoms[0].x;
        for(unsigned int i=1; i<fatoms.size(); i++)
        {
            if constexpr(sortable)
            {
                if(sorted && fatoms[i].x < last)
                   throw exception("ldistribution initialized by non-sorted data");
                last = fatoms[0].x;
            }
        }
#endif

    }
private:

    virtual unsigned int natoms_is(const nothing&) const
    {
        return fatoms.size();
    }

    virtual atom<I> atom_is(unsigned int i, const nothing&) const
    {
        return fatoms[i];
    }
    std::vector<atom<I>> fatoms;
    bool fequiprobable;
    bool fsorted;
    virtual bool is_sorted() const { return fsorted; }
    virtual bool is_equiprobable(const nothing&) const { return fequiprobable; }
};


template <typename I, bool sortable=false>
class equipdistribution: public fdistribution<I,nothing>
{

public:
    equipdistribution(const std::vector<I>& values, bool sorted) :
        fvalues(values), fsorted(sorted)
    {
/*#ifndef NDEBUG
        auto last = values[0];
        for(unsigned int i=1; i<values.size(); i++)
        {
            if constexpr(sortable)
            {
                if(sorted && values[i] < last)
                   throw exception("ldistribution initialized by non-sorted data");
                last = values[0];
            }
        }
#endif*/

    }
protected:
    std::vector<I>& values()
    {
        return fvalues;
    }

private:

    virtual unsigned int natoms_is(const nothing&) const
    {
        return fvalues.size();
    }

    virtual atom<I> atom_is(unsigned int i, const nothing&) const
    {
        return { fvalues[i], 1.0 / fvalues.size() };
    }
    std::vector<I> fvalues;
    bool fsorted;
    virtual bool is_sorted() const { return fsorted; }
    virtual bool is_equiprobable(const nothing&) const { return true; }
};


/// \brief Alternative lists defined distribution
/// \tparam I type of the values
///
/// The atoms sets are based on the condition (of type <tt>unsigned int</tt>).
template <typename I>
class altldistribution: public fdistribution<I,unsigned int>
{
    std::vector<ptr<ldistribution<I>>> makedists(const std::vector<std::vector<I>>& values)
    {
        std::vector<ptr<ldistribution<I>>> ret(values.size());
        for(unsigned int i=0; i< values.size(); i++)
            ret[i].reset(new ldistribution<I>(values[i]));
        return ret;
    }
public:
/*    altldistribution(const std::vector<std::vector<atom<I>>>& atoms) :
        fdists(makedists(atoms))
    {
    }*/

private:
    virtual void atoms_are(std::vector<atom<I>>& a, const unsigned int& c) const
    {
        assert(c<fdists.size());
        fdists[c]->atoms(a,na);
    }

    virtual unsigned int natoms_is(const unsigned int& c) const
    {
        assert(c<fdists.size());
        return fdists[c]->natoms(na);
    }

    virtual atom<I> atom_is(unsigned int i,  unsigned int& c) const
    {
        assert(c<fdists.size());
         return fdists[c]->atom(i,na);
    }

    virtual bool is_equiprobable(const unsigned int& c) const
    { return fdists[c]->is_equiprobable(); }

    std::vector<ptr<ldistribution<I>>> fdists;

};


/// \brief List defined std::vector distribution
/// \tparam X type of the values' components
///
template <typename X>
class lvdistribution: public ldistribution<std::vector<X>>,
        virtual public vdistribution<X,nothing>
{
public:
    lvdistribution(const std::vector<atom<std::vector<X>>>& atoms) :
        ldistribution<std::vector<X>>(atoms)
    {
        assert(atoms.size());
        assert(atoms.size());
        for(unsigned int i=1; i<atoms.size(); i++)
            assert(atoms[i].x.size()==atoms[0].x.size());
    }

    lvdistribution(const std::vector<std::vector<X>>& values) :
        ldistribution<std::vector<X>>(values)
    {
        assert(values.size());
        assert(values[0].size());
        for(unsigned int i=1; i<values.size(); i++)
            assert(values[i].size()==values[0].size());
    }
private:
    virtual unsigned int dim_is() const
    {
        return (*this)(0).x.size();
    }
};


/// \brief Dirac distribution
/// \tparam I type of the value
///
template <typename I>
class diracdistribution: public ldistribution<I>
{
    std::vector<I> point(const I& a) const
    {
       std::vector<I> ret = {a};
       return ret;
    }
public:
    diracdistribution(const I& a)
        : ldistribution<I>(point(a)) {}
    I x() const
    {
        atom<I> a = (*this)(0);
        return a.x;
    }
};

/// \brief Product distribution
template <typename X>
lvdistribution<X> operator *(const fdistribution<X,nothing>& x,
                             const fdistribution<X,nothing>& y)
{
    assert(x.natoms());
    assert(y.natoms());
    std::vector<atom<std::vector<X>>> a(x.natoms()*y.natoms());

    unsigned int dst = 0;
    for(unsigned int i=0; i<x.natoms(); i++)
        for(unsigned int j=0; j<y.natoms(); j++)
        {
            atom<X> xa = x(i);
            atom<X> ya = y(j);
            std::vector<X> r = { xa.x, ya.x };
            a[dst++]= { r, xa.p*ya.p};
        }
    return lvdistribution<X>(a);
}

template <typename D, typename C>
class discretization : public fdistribution<double, C>
{
public:
  discretization(const D& d, unsigned int n) : fd(d), fn(n)
  {
      assert(fn);
  }
private:
  virtual unsigned int natoms_is(const C& c) const
  {
      return fn;
  }
  virtual atom<double> atom_is(unsigned int i, const C& c) const
  {
      probability p= static_cast<double> (2*i+1) / static_cast<double>(2*fn);
      return { fd.quantile(p,c), 1.0 / static_cast<double> (fn) };
  }
  virtual bool is_equiprobable(const C& c) const { return true; }
  D fd;
  unsigned int fn;
};

/// @} - discrete distributions

/// \addtogroup multidim Multidimensional distributions
/// \ingroup Distributions
/// @{

template <typename  D>
class meanvardistribution:
        virtual public mdistribution<std::vector<double>,nothing>,
        virtual public vdistribution<double,nothing>
{
    void check()
    {
        static_assert(std::is_base_of
          <mdistribution<typename D::I_t,nothing>,D>::value);
        static_assert(std::is_convertible
           <typename D::I_t,double>::value);
        assert(fm.size() == fsqV.size());
        for(unsigned int i=0; i<fsqV.size(); i++)
            assert(fsqV[i].size()==fm.size());
    }
public:
    meanvardistribution(const std::vector<double>& m,
                        const std::vector<std::vector<double>>& sqV
                        )
        : fm(m), fsqV(sqV)
    {
        check();
    }

    meanvardistribution(const std::vector<double>& m,
                        const std::vector<std::vector<double>>& sqV,
                        const D& d
                        )
        : fm(m), fsqV(sqV),fd(d)
    {
        check();
    }
    const D& d() const { return fd; }
    const std::vector<double>& m() const {return fm; }
    const std::vector<std::vector<double>>& sqV() const { return fsqV; }

protected:
    virtual unsigned int dim_is() const { return fm.size(); }
    virtual std::vector<double> do_draw(const nothing&) const
    {
        unsigned int d = this->dim();
        std::vector<double> r(fm);
        for(unsigned i=0; i<d; i++)
        {
            double u = fd.draw();
            for(unsigned int j=0; j<d; j++)
                r[j] += fsqV[j][i] * u;
        }
        return r;
    }
private:
    std::vector<double> fm;
    std::vector<std::vector<double>> fsqV;
    D fd;
};

template <typename  D>
class fmeanvardistribution:
        public meanvardistribution<D>,
        public fdistribution<std::vector<double>,nothing>
{
    void init(unsigned int dim)
    {
        static_assert(std::is_base_of
          <fdistribution<double,nothing>,D>::value);
        fn = this->d().natoms(na);
        fN = 1;
        for(unsigned int i=0; i<dim; i++)
            fN *= fn;
    }
public:
    fmeanvardistribution(const std::vector<double>& m,
                        const std::vector<std::vector<double>>& sqV
                        ) :
        meanvardistribution<D>(m,sqV)
    {
        init(m.size());
    }

    fmeanvardistribution(const std::vector<double>& m,
                        const std::vector<std::vector<double>>& sqV,
                        const D& d
                        ):
        meanvardistribution<D>(m,sqV,d)
    {
        init(m.size());
    }

private:
    virtual unsigned int dim_is() const { return meanvardistribution<D>::dim_is(); }
    virtual atom<std::vector<double>> atom_is(unsigned int a, const nothing&) const
    {
        unsigned int d=this->dim();
        std::vector<double> r(this->m());
        unsigned int ind = a;
//sys::log() << a ;
        for(unsigned int i=0; i<d; i++)
        {
            unsigned int k = ind % fn;
            double u=this->d()(k).x;
//sys::log() << "," << u;
            for(unsigned int j=0; j<d; j++)
                r[j] += this->sqV()[j][i] * u;
            ind /= fn;
        }
//for(unsigned int i=0; i<r.size(); i++)
//   sys::log() << "," << r[i];
//sys::log() << endl;
        return { r,1.0 / (double) fN};
    }
    virtual unsigned int natoms_is(const nothing& ) const { return fN;}
    virtual std::vector<double> do_draw(const nothing&) const
    {
        return meanvardistribution<D>::do_draw(na);
    }
    unsigned int fn;
    unsigned int fN;
};



/// @} - multidimensional distributions


/// \addtogroup iterativedists Iteratively defined distributions
/// \ingroup Distributions
/// @{


/// \brief Iterative joint distribution,
/// \tparam F the first distribution
/// \tparam S the second distribution
/// \tparam M a mapping om the
/// <tt>pair<condition of F, condition of S></tt> into
/// <tt>condition of \p F</tt>
///
/// Tne condtition of the distribution is the same type as that of \p F.
/// The distribution itself is defined iteratively by <em>F</em>|(condition of \e F) and
/// <em>S</em>|<em>M</em>(value of \e F,condition of <em>F</em>)

template <typename F, typename S, typename M>
class ijdistribution:
     public jdistribution<typename F::I_t, typename S::I_t, typename F::C_t>
{
public:
    using F_t = F;
    using S_t = S;
    using M_t = M;

    ijdistribution(const F& f, const S& s) : ff(f), fs(s)
    {
        static_assert(std::is_base_of<mapping<typename M::D_t,typename M::R_t>,M>::value);
        static_assert(std::is_same<typename M::D_t,
               std::pair<typename F::C_t,typename F::I_t>>::value);
        static_assert(std::is_convertible<typename M::R_t&,typename S::C_t&>::value);
    }
    const F& first() const { return ff; }
    const S& second() const { return fs; }
private:
    F ff;
    S fs;
};

/// Helper class, used as \p M in \p uijdistribution
template <typename D,typename R, typename M>
class uijmapping : public mapping<std::pair<nothing,D>,R>
{
public:
    virtual R operator() (const std::pair<nothing,D>& p) const
    {
        return M()(p.second);
    }
};

/// \brief Unconditional iterative distribution
/// \tparam F the first distribution
/// \tparam S the second distribution
/// \tparam M s mspping from the value of \p F into the condition of \f S
///
/// By default, \p M is the identity from value of \p F into condition of \p S.
/// If those types differ, then \p M should be explicitly specified.

template <typename F, typename S, typename M=idmapping<typename F::I_t>>
using uijdistribution = ijdistribution<F,S,uijmapping<typename F::I_t,typename S::C_t,M>>;

/// \brief MC iterative distribution
/// \tparam F the first distribution
/// \tparam S the second distribution
/// \tparam M a mapping from the condition of \p F into
/// <tt>(condition of F, condition of S)</tt>
///
/// See \ref ijdistribution
template <typename F, typename S, typename M>
class mijdistribution: public ijdistribution<F,S,M>,
                       virtual public mdistribution<std::pair<typename F::I_t,typename S::I_t>,typename F::C_t>
{
public:
    using I_t = typename ijdistribution<F,S,M>::I_t;
    using C_t = typename ijdistribution<F,S,M>::C_t;
    mijdistribution(const F& d, const S& e) :
      ijdistribution<F,S,M>(d,e)
    {
    }
protected:
    virtual I_t do_draw(const C_t& c) const
    {
        typename F::I_t i = this->first().draw(c);
        std::pair<C_t,typename F::I_t> p = { c, i };
        typename S::I_t j = this->second().draw(M()(p));
        return {i,j};
    }
};

/// \brief Monte Carlo iterative distribution with finite support
/// \tparam F the first distribution
/// \tparam S the second distribution
/// \tparam M a mapping from the condition of \p F into
/// <tt>(condition of F, condition of S)</tt>
///
/// See \ref ijdistribution
template <typename F, typename S, typename M>
class fijdistribution: public mijdistribution<F,S,M>,
        public fdistribution<std::pair<typename F::I_t,typename S::I_t>,
                                          typename F::C_t>
{
public:
    using I_t = typename ijdistribution<F,S,M>::I_t;
    using C_t = typename ijdistribution<F,S,M>::C_t;
    fijdistribution(const F& d, const S& e) :
      mijdistribution<F,S,M>(d,e)
    {
        static_assert(std::is_base_of<
            fdistribution<typename F::I_t, typename F::C_t>,F>::value);
        static_assert(std::is_base_of<
            fdistribution<typename S::I_t, typename S::C_t>,S>::value);
    }
private:
    virtual void atoms_are(std::vector<atom<I_t>>& a, const C_t& c) const
    {
        a.resize(0);
        std::vector<atom<typename F::I_t>> fa;
        this->first().atoms(fa,c);

        for(unsigned int i=0; i<fa.size(); i++)
        {
            std::vector<atom<typename S::I_t>> sa;
            this->second().atoms(sa,M()(std::pair<C_t,typename F::I_t>(c,fa[i].x)));
            for(unsigned int j=0; j<sa.size(); j++)
                a.push_back({ { fa[i].x, sa[j].x}, fa[i].p * sa[j].p});
        }
    }
    virtual I_t do_draw(const C_t& c) const
    {
        return mijdistribution<F,S,M>::do_draw(c);
    }

};

/// \brief Unconditional MC iterative distribution
/// \tparam F the first distribution
/// \tparam S the second distribution
/// \tparam M s mspping from the value of \p F into the condition of \f S
///
/// See \ref uijdistribution

template <typename F, typename S, typename M>
using umijdistribution = mijdistribution<F,S,uijmapping<typename F::I_t,typename S::C_t,M>>;

/// Scenario (history) of a stochastic process with values in \p X
template <typename X=double>
using scenario = std::vector<X>;


/// \ingroup Processes
/// \brief Base class for scenario restrictions
/// \tparam X type of the scenario components
/// \tparam R type into which the scenario is transformed
///
/// Scenario restrictions are mappings transforming scenarios
/// into some other space, which server either as conditions
/// of stochastic processes or as random parameters of decision problems.
///
template <typename X,typename R>
class zeta: public mapping<scenario<X>,R>
{
public:
    using R_t=R;
    using X_t=X;
};

/// \brief Iteratively defined std::vector distribution
/// \tparam D conditional distribution of the std::vector's component
/// \tparam Z mapping from a history of the distribution into the condition of \p D
///
/// The distribution is defned iteratively, the \e i-th component being
/// conditioned by <tt>{Z}(s)</tt>$ where \p s is the scenario
/// containing previous values of the std::vector.
/// Template parameter \p Z has to be a descendant of \ref zeta<\p X, \p C>
/// where \p X is the type of \p D's values and \p C is the type
/// of its condition.
///

template <typename E, typename D, typename Z>
class ivdistribution:
    virtual public vdistribution<typename D::I_t,nothing>
{
    static void constexpr check()
    {
        static_assert(std::is_base_of<zeta<typename Z::X_t, typename Z::R_t>,Z>::value);
        static_assert(std::is_base_of<
                      distribution<typename D::I_t,typename D::C_t>,D>::value);
        static_assert(std::is_base_of<
                   distribution<typename E::I_t,nothing>,E>::value);
        static_assert(std::is_same<typename Z::R_t,typename D::C_t>::value);
        static_assert(std::is_same<typename Z::X_t,typename D::I_t>::value);
        static_assert(std::is_same<typename E::I_t,typename D::I_t>::value);
    }
public:
    using X_t = typename D::I_t;
    using E_t = E;
    using Z_t = Z;
    using D_t = D;
    ivdistribution(const E& e, const std::vector<D>& d)
       : fd(d), fe(e)
    {
       check();
    }
    ivdistribution(const E& e, const D& d, unsigned int dim)
       : fe(e), fd(dim-1,d)
    {
        check();
        assert(dim);
    }

    /// usable only if \p E is \ref diracdistribution
/*    ivdistribution(const X_t& x0, const D& d, unsigned int dim)
       : fe(diracdistribution<X_t>(x0)), fd(dim-1,d)
    {
        check();
        assert(dim);
    }
*/
    typename Z::R_t z(const scenario<X_t>& s) const
    {
        return Z_t()(s);
    }

    /// not sure whether I use it
    std::vector<X_t> draw() const
    {
        std::vector<X_t> ret;
        X_t x = fe.draw(na);
        ret.push_back(x);
        typename Z::R_t c = z(ret);
        for(unsigned int i=0; i<this->dim()-1; i++)
        {
            x = this->fd[i].draw(c);
            ret.push_back(x);
            c = z(ret);
        }
        return ret;
    }

    /// may be used only if \p D is a descendant of \ref fdistribution
    ///
    /// TBD maybe adding this functionality to fdistribution wwould be more clean...
    void atoms(const std::vector<X_t>& s, std::vector<atom<X_t>>& a) const
    {
        if(!s.size())
            fe.atoms(a);
        else
        {
            assert(s.size() < this->dim());
            fd[s.size()-1].atoms(a,Z()(s));
        }
    }
    const D& d(unsigned int i) const
    {
        assert(i>0);
        assert(i-1 < fd.size());
        return fd[i-1];
    }
    const E& e() const
    {
        return fe;
    }
private:
    E fe;
    std::vector<D> fd;
protected:
    virtual unsigned int dim_is() const
    {
        return fd.size()+1;
    }
};

/// \brief MC version of \ref ivdistribution

template <typename E, typename D, typename Z, typename M>
class mivdistribution:
    public ivdistribution<E,D,Z>
{
    static void constexpr check()
    {
        static_assert( std::is_base_of<distribution<typename M::I_t,nothing>, M>::value);
        static_assert( std::is_same<typename M::I_t,typename D::I_t>::value);
    }
public:
    using E_t = E;
    using M_t = M;
    using D_t = D;
    using Z_t = Z;
    mivdistribution(const E& e, const D& d, unsigned int dim)
       : ivdistribution<E,D,Z>(e,d,dim)

    {
        check();
        assert(dim);
    }

    mivdistribution(const typename E::I_t& x0, const D& d, unsigned int dim)
        : ivdistribution<E,D,Z>(x0,d,dim)
    {
        check();
        assert(dim);
    }
    // marginal distribution
    M md(unsigned int i) const
    {
        assert(i>0);
        assert(i<this->dim());
        return md_is(i);
    }
private:
    virtual M md_is(unsigned int i) const = 0;
};


/*
template <typename D, typename E, typename Z>
class dipdistribution:
        public ijdistribution<D,E,Z>,
        public ddistribution<pair<typename D::I_t,typename E::I_t>,
                                               typename D::C_t>
{
public:
    using I_t = typename iterativedistribution<D,E>::I_t;
    using C_t = typename iterativedistribution<D,E>::C_t;
    diterativedistribution(const D& d, const E& e) :
        iterativedistribution<D,E>(d,e)
    {
    }
private:
    virtual void atoms_are(std::vector<atom<I_t>>& a, const C_t& c) const
    {
        a.clear();
        std::vector<atom<typename D::I_t>> da;
        this->d().atoms(da,c);
        unsigned int nd=da.size();
        for(unsigned int i=0; i<nd; i++)
        {
            atom<typename D::I_t> x = da[i];

            std::vector<atom<typename E::I_t>> ea;
            this->e().atoms(ea,x.x);
            unsigned int ne=ea.size();

            for(unsigned int j=0; j<ne; j++)
            {
                atom<typename E::I_t> y = ea[j];
                a.push_back({{x.x,y.x},x.p*y.p});
            }
        }
    }

    virtual atom<I_t> atom_is(unsigned int i, const C_t& c) const
    {
        std::vector<atom<I_t>> a;
        atoms_are(a,c);
        assert(i<a.size());
        return a[i];
    }

    virtual unsigned int natoms_is(const C_t& c) const
    {
        std::vector<atom<I_t>> a;
        atoms_are(a,c);
        return a.size();
    }

    virtual I_t do_draw(const C_t& c) const
    {
        I_t r;
        r.i = this->d().draw(c);
        r.j = this->e().draw(r.i);
        return r;
    }
};

*/

/// @} - iterative distributions

/// @} - distribution


/// \addtogroup Samples
/// @{


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

/// @} - Samples


/// \addtogroup Risk Measures
/// @{

template <typename Distribution>
class riskmeasure : public object
{
public:
    virtual double operator () (const Distribution&,
                                const typename Distribution::C_t&
                                ) const = 0;
    double operator () (const Distribution& d) const
    {
        return operator()(d,nothing());
    }
};

template <typename Distribution>
class estimableriskmeasure : public riskmeasure<Distribution>
{
public:
    virtual realestimate operator ()
       (const empiricaldistribution& d) const = 0;
};


using empriricalriskmeasure = riskmeasure<empiricaldistribution>;

template <typename Distribution, bool negative = false>
class MeanCVaR : public estimableriskmeasure<Distribution>
{
public:
    MeanCVaR(probability alpha, double lambda) : falpha(alpha),
        flambda(lambda)
    {         assert(falpha < 1);  assert(flambda <=1 && flambda >=0); }

    /// tbd can be done for sorted ddistribution
    virtual double operator () (
         const Distribution& d,
         const typename Distribution::C_t& c) const
    {
        // tbd static assert here
        std::vector<atom<typename Distribution::I_t>> a;
        if(!d.sorted())
            throw("The distribution has to be sorted.");
        d.atoms(a,c);
        double p = 1-falpha;
        double s = 0;
        double m = 0;
        double mean = 0;
        assert(a.size());
        for(int i = 0; i<a.size(); i++)
        {
            int index;
            if constexpr(negative)
                 index = i;
            else
                 index = a.size()-1-i;
            double delta = std::min(a[index].p, p-s);
            mean += a[index].x * a[index].p;
            if(delta > 0)
            {
                m+= a[index].x * delta;
                s += delta;
            }
        }
        double cvar = m / p;
        return (1-flambda)*mean + flambda*cvar;
    }

    virtual realestimate operator () (const empiricaldistribution& d) const
    {
        statcounter sc;
        unsigned n = d.natoms();
        assert(n);


        double p = 1-falpha;
        double s = 0;
        for(int i = 0; i< d.natoms() ;i++)
        {
            int index;
            if constexpr(negative)
                index = i;
            else
                index = d.natoms()-1-i;

            double x = d.value(index);
            double delta = std::min(1.0 / n, p-s);
            double inc = x * (1-flambda);
            if(delta >= 0)
            {
               inc += d.value(index) / p * delta * n * flambda;
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


template <typename Distribution, bool negative = false>
class CVaR : public MeanCVaR<Distribution, negative>
{
public:
    CVaR(probability alpha) : MeanCVaR<Distribution,negative>(alpha,1)
     {}
};


/// @} - Risk Measures



} // namespace


#endif // RANDOM_H

