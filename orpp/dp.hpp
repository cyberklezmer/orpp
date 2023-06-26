#ifndef DP_HPP
#define DP_HPP

//#include "orpp/optimization.hpp"
#include "orpp/simulations.hpp"
#include "orpp/random.hpp"

//#include <fstream>

namespace orpp
{

template <typename Element>
class space : public object
{
public:
    using Element_t = Element;
};



template <typename Element, typename Condition = nothing>
class iteratedspace : public space<Element>
{
public:
    bool first(Element& e) const
    {
        return getfirst(e);
    }
    bool next(Element& e) const
    {
        return getnext(e);
    }
    bool firstfeasible(Element& e, const Condition& c) const
    {
        if(first(e))
        {
           if(isfeasible(e,c))
               return true;
           else
               return nextfeasible(e,c);
        }
        else
            return false;
    }
    bool nextfeasible(Element& e, const Condition& c) const
    {
        for(;;)
        {
            if(!getnext(e))
                return false;
            if(feasible(e,c))
                return true;
        }
        return false;
    }
    virtual bool isfeasible(const Element&, const Condition&) const { return true; }
    bool ith(Element& e, index i) const { return getith(e,i);}
    Element operator[] (index& i) const
    {
        Element res;
        if(!getith(res,i))
            throw exception("Index out of bounds");
        return res;
    }
    index num() const { return getnum(); }
private:
    virtual bool getfirst(Element&) const = 0;
    virtual bool getnext(Element& e) const = 0;
    virtual bool getith(Element& e, index i) const
    {
        if(!first(e))
            return false;
        index j=0;
        for(;j<i;j++)
        {
            if(!next(e))
                return false;
        }
        return true;
    }
    virtual index getnum() const
    {
        Element e;
        if(!first(e))
            return 0;
        for(index i=1;;i++)
        {
            if(!next(e))
                return i;
        }
        throw "souhld not get here";
    }
};

/*template <typename Element>
class indexedspace : public space<Element>
{
public:
private:
    virtual bool find(Element&, index) const = 0;
};*/


template <typename Condition = nothing>
class integeriteratedspace : public iteratedspace<int, Condition>
{
public:
    integeriteratedspace(int amin, int amax) : fmin(amin), fmax(amax)
    {}

private:
    virtual bool getfirst(int& e) const
    {
        e=fmin;
        return fmin <= fmax;
    }
    virtual bool getnext( int& e) const
    {
        if(e >= fmax)
            return false;
        e++;
        return true;
    }
    virtual bool getith(int& e, index i) const
    {
        e = i+fmin;
        return e<=fmax;
    }
    virtual index getnum() const
    {
        return fmax - fmin + 1;
    }
    int fmin;
    int fmax;
};



/*
template <typename Element>
class actionspace: public feasibilityset<Element>
{
public:
    using Element_t = Element;
};



template <typename E>
class iteratedstatespace: public iteratedspace<E> //, statespace<int>
{

};

class integerstatespace: public integeriteratedspace //, statespace<int>
{
public:
    integerstatespace(int amin, int amax) : integeriteratedspace(amin,amax)
    {}

};

template <typename Element, typename Condition>
class iteratedactionspace: public iteratedspace<Element,Condition> //, actionspace<int>
{

};

*/

template <typename SElement, typename AElement>
struct dpcondition
{
    AElement a;
    SElement s;
};

template <typename S, typename A>
using reward=function<dpcondition<typename S::Element_t,typename A::Element_t>>;



/// \brief Base abstract class - dp problem
/// \tparam O risk criterion
/// \tparam R reward
/// \tparam A action space
/// \tparam S state space
/// \tparam T transition operator
template<typename Criterion,
          typename Statespace,
          typename Actionspace,
         typename Policy,
         typename Transition,
         typename Reward = reward<Statespace,Actionspace>
          >
class dpproblem : public object
{
public:
    using Statespace_t = Statespace;
    using Actionspace_t = Actionspace;
public:
    dpproblem(const Criterion& crit,
              const Statespace& state,
              const Actionspace& action,
              const Transition& transition,
               const Reward& reward,
              double gamma ) :
         fcrit(crit), fstatespace(state),
         faction(action), ftransition(transition), freward(reward), fgamma(gamma)
    {}
    double gamma() const { return fgamma; }
protected:
    Criterion fcrit;
    Statespace fstatespace;
    Actionspace faction;
    Transition ftransition;
    Reward freward;
    double fgamma;
};

template <typename Statespace, typename Actionspace>
using finitetransition
 = ddistribution<typename Statespace::Element_t, dpcondition<Statespace,Actionspace>>;

using finitepolicy = std::vector<index>;

template<typename Criterion,
         typename Statespace,
         typename Actionspace,
         typename Transition = finitetransition<Statespace,Actionspace>,
         typename Reward = reward<Statespace,Actionspace>>
class finitedpproblem : public dpproblem<Criterion,Statespace,Actionspace,
                                 finitepolicy, Transition,Reward>
{
public:
//    using CritDistribution = fdistribution<double,nothing>;

    finitedpproblem(const Criterion& crit,
              const Statespace& state,
              const Actionspace& action,
              const Transition& transition,
              const Reward& reward,
              double gamma ) :
         dpproblem<Criterion, Statespace,Actionspace, finitepolicy,
                          Transition,Reward>
           (crit, state,action, transition, reward, gamma)
    {
    }

    statcounter evaluateraw(index s0index,
                         const std::vector<finitepolicy>& p,
                         unsigned timehorizon,
                         double accuracy,
                         unsigned miniters =100,
                         unsigned maxiters = 100000) const
    {
        statcounter sc;
        for(unsigned j=0; j<maxiters; j++)
        {
            double discount = 1;
            double sum = 0;
            index sindex = s0index;
            for(unsigned i=0; i<=timehorizon; i++, discount *= this->fgamma)
            {
                dpcondition<int,int> c;
                c.s = this->fstatespace[sindex];
                c.a = i<p.size() ? p[i][sindex] : p[p.size()-1][sindex];
                sum += discount * this->freward(c);
                sindex = this->ftransition.draw(c);
                discount *= this->fgamma;

/* test of distribution                if(i==0 && j==0)
                {
                    std::ofstream o("dist.csv");
                    std::cout << "c.s=" << c.s << ", c.a=" << c.a << std::endl;
                    for(unsigned k=0; k<50000; k++)
                        o << k << "," << this->ftransition.draw(c) << std::endl;
                    throw;
                } */

//std::cout << j << " - " << i << ", s=" << sindex << ", sum="  << sum << std::endl;
            }
            sc.add(sum);
            if(j>miniters && sc.averagestdev() < accuracy)
                break;
        }
        return sc;
    }

    using vapprox = std::vector<double>;



    valuewitherror<vapprox> evaluate(const vapprox& initialV,
                    double initialerr,
                    finitepolicy& p,
                         double accuracy,
                         unsigned maxiters = 100000) const
    {
        vapprox V = initialV;

        assert(V.size() == this->fstatespace.num());

        double error = initialerr;

        for(unsigned j=0; j<maxiters; j++)
        {
std::cout << error << std::endl;
            vapprox newV(V.size());
            typename Statespace::Element_t e;
            this->fstatespace.first(e);
            for(index i = 0; i < V.size(); i++, this->fstatespace.next(e))
            {
                dpcondition<int,int> c;
                c.s = e;
                assert(i < p.size());
                c.a = p[i];
                double r =  this->freward(c);
                double v=0;
                std::vector<atom<double>> atoms;
                for(unsigned k=0; k<this->ftransition.natoms(c); k++)
                {
                    atom<typename Transition::I_t> a =
                            this->ftransition(k,c);
                    assert(a.x == k);
                    atom<double> newa = {a.p,this->fgamma*V[k]+r};

                    auto insertionPos = std::lower_bound(atoms.begin(), atoms.end(), newa,
                                                         atom<double>::comparator);
                    atoms.insert(insertionPos, newa);
                }
                ldistribution<double> d(atoms, false, false);
                newV[i]=this->fcrit(d,nothing());
            }
            V = newV;
            error *= this->fgamma;
            if(error < accuracy)
                break;
        }
        return { V, error };
    }


};



template <typename P, typename D>
class valueitersolution;

template <typename P, typename D, typename SolutionInfo = nothing>
class dpsolution
{
    friend class valueitersolution<P,D>;

    dpsolution(const std::pair<ptr<typename P::Policy>,SolutionInfo>& s) :
            fa(s.first), fo(s.second)
    {}
private:
    SolutionInfo fo;
    ptr<typename P::Policy> fa;
};

using VIInfo = nothing;

template <typename P, typename D>
class valueitersolution: public dpsolution<P,D,VIInfo>
{
public:
    struct viparams
    {
        double initialV = 0;
//        typename P::Statespace_t::Element initialState;
        double epsilon = 0.1;
    };
    valueitersolution(const P& p, viparams vp) :
        dpsolution<P,D>(solve(p,vp))
    {}

    using Result = std::pair<ptr<typename P::Policy>,VIInfo>;
private:

    static Result solve(const P& p, viparams vp)
    {
       std::vector<double> V;
       typename P::Statespace_t::Element s;
       index sindex;
       bool lasts = p.fstate.first(s,sindex);
       for(;lasts;)
       {
           V.push_back(vp.initialV);
           lasts = p.fstate.next(s,sindex);
       }
       auto newV = V;
       typename P::Policy policy;
       for(;;)
       {
           bool lasts = p.fstate.first(s,sindex);
           policy.clear();
           for(;lasts; lasts = p.fstate.next(s,sindex))
           {
               typename P::Actionspace_t::Element a, besta;
               bool lasta = p.faction.firstfeasible(a,s);
               double mincrit = infinity<double>;
               for(;lasta;lasta = p.faction.nextfeasible(a,s))
               {
                   double r = p.freward(s,a);

                   typename P::Statespace_t::Element news;

                   std::vector<atom<double>> dstatoms;

                   for(index i = 0; i<V.size(); i++)
                   {
                       atom<index> srca = ftransition(i,{a,s});
                       assert(i<V.size());
                       dstatoms.push_back({-(r+p.fgamma*V[i]),srca.p});
                   }
                   ldistribution<double> dstd(dstatoms);
                   double crit = p.fcrit(dstd);
                   if(crit < mincrit)
                   {
                       mincrit = crit;
                       a = besta;
                   }
               }
               newV[sindex] = -mincrit;
               policy.push_back(besta);
           }
           double m = 0;
           for(unsigned i=0; i<V.size(); i++)
           {
               double d = fabs(V[i]-newV[i]);
               if(d>m)
                   m=d;
           }
           if(m<vp.epsilon)
               break;
           V = newV;
       }
       ptr<Result> res(new ptr<Result>({policy, nothing()}));
       return res;
    }
};



} // namespace

#endif // DP_HPP
