#ifndef DP_HPP
#define DP_HPP

#include "orpp/optimization.hpp"
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

template <typename Element>
class indexedspace : public space<Element>
{
public:
    Element operator[] (index& i) const
    {
        Element res;
        if(!find(res,i))
            throw exception("Index out of bounds");
        return res;
    }
private:
    virtual bool find(Element&, index) const = 0;
};


template <typename Element, typename Condition = nothing>
class iteratedspace : public indexedspace<Element>
{
public:
    bool first(Element& e, const Condition& c, index& i)
    {
        if(getfirst(e,c))
        {
            i=0;
            return true;
        }
        else
            return false;
    }
    bool next(Element& e, const Condition& c, index& i)
    {
        if(getnext(e,c))
        {
            i++;
            return true;
        }
        else
            return false;
    }
    bool firstfeasible(Element& e, const Condition& c) const
    {
        if(getfirst(e,c))
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
            if(!getnext(e,c))
                return false;
            if(feasible(e,c))
                return true;
        }
        return false;
    }
    virtual bool isfeasible(const Element&, const Condition&) const { return true; }
private:
    virtual bool getfirst(Element&, const Condition& c) const = 0;
    virtual bool getnext(Element& e,const Condition& c) const = 0;
    virtual bool find(Element&, index) const { return false; }
};

template <typename Condition = nothing>
class integeriteratedspace : public iteratedspace<int, Condition>
{
public:
    integeriteratedspace(int amin, int amax) : fmin(amin), fmax(amax)
    {}
    virtual bool getfirst(int& e, const Condition&) const
    {
        e=fmin;
        return fmin <= fmax;
    }
    virtual bool getnext( int& e,const Condition&) const
    {
        if(e >= fmax)
            return false;
        e++;
        return true;
    }
    /// isfeasible colides with find

private:
    int fmin;
    int fmax;
    virtual bool find(int& e, index i) const
    {
        if(i <= fmax-fmin)
        {
            e = i-fmin;
            return true;
        }
        else
            return false;
    }
};



/*template <typename Element>
class statespace: public object
{
public:
    using Element_t = Element;
};

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
         fcrit(crit), fstate(state),
         faction(action), ftransition(transition), freward(reward), fgamma(gamma)
    {}
    double gamma() const { return fgamma; }
protected:
    Criterion fcrit;
    Statespace fstate;
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

    statcounter evaluate(index s0index, const finitepolicy& p,
                    unsigned timehorizon, double accuracy, unsigned miniters =100,
                         unsigned maxiters = 100000)
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
                c.s = this->fstate[sindex];
                c.a = p[sindex];
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
