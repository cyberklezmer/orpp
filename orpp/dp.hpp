#ifndef DP_HPP
#define DP_HPP

//#include "orpp/optimization.hpp"
#include "orpp/random.hpp"
#include <thread>
//#include <fstream>

namespace orpp
{

template <typename Element>
class space : public object
{
public:
    using Element_t = Element;
};


/// \brief Base abstract class - dp problem
/// \tparam O risk criterion
/// \tparam R reward
/// \tparam A action space
/// \tparam S state space
/// \tparam T transition operator
template<typename Criterion,
          typename Statespace,
          typename ConstrainedActionSpace,
          typename Transition,
          typename Reward
          >
class dpproblem : public object
{
public:
    using Statespace_t = Statespace;
    using ConstrainedActionSpace_t = ConstrainedActionSpace;
public:
    dpproblem(const Criterion& crit,
              const Statespace& state,
              const ConstrainedActionSpace& constraint,
              const Transition& transition,
               const Reward& reward,
              double gamma ) :
         fcrit(crit), fstatespace(state),
         fconstraint(constraint),
         ftransition(transition), freward(reward), fgamma(gamma)
    {}
    double gamma() const { return fgamma; }
    const ConstrainedActionSpace& constraint() const { return fconstraint; }
    const Statespace& statespace() const { return fstatespace; }
protected:
    Criterion fcrit;
    Statespace fstatespace;
    ConstrainedActionSpace fconstraint;
    Transition ftransition;
    Reward freward;
    double fgamma;
};

template <typename Element>
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



template <typename SElement, typename AElement>
class constrainedspace : virtual public iteratedspace<AElement>
{
public:
    bool firstfeasible(AElement& e, const SElement& c) const
    {
        if(this->first(e))
        {
           if(isfeasible(e,c))
               return true;
           else
               return nextfeasible(e,c);
        }
        else
            return false;
    }
    bool nextfeasible(AElement& e, const SElement& c) const
    {
        for(;;)
        {
            if(!this->next(e))
                return false;
            if(isfeasible(e,c))
                return true;
        }
        return false;
    }
    bool feasible(const AElement& a, const SElement& s) const
    {
        return isfeasible(a,s);
    }
private:
    virtual bool isfeasible(const AElement&, const SElement& ) const
      { return true; }
};


template <typename SElement, typename AElement>
struct dpcondition
{
    AElement a;
    SElement s;
};



template <typename StateSpace, typename ActionSpace>
using dpreward=function<dpcondition<typename StateSpace::Element_t,
                                  typename ActionSpace::Element_t>>;


//template <typename StateSpace, typename ActionSpace>
//using =mapping<typename StateSpace::Element_t,
//                     typename ActionSpace::Element_t>;



class integerspace : virtual public iteratedspace<unsigned int>
{
public:
    integerspace(unsigned int amin, unsigned int amax) : fmin(amin), fmax(amax)
    {}

private:
    virtual bool getfirst(unsigned int& e) const
    {
        e=fmin;
        return fmin <= fmax;
    }
    virtual bool getnext(unsigned int& e) const
    {
        if(e >= fmax)
            return false;
        e++;
        return true;
    }
    virtual bool getith(unsigned int& e, index i) const
    {
        e = i+fmin;
        return e<=fmax;
    }
    virtual index getnum() const
    {
        return fmax - fmin + 1;
    }
    unsigned int fmin;
    unsigned int fmax;
};



template <typename Statespace, typename Actionspace>
using finitetransition
 = fdistribution<typename Statespace::Element_t,
   dpcondition<typename Statespace::Element_t,typename Actionspace::Element_t>>;

/*using finitedpstatespace = integerspace;
using finitedpactionspace =  constrainedspace<unsigned int,unsigned int>;
using finitedptransition = finitetransition<finitedpactionspace, finitedpactionspace>;
using finitedpreward = dpreward<finitedpactionspace, finitedpactionspace>;
*/


class finitepolicy : public std::vector<index>, public mapping<index,index>
{
public:
    index operator () (const index& i) const
    {
        assert(i < this->size());
        return (*this)[i];
    }
    bool operator==(const finitepolicy& b) const
    {
        assert(b.size() == this->size());
        for(unsigned i=0; i<b.size(); i++)
        {
            if((*this)[i] != b[i])
                return false;
        }
        return true;
    }
    finitepolicy(const std::vector<index> v) : std::vector<index>(v) {}

    finitepolicy(unsigned size) : std::vector<index>(size) {}
//    template <typename Statespace>
//    finitepolicy(const Statespace& s, index initial = 0) : std::vector<index>(s.num(),initial) {}
    template <typename Problem>
    finitepolicy(const Problem& p, index initial = 0) : std::vector<index>(p.statespace().num(),initial) {}
    finitepolicy(const finitepolicy& p): std::vector<index>(p) {}
};

class finitevaluefunction : public std::vector<double>, public mapping<index,double>
{
public:
    double operator () (const index& i) const
    {
        assert(i < this->size());
        return (*this)[i];
    }
    template <typename Problem>
    finitevaluefunction(const Problem& p, double initial = 0)
        : std::vector<double>(p.statespace().num(),initial) {}
    finitevaluefunction(unsigned num, index initial = 0 )
        : std::vector<double>(num,initial) {}
};


inline double dist(const finitevaluefunction& a, const finitevaluefunction& b)
{
    double d=0;
    assert(a.size() == b.size());
    for(unsigned i=0; i<a.size(); i++)
        d = std::max(d, fabs(a[i]-b[i]) );
    return d;
}



template<typename Criterion,
         typename Statespace,
         typename ConstrainedActionSpace,
         typename Transition,
         typename Reward>
class finitedpproblem : public dpproblem<Criterion,Statespace,ConstrainedActionSpace,
                                 Transition,Reward>
{
public:
//    using CritDistribution = fdistribution<double,nothing>;

    struct computationparams
    {
        computationparams() : fmaxevaliterations(1000000),
            fthreadstouse(0), fthreadbatch(10) {}
        unsigned fmaxevaliterations;
        unsigned fthreadstouse;
        unsigned fthreadbatch;
    };

    finitedpproblem(const Criterion& crit,
              const Statespace& state,
              const ConstrainedActionSpace& constraint,
              const Transition& transition,
              const Reward& reward,
              double gamma,
              double maxreward) :
         dpproblem<Criterion, Statespace,ConstrainedActionSpace,
                          Transition,Reward>
           (crit, state, constraint, transition, reward, gamma),
         fmaxreward(maxreward)
    {
        static_assert(std::is_base_of<iteratedspace<typename Statespace::Element_t>,Statespace>::value);
        static_assert(std::is_base_of<constrainedspace<
                          typename Statespace::Element_t,
                          typename ConstrainedActionSpace::Element_t>,ConstrainedActionSpace>::value);
        static_assert(std::is_base_of<finitetransition<Statespace,ConstrainedActionSpace>,Transition>::value);
        static_assert(std::is_base_of<dpreward<Statespace,ConstrainedActionSpace>,Reward>::value);
    }


    unsigned requiredhorizon(double accuracy) const
    {
        return (log(accuracy * (1-this->gamma())) - log(maxreward())) / log( this->gamma());
    }

    valuewitherror<double> evaluatecrit(index s0ind,
                                         const std::vector<finitepolicy>& ps,
                                         double accuracy,
                                         const computationparams& params) const
    {
        statcounter sc = this->evaluateraw(s0ind, ps, accuracy, params);

        return this->fcrit(sc.dist);
    }
private:
     void computepath(orpp::index s0index,
                     const std::vector<finitepolicy>& p,
                     unsigned timehorizon, std::vector<double>* results) const
    {
        std::vector<double>& rs = *results;
        for(unsigned n=0; n<rs.size(); n++)
        {
            double discount = 1;
            double sum = 0;
            index sindex = s0index;
            for(unsigned i=0; i<=timehorizon; i++)
            {
                dpcondition<unsigned int,unsigned int> c;
                c.s = this->fstatespace[sindex];
                c.a = i<p.size() ? p[i][sindex] : p[p.size()-1][sindex];
                sum += discount * this->freward(c);
                sindex = this->ftransition.draw(c);
                discount *= this->fgamma;
            }
           rs[n] = sum;
        }
    }
public:
    statcounter evaluateraw(index s0index,
                         const std::vector<finitepolicy>& p,
                         double accuracy,
                         const computationparams& params) const
    {
        unsigned timehorizon = requiredhorizon(accuracy / 2);
        statcounter sc;
        for(unsigned j=0;j<params.fmaxevaliterations;)
        {
            if(params.fthreadstouse==0)
            {
                std::vector<double> sum(1);
                computepath(s0index,p, timehorizon, &sum);
                sc.add(sum[0]);
                j++;
            }
            else
            {
                std::vector<std::thread> ts;
                std::vector<std::vector<double>> rs(params.fthreadstouse,std::vector<double>(params.fthreadbatch));
                for(unsigned k=0; k<params.fthreadstouse; k++)
                    ts.push_back(std::thread(&finitedpproblem<Criterion, Statespace,ConstrainedActionSpace,
                                 Transition,Reward>::computepath,
                                 this, s0index,p, timehorizon, &rs[k]));
                for(unsigned k=0; k<params.fthreadstouse; k++)
                    ts[k].join();
                for(unsigned k=0; k<params.fthreadstouse; k++)
                    for(unsigned n=0; n<params.fthreadbatch; n++)
                       sc.add(rs[k][n]);
                j+=params.fthreadstouse * params.fthreadbatch;
             }

            if(j>10 && sc.averagestdev() < accuracy / 4)
                break;
        }
        return sc;
    }


    double maxreward() const { return fmaxreward; }
protected:
    double fmaxreward;

};

template<typename NestingCriterion>
class nestedcriterion {
public:
    using NestingCriterion_t = NestingCriterion;
    const NestingCriterion& nesting() const { return fnestingcrit; }
    NestingCriterion& nesting() { return fnestingcrit; }
    nestedcriterion(const NestingCriterion& nestingcrit) :
        fnestingcrit(nestingcrit) {}
private:
    NestingCriterion fnestingcrit;
};

template<typename Criterion,
         typename Statespace,
         typename ConstrainedActionSpace,
         typename Transition,
         typename Reward>
class finitehomodpproblem : public finitedpproblem<Criterion,Statespace,ConstrainedActionSpace,
                                 Transition,Reward>
{
public:
    finitehomodpproblem(const Criterion& crit,
              const Statespace& state,
              const ConstrainedActionSpace& constraint,
              const Transition& transition,
              const Reward& reward,
              double gamma,
              double maxreward) :finitedpproblem<Criterion,Statespace,ConstrainedActionSpace,
                                 Transition,Reward>(crit, state,constraint,
                                                 transition, reward, gamma,maxreward)
    {
        static_assert(std::is_base_of<nestedcriterion<typename Criterion::NestingCriterion_t>,Criterion>::value);
    }

    struct computationparams:
          public finitedpproblem<Criterion,Statespace,ConstrainedActionSpace,
            Transition,Reward>::computationparams
    {
        computationparams() : fvaluemaxiterations(100) {}
        unsigned fvaluemaxiterations;
    };

    struct viresult
    {
        viresult(const finitedpproblem<Criterion,Statespace,ConstrainedActionSpace,
                 Transition,Reward>& pr) : p(pr),v(pr), e(0) {}
        finitepolicy p;
        finitevaluefunction v;
        double e;
    };

    template <bool optimize>
    void iterate(viresult& params,
                            double accuracy,
                             const computationparams& evalparams) const
    {
        finitevaluefunction V = params.v;

        assert(V.size() == this->fstatespace.num());

        double error = infinity<double>;

        for(unsigned j=0; j<evalparams.fvaluemaxiterations; j++)
        {
    /* std::cout << error;
    for(unsigned k=0; k<V.size(); k++)
    std::cout << "," << V[k];
    std::cout << std::endl;*/
            finitevaluefunction newV(*this);
            typename Statespace::Element_t e;
            this->fstatespace.first(e);
            for(index i = 0; i < V.size(); i++, this->fstatespace.next(e))
            {
                dpcondition<unsigned int,unsigned int> c;
                typename ConstrainedActionSpace::Element_t a;
                double bestv = minfinity<double>;
                if constexpr(optimize)
                {
                   if(!this->fconstraint.firstfeasible(a,e))
                       throw exception("At least one feasible value missing in value iteration");
                }
                else
                {
                    a = params.p[i];
#ifndef NDEBUG
                    if(!this->fconstraint.feasible(a,e))
                        throw exception("Cannot evaluate an infeasible policy.");
#endif
                }
                for(index k=0; ; k++)
                {
                    c.s = e;
                    // assert(i < p.size());
                    c.a = a; // p[i];
                    double r =  this->freward(c);
                    double v=0;
                    std::vector<atom<double>> atoms;
                    for(unsigned k=0; k<this->ftransition.natoms(c); k++)
                    {
                        atom<typename Transition::I_t> a =
                              this->ftransition(k,c);
                        assert(a.x == k);
                        atom<double> newa = {this->fgamma*V[k]+r,a.p};

                        auto insertionPos = std::lower_bound(atoms.begin(), atoms.end(), newa,
                                                         atom<double>::comparator);
                        atoms.insert(insertionPos, newa);
                    }
                    ldistribution<double> d(atoms, false, true);
                    //newV[i]=this->fcrit(d,nothing());
                    const auto& crit = this->fcrit.nesting();
                    double x = crit(d,nothing());
                    if constexpr(optimize)
                    {
                        if(x > bestv)
                        {
                            bestv = x;
                            params.p[i]= a;
                        }
                        if(!this->fconstraint.nextfeasible(a,e))
                            break;
                    }
                    else
                    {
                        bestv = x;
                        break;
                    }
                }
                newV[i] = bestv;
            }
            if(j==0)
                error = dist(V,newV) / (1.0-this->fgamma) ;
            else
               error *= this->fgamma;
            V = newV;
            if(error < accuracy)
                break;
        }
        params.v = V;
        params.e = error;
    }

    viresult valueiteration(const finitevaluefunction& initialV,
                            double accuracy,
                            const computationparams& evalparams) const
    {
        viresult vr(*this);
        vr.v = initialV;
        iterate<true>(vr, accuracy, evalparams);
        return vr;
    }

    valuewitherror<finitevaluefunction> evaluate(
                         const finitevaluefunction& initialV,
                         finitepolicy& p,
                         double accuracy,
                         const computationparams& evalparams) const
    {
        viresult vr(*this);
        vr.p = p;
        vr.v = initialV;
        iterate<false>(vr, accuracy, evalparams);

        return { vr.v, vr.e };
    }


};


} // namespace

#endif // DP_HPP
