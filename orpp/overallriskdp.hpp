#ifndef OVERALLRISKDP_HPP
#define OVERALLRISKDP_HPP

#include <boost/math/tools/roots.hpp>
#include "orpp/dp.hpp"
#include "orpp/random.hpp"

namespace orpp
{


template<typename Criterion,
          typename Statespace,
          typename ConstrainedActionSpace,
          typename Transition,
          typename Reward
          >
class overallriskproblem : public finitedpproblem<Criterion,
        Statespace,ConstrainedActionSpace,Transition,Reward>
{
public:
    using Statespace_t = Statespace;
    using ConstrainedActionSpace_t = ConstrainedActionSpace;

    template <typename A, typename R>
    class nestedproblembase : public finitehomodpproblem<nestedcriterion<Criterion>,Statespace,A,
            Transition, R>
    {
    public:
        nestedproblembase(const Criterion& crit,
                  const Statespace& state,
                  const A& constraint,
                  const Transition& transition,
                  const R& reward,
                  double gamma,
                  double maxreward) :
             finitehomodpproblem<nestedcriterion<Criterion>, Statespace,A,
                              Transition,R>
               (crit, state, constraint, transition, reward, gamma, maxreward)
        {
            static_assert(std::is_base_of<parametricriskmeasure<ldistribution<double>>,
                                                              Criterion>::value);
            static_assert(std::is_base_of<iteratedspace<typename Statespace::Element_t>,Statespace>::value);
            static_assert(std::is_base_of<constrainedspace<
                              typename Statespace::Element_t,
                              typename A::Element_t>,A>::value);
            static_assert(std::is_base_of<finitetransition<Statespace,A>,Transition>::value);
            static_assert(std::is_base_of<dpreward<Statespace,A>,R>::value);
        }

//        using computationparams = typename finitehomodpproblem<nestedcriterion<Criterion>,
//                           Statespace,ConstrainedActionSpace,Transition,Reward>::computationparams;
        void setriskaversion(double ra) { this->fcrit.nesting().setparam(ra); }
    };

    using nestedproblem = nestedproblembase<ConstrainedActionSpace,Reward>;


    class taylorreward : public Reward
    {
    public:
        taylorreward(const finitepolicy& base,
                     const std::vector<double>& gradient,
                     const Reward& orig) : Reward(orig),
            fbase(base), fgradient(gradient)
        {
            assert(fbase.size() == fgradient.size());
        }
        double operator() (const dpcondition<unsigned int,
                                        typename ConstrainedActionSpace::Element_t>& x) const
        {
            assert(x.s < fgradient.size());
            assert(x.s < fbase.size());
            return std::max(0.0,Reward::operator()(x)+(x.a-fbase[x.s]) * fgradient[x.s]);
        }
        double maxgradadd(const ConstrainedActionSpace& space) const
        {
            double m = 0;
            for(unsigned i=0; i<fgradient.size(); i++)
            {
                double x = 0;
                if(fgradient[i] < 0)
                {
                    typename ConstrainedActionSpace::Element_t a;
                    bool exists = space.firstfeasible(a,i);
                    assert(exists);
                    if(!exists)
                        throw exception("maxgradadd: At least one feasible element needed");
                    assert(a <= fbase[i]);
                    x = - (fbase[i]-a) * fgradient[i];
                }
                else if(fgradient[i] > 0)
                {
                    typename ConstrainedActionSpace::Element_t a;
                    bool exists = space.lastfeasible(a,i);
                    assert(exists);
                    if(!exists)
                        throw exception("maxgradadd: At least one feasible element needed");
                    assert(a >= fbase[i]);
                    x = (fbase[i]-a) * fgradient[i];
                }
                m = std::max(m,x);
            }
            return m;
        }

    private:
        finitepolicy fbase;
        std::vector<double> fgradient;
    };


    class onedreward : public Reward
    {
    public:
        onedreward(std::vector<double> additions) :
            fadditions(additions) {}
        double operator() (const dpcondition<unsigned int, unsigned int>& x) const
        { return Reward::operator()(x)+fadditions[x.a]; }
    private: 
        std::vector<double> fadditions;
    };

    class onedactionspace : public ConstrainedActionSpace
    {
    public:
        onedactionspace(const finitepolicy& maskpolicy, index freeindex,
                        const ConstrainedActionSpace& space)
            : ConstrainedActionSpace(space),
              ffreeindex(freeindex), fmaskpolicy(maskpolicy) {}
    private:
        index ffreeindex;
        finitepolicy fmaskpolicy;
    public:
        virtual bool isfeasible(const typename ConstrainedActionSpace::Element_t& a,
                                const typename Statespace::Element_t& s ) const
        {
            assert(s < fmaskpolicy.size());
            if(s == ffreeindex)
                return ConstrainedActionSpace::isfeasible(a,s);
            else
                return a == fmaskpolicy[s];
        }
    };

    using nestedonedproblem = nestedproblembase<onedactionspace,onedreward>;

    using nestedtaylorproblem = nestedproblembase<ConstrainedActionSpace,taylorreward>;

    struct computationparams:
          public finitedpproblem<Criterion,Statespace,ConstrainedActionSpace,
            Transition, Reward>::computationparams
    {
        computationparams() : fpseudogradientmaxiters(100),
            fheuristicmaxiters(100), fopttimelimit(std::numeric_limits<timems>::max()) {}
        typename nestedtaylorproblem::computationparams fnestedtaylorparams;
        typename nestedproblem::computationparams fnestedparams;
        typename nestedonedproblem::computationparams fnestedonedparams;
        unsigned fpseudogradientmaxiters;
        unsigned fheuristicmaxiters;
        timems fopttimelimit;
    };


    overallriskproblem(const Criterion& crit,
              const Statespace& state,
              const ConstrainedActionSpace& constraint,
              const Transition& transition,
              const Reward& reward,
              double gamma,
              double maxreward ,
              double rmlipshitzfactor) :
        finitedpproblem<Criterion, Statespace,ConstrainedActionSpace,
                         Transition,Reward>
          (crit, state, constraint, transition, reward, gamma, maxreward),
           frmlipshitzfactor(rmlipshitzfactor)
    {}

private:
    struct rmtermination  {
      bool operator() (double min, double max)  {
        return abs(min - max) <= faccuracy / 2;
      }

      rmtermination(double accuracy) : faccuracy(accuracy) {}
      double accuracy() const { return faccuracy; }
    private:
        double faccuracy;

    };
private:
    template <typename Problem>
    struct rmdifference  {
      double operator() (double iota)
      {
          fproblem.setriskaversion(iota);
          auto r = fproblem.evaluate(finitvalue,fpolicy,faccuracy / 2, fparams);
          if(r.sd > faccuracy / 2)
              throw exception("rmdifference: failed to achieve desired accuracy.");
          return r.x[finitindex] - fconstant;
      }
      rmdifference(const Problem& problem,
                   const finitepolicy& policy,
                   double constant,
                   const finitevaluefunction& initvalue,
                   orpp::index initindex,
                   const rmtermination& term,
                   const typename nestedproblem::computationparams& pars ):
          fproblem(problem), fpolicy(policy), fconstant(constant),
          finitvalue(initvalue), finitindex(initindex), faccuracy(term.accuracy()),
          fparams(pars) {}
      Problem fproblem;
      finitepolicy fpolicy;
      double fconstant;
      finitevaluefunction finitvalue;
      orpp::index finitindex;
      double faccuracy;
      typename nestedproblem::computationparams fparams;
    };

public:

    struct heuristicresult
    {
        heuristicresult(const overallriskproblem<Criterion,Statespace,ConstrainedActionSpace,
                 Transition,Reward>& pr) : p(pr),v(0), e(0) {}
        finitepolicy p;
        double v;
        double e;
        double iota;
        void output(std::ostream& o)
        {
            o << "Policy:";
            for(unsigned i=0; i<p.size(); i++)
                o << " " << p[i];
            o << std::endl;
            o << "Value: " << v;
            o << std::endl;
            o << "Error: " << e << std::endl;
            o << "iota: " << iota << std::endl;
        }
    };

    double findiota(const nestedproblem& problem,
                    const finitepolicy& policy,
                    double constant,
                    const finitevaluefunction initvalue,
                    orpp::index initindex,
                    double accuracy,
                    const typename nestedproblem::computationparams& pars
                    ) const
    {
        rmtermination term(accuracy / 2);

        rmdifference<nestedproblem> dif(problem,policy,constant,initvalue, initindex ,term, pars);

        using boost::math::tools::bisect;
        double from = 0;  // The solution must lie in the interval [from, to], additionally f(from) <= 0 && f(to) >= 0
        double to = this->fcrit.getparam();
        assert(from <= to);
        double dfrom = dif(from);
        double dto = dif(to);
        if(fabs(dfrom)<accuracy)
            return from;
        else if(fabs(dto)<accuracy)
            return to;
        else if(from < to)
        {
            if(dfrom < 0)
                throw exception("Error finding nested r.a. parameter: dif(from) < 0");
            else if(dto > 0)
                throw exception("Error finding nested r.a. parameter: dif(to) > 0");
            else
            {
                try
                {
                    std::pair<double, double> result =
                       bisect(dif, from, to, term);
                    return (result.first + result.second) / 2;
                }
                catch(...)
                {
                    sys::logline() << std::endl << "error finding dif(iota)=0, constant=" << constant << std::endl;
                    const unsigned nsteps = 5;
                    for(unsigned i=0; i<=nsteps; i++)
                    {
                        double x = from+ (to-from)/nsteps * i;
                        sys::log() << "dif(" << x << ")=" << dif(x) << std::endl;
                    }
                    throw;
                }
             }
        }
        else
            throw exception("Error finding nested r.a. parameter: from > to");
    }

    heuristicresult heuristic(index s0ind, double accuracy,
                              const computationparams& params) const
    {
        sys::logline() << "overallriskproblem::heuristic" << std::endl;

        auto st = sys::gettimems();

        nestedproblem problem(this->fcrit,
                              this->fstatespace,
                              this->fconstraint,
                              this->ftransition,this->freward,
                              this->fgamma,
                              this->fmaxreward);
        problem.setriskaversion(0);

        auto vires = problem.valueiteration(finitevaluefunction(problem,0),accuracy,params.fnestedparams);

        finitepolicy bestp = vires.p;
        finitepolicy lastp = bestp;

        sys::logline() << "RN iteration: p="
                       << vires.p << ", expectation=" << vires.v[s0ind]<< std::endl;


        double iota = this->fcrit.getparam();

        finitevaluefunction initV(this->fstatespace.num(),0);

        double lastvalueofcrit = 0;
        double resultingvalueofcrit = 0;



        for(unsigned i=0; ; i++)
        {


            problem.setriskaversion(iota);

            auto vires = problem.valueiteration(initV,accuracy/3,params.fnestedparams);

            initV = vires.v;

            statcounter cs = problem.evaluateraw( s0ind, {vires.p}, accuracy / 2.0, params.fnestedparams);

            double valueofcrit = this->fcrit(cs.dist).x;

            sys::logline() << "iteration " << i << " iota=" << iota << ", p="
                           << vires.p << ", crit=" << valueofcrit << std::endl;


            if(i > 0)
            {
                if(valueofcrit <= lastvalueofcrit + accuracy
                        || i==params.fheuristicmaxiters-1)
                {
                    resultingvalueofcrit = lastvalueofcrit;
                    bestp = lastp;
                    break;
                }
            }
            lastp=vires.p;
            lastvalueofcrit = valueofcrit;

            iota = findiota(problem,vires.p,valueofcrit,vires.v, s0ind, accuracy, params.fnestedparams);
            auto t = sys::gettimems();
            if(t - st > params.fopttimelimit)
            {
                sys::logline() << "Exiting for time reasons" << std::endl;
                throw timelimitexception(params.fopttimelimit);
            }
        }
        heuristicresult res(*this);
        res.p = bestp;
        res.v = resultingvalueofcrit;
        res.e = accuracy;
        res.iota = iota;
        sys::logline() << "huristic ended: ";
        res.output(sys::log());
        return res;
    }

    heuristicresult taylorheuristic(index s0ind, double accuracy,
                                    finitepolicy& initp,
                              const computationparams& params) const
    {
        sys::logline() << "overallriskproblem::taylorheuristic" << std::endl;

        finitepolicy candp(initp);
        finitepolicy lastp(initp);
        finitevaluefunction initV(*this,0);
        finitevaluefunction taylorinitV(*this,0);
        double lastvalueofcrit = 0;
        double iota;
        auto st = sys::gettimems();

        for(unsigned i=0; ; i++)
        {
            nestedproblem nproblem(this->fcrit,
                                  this->fstatespace,
                                  this->fconstraint,
                                  this->ftransition,
                                  this->freward,
                                  this->fgamma,
                                  this->fmaxreward);

            double valueofcrit =  this->evaluatecrit( s0ind, candp, accuracy / 2.0, params).x;

            iota = findiota(nproblem,candp,valueofcrit,initV, s0ind, accuracy, params.fnestedparams);

            sys::logline() << "iteration " << i << " iota=" << iota << ", p="
                         << candp << ", crit=" << valueofcrit << std::endl;

            if((i>0 && valueofcrit <= lastvalueofcrit + accuracy)
                    || i==params.fheuristicmaxiters-1)
                break;

            nproblem.setriskaversion(iota);

            std::vector<double> grad;
            for(unsigned s=0; s<candp.size(); s++)
            {
                auto p=candp;
                typename ConstrainedActionSpace::Element_t a0 = p[s];
                typename ConstrainedActionSpace::Element_t a1 = a0;
                if(this->fconstraint.nextfeasible(a1,s))
                {
                    assert(a1 > a0);
                }
                else if(this->fconstraint.previousfeasible(a1,s))
                {
                    assert(a1 < a0);
                }
                else
                {
                    grad.push_back(0);
                    continue;
                }
                assert(a0 != a1);
                p[s] = a1;
                double rhoalpha = this->evaluatecrit( s0ind, p, accuracy / 2.0, params).x;
                auto taylorV = nproblem.evaluate(initV,p,accuracy / 2, params.fnestedparams).x;
                double rhoeta = taylorV[s0ind];
 //std::cout << "s=" << s << " r0=" << rhoeta << " r1=" << rhoalpha
//          << " d=" << (static_cast<double>(a1) - static_cast<double>(a0))
//          << " p:" << p << std::endl;
                grad.push_back((1-this->fgamma)*(rhoalpha-rhoeta) / (a1 - a0));
                taylorinitV = taylorV;
            }
            taylorreward r(candp,grad, this->reward());
            double addition = (1-this->fgamma) * r.maxgradadd(this->fconstraint);

            sys::logline(0) << "grad=";
            for(unsigned k=0; k<grad.size(); k++)
                sys::log() << grad[k] << " ";
            sys::log() << std::endl;

            nestedtaylorproblem taylorproblem(this->fcrit,
                                  this->fstatespace,
                                  this->fconstraint,
                                  this->ftransition,r,
                                  this->fgamma,
                                  this->fmaxreward + addition); // should be reset

            taylorproblem.setriskaversion(iota);

            auto vires = taylorproblem.valueiteration(initV,accuracy/3,params.fnestedtaylorparams);




            lastp = candp;
            candp = vires.p;
            taylorinitV = vires.v;
            lastvalueofcrit = valueofcrit;
            auto t = sys::gettimems();
            if(t - st > params.fopttimelimit)
            {
                sys::logline() << "Exiting for time reasons" << std::endl;
                throw timelimitexception(params.fopttimelimit);
            }

        }
        heuristicresult res(*this);
        res.p = lastp;
        res.v = lastvalueofcrit;
        res.e = accuracy;
        res.iota = iota;
        sys::logline() << "taylor huristic ended: ";
        res.output(sys::log());
        return res;
    }


    struct heuristicplusresult
    {
        heuristicresult hres;
        typename finitedpproblem<Criterion,
                Statespace,ConstrainedActionSpace,Transition,Reward>::
             pgdhomoresult pgres;
    };

    heuristicplusresult heuristicplus(index s0ind, double accuracy,
                              const computationparams& params) const
    {
        // enumeration

        sys::logline() << "overallriskdpproblem::heuristicplus" << std::endl;

        heuristicresult hres = this->heuristic(s0ind,accuracy,params);

        auto pgres = this->pseudogradientdescenthomo(s0ind, hres.p, accuracy, params,1U);

        return { hres, pgres };
    }



    double lipschitzconstant() const
    {
        return frmlipshitzfactor * this->maxreward() / (1 - this->gamma());
    }

    void setriskaversion(double ra) { this->crit().setparam(ra); }

private:
    double frmlipshitzfactor;
};





} // namespace

#endif // OVERALLRISKDP_HPP

