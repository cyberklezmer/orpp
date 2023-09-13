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
        onedactionspace(const finitepolicy& maskpolicy, index freeindex)
            : ffreeindex(freeindex), fmaskpolicy(maskpolicy) {}
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


    struct computationparams:
          public finitedpproblem<Criterion,Statespace,ConstrainedActionSpace,
            Transition, Reward>::computationparams
    {
        computationparams() : fpseudogradientmaxiters(100),
            fheuristicmaxiters(100) {}
        typename nestedproblem::computationparams fnestedparams;
        typename nestedonedproblem::computationparams fnestedonedparams;
        unsigned fpseudogradientmaxiters;
        unsigned fheuristicmaxiters;
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
          (crit, state, constraint, transition, reward, gamma, maxreward)
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
              throw exception("rmdifference: failed to achieve desired accuracy");
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

    template <bool gradientdescent=false>
    heuristicresult heuristic(index s0ind, double accuracy, const computationparams& params) const
    {
        sys::logline() << "overallriskproblem::heuristic" << std::endl;

        double iota = this->fcrit.getparam();

        double valueofcrit = 0;
        finitepolicy bestp(this->fstatespace.num());
        finitevaluefunction initV(this->fstatespace.num(),0);

        double error;
        double lastvalue = 0;
        for(unsigned i=0; i<params.fheuristicmaxiters; i++)
        {
            nestedproblem problem(this->fcrit,
                                  this->fstatespace,
                                  this->fconstraint,
                                  this->ftransition,this->freward,
                                  this->fgamma,
                                  this->fmaxreward);

            if constexpr(gradientdescent)
            {
                if(i==0)
                {
                   sys::logline() << "Initial iteration" << std::endl;
                   problem.setriskaversion(iota);
                   auto vires = problem.valueiteration(initV,accuracy/3,params.fnestedparams);
                   bestp = vires.p;
                   initV = vires.v;
                   statcounter cs = problem.evaluateraw( s0ind, {vires.p}, accuracy / 2.0, params.fnestedparams);
                   valueofcrit = this->fcrit(cs.dist).x;
                   iota = findiota(problem,vires.p,valueofcrit,vires.v, s0ind, accuracy, params.fnestedparams);
                }

                for(unsigned j=0; j<this->fstatespace.num(); j++)
                {
                    onedactionspace as(bestp,j);
                    std::vector<double> adds;
                    auto p = bestp;
                    double maxadd = 0;
                    for(index k=0; k<as.num(); k++)
                    {
                        if(as.isfeasible(k,j))
                        {
                            p[j] = k;
                            statcounter cs = problem.evaluateraw( s0ind, {p}, accuracy / 2.0, params.fnestedparams);
                            auto c = this->fcrit(cs.dist).x;

                            auto lambda = findiota(problem,p,c,initV, s0ind, accuracy, params.fnestedparams);
                            double addition = 2 * this->maxreward() * this->frmlipshitzfactor * (this->fcrit.getparam() - fabs(lambda-iota) );
                            adds.push_back(addition);
                            if(addition > maxadd)
                                maxadd = addition;
                        }
                        else
                            adds.push_back(0);
                    }
                    onedreward r(adds);
                    nestedonedproblem onedproblem(this->fcrit,
                                          this->fstatespace,
                                          as,
                                          this->ftransition,r,
                                          this->fgamma,
                                          this->fmaxreward + maxadd);
                    onedproblem.setriskaversion(iota);
                    sys::logline() << "iteration " << i << "(" << j << ") " << " iota=" << iota << " cirt=" << valueofcrit;

                    auto vires = onedproblem.valueiteration(initV,accuracy/3,params.fnestedonedparams);


                    initV = vires.v;
                    bestp = vires.p;

                    for(unsigned i=0; i<bestp.size(); i++)
                        sys::log() << " " << bestp[i];
                    sys::log() << std::endl;

                    statcounter cs = onedproblem.evaluateraw( s0ind, {vires.p}, accuracy / 2.0, params.fnestedonedparams);
                    valueofcrit = this->fcrit(cs.dist).x;
                }
                if(i > 0)
                {
                    error  = valueofcrit - lastvalue;
                    if(error < accuracy)
                    {
                        break;
                    }
                }

                lastvalue = valueofcrit;

                iota = findiota(problem,bestp,valueofcrit,initV, s0ind, accuracy, params.fnestedparams);
            }
            else //constexpr !gradientdescent
            {
                problem.setriskaversion(iota);
                sys::logline() << "iteration " << i << " iota=" << iota << " cirt=" << valueofcrit << std::endl;

                auto vires = problem.valueiteration(initV,accuracy/3,params.fnestedparams);

                initV = vires.v;
                bestp = vires.p;

                statcounter cs = problem.evaluateraw( s0ind, {vires.p}, accuracy / 2.0, params.fnestedparams);

                valueofcrit = this->fcrit(cs.dist).x;

                if(i > 0)
                {
                    error  = valueofcrit - lastvalue;
                    if(error < accuracy)
                        break;
                }

                lastvalue = valueofcrit;
                iota = findiota(problem,vires.p,valueofcrit,vires.v, s0ind, accuracy, params.fnestedparams);
            }
        }
        heuristicresult res(*this);
        res.p = bestp;
        res.v = valueofcrit;
        res.e = accuracy;
        res.iota = iota;
        sys::logline() << "huristic ended: ";
        res.output(sys::log());
        return res;
    }

    valuewitherror<double> pseudogradientdescent(
                          orpp::index s0ind,
                          std::vector<finitepolicy>& ps,
                          double accuracy,
                          const computationparams& params
                           )
    {
        sys::logline() << "overallriskproblem::pseudogradientdescent,";

        unsigned horizon = this->requiredhorizon(accuracy / 2);

        sys::log() << "Required horizon: " << horizon << std::endl;
        auto bestv = this->evaluatecrit(s0ind, ps, accuracy,params);

        for(unsigned i=0; i<params.fheuristicmaxiters; i++)
        {
            std::vector<finitepolicy> bestps = ps;
            sys::logline() << "iteration " << i << " value=" << bestv.x << " ";
            for(unsigned x=0; x<ps.size(); x++)
                sys::log() << ps[x] << ",";
            sys::log() << std::endl;
            for(unsigned j=0; j<ps.size(); j++)
            {
                sys::logline() << "time=" << j << std::endl;
                for(unsigned k=0; k<ps[j].size(); k++)
                {
                    for(unsigned n = 0; n <= 1; n++ )
                    {
                        std::vector<finitepolicy> p = ps;
                        if(n == 0 && p[j][k]>0)
                            p[j][k]--;
                        else if(n == 1 && p[j][k] < this->fconstraint.num()-1)
                            p[j][k]++;
                        else
                            continue;
                        sys::logline() << "k=" << k << "," << "p:";
                        if(!this->fconstraint.isfeasible(p[j][k],k))
                        {
                            for(unsigned x=0; x<ps.size(); x++)
                              sys::log() <<  p[x] << ",";
                            sys::log() << " infeasible" << std::endl;
                            continue;
                        }
                        auto v = this->evaluatecrit(s0ind, p, accuracy,params);

                        for(unsigned x=0; x<ps.size(); x++)
                            sys::log() <<  p[x] << ",";
                        sys::log() << " = " << v.x << "(" << v.sd << ")";

                        if(v.x > bestv.x + bestv.sd)
                        {
                            sys::log() << "*";
                            bestv = v;
                            bestps = p;
                        }

                        sys::log() << std::endl;
                    }
                }
            }
            bool differs = false;
            for(unsigned j=0; j<ps.size(); j++)
            {
                if(!(ps[j] == bestps[j]))
                {
                    differs = true;
                    break;
                }
            }
            if(!differs)
                return bestv;
            ps = bestps;

            sys::log() << "bestv=" << bestv.x << std::endl;
        }
        return bestv;
    }

private:
    double frmlipshitzfactor;
};





} // namespace

#endif // OVERALLRISKDP_HPP
