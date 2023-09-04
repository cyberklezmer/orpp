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

    class nestedproblem : public finitehomodpproblem<nestedcriterion<Criterion>,Statespace,ConstrainedActionSpace,
            Transition, Reward>
    {
    public:
        nestedproblem(const Criterion& crit,
                  const Statespace& state,
                  const ConstrainedActionSpace& constraint,
                  const Transition& transition,
                  const Reward& reward,
                  double gamma,
                  double maxreward) :
             finitehomodpproblem<nestedcriterion<Criterion>, Statespace,ConstrainedActionSpace,
                              Transition,Reward>
               (crit, state, constraint, transition, reward, gamma, maxreward)
        {
            static_assert(std::is_base_of<parametricriskmeasure<ldistribution<double>>,
                                                              Criterion>::value);
            static_assert(std::is_base_of<iteratedspace<typename Statespace::Element_t>,Statespace>::value);
            static_assert(std::is_base_of<constrainedspace<
                              typename Statespace::Element_t,
                              typename ConstrainedActionSpace::Element_t>,ConstrainedActionSpace>::value);
            static_assert(std::is_base_of<finitetransition<Statespace,ConstrainedActionSpace>,Transition>::value);
            static_assert(std::is_base_of<dpreward<Statespace,ConstrainedActionSpace>,Reward>::value);
        }
        using value = typename finitedpproblem<nestedcriterion<Criterion>,
             Statespace,ConstrainedActionSpace,Transition,Reward>::value;
        using computationparams = typename finitehomodpproblem<nestedcriterion<Criterion>,
           Statespace,ConstrainedActionSpace,Transition,Reward>::computationparams;

        void setriskaversion(double ra) { this->fcrit.nesting().setparam(ra); }
    };


    struct computationparams:
          public finitedpproblem<Criterion,Statespace,ConstrainedActionSpace,
            Transition, Reward>::computationparams
    {
        computationparams() : fpseudogradientmaxiters(100),
            fheuristicmaxiters(100) {}
        typename nestedproblem::computationparams fnestedparams;
        unsigned fpseudogradientmaxiters;
        unsigned fheuristicmaxiters;
    };


    overallriskproblem(const Criterion& crit,
              const Statespace& state,
              const ConstrainedActionSpace& constraint,
              const Transition& transition,
              const Reward& reward,
              double gamma,
              double maxreward ) :
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
                   const typename Problem::value& initvalue,
                   orpp::index initindex,
                   const rmtermination& term,
                   const typename nestedproblem::computationparams& pars ):
          fproblem(problem), fpolicy(policy), fconstant(constant),
          finitvalue(initvalue), finitindex(initindex), faccuracy(term.accuracy()),
          fparams(pars) {}
      Problem fproblem;
      finitepolicy fpolicy;
      double fconstant;
      typename Problem::value finitvalue;
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
                    const typename nestedproblem::value& initvalue,
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
                    sys::log() << std::endl << "error finding dif(iota)=0, constant=" << constant << std::endl;
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

    heuristicresult heuristic(index s0ind, double accuracy, const computationparams& params) const
    {
        sys::log() << "overallriskproblem::heuristic" << std::endl;

        double iota = this->fcrit.getparam();

        double valueofcrit = 0;
        nestedproblem problem(this->fcrit,this->fstatespace,
                              this->fconstraint,
                              this->ftransition,this->freward,
                              this->fgamma,
                              this->fmaxreward);
        finitepolicy bestp(problem);
        typename nestedproblem::value initV(problem,0);

        double error;
        double lastvalue = 0;
        for(unsigned i=0; i<params.fheuristicmaxiters; i++)
        {
            problem.setriskaversion(iota);
            sys::log() << "iteration " << i << " iota=" << iota << " cirt=" << valueofcrit << std::endl;

            auto vires = problem.valueiteration(initV,accuracy/3,params.fnestedparams);

            initV = vires.v;
            bestp = vires.p;

            statcounter cs = problem.evaluateraw( s0ind, {vires.p}, accuracy / 2.0, params.fnestedparams);

            valueofcrit = this->fcrit(cs.dist).x;

            if(i > 0)
            {
                error  = valueofcrit - lastvalue;
                if(error < accuracy)
                {
                    break;
                }
            }

            lastvalue = valueofcrit;

            iota = findiota(problem,vires.p,valueofcrit,vires.v, s0ind, accuracy, params.fnestedparams);
        }
        heuristicresult res(*this);
        res.p = bestp;
        res.v = valueofcrit;
        res.e = accuracy;
        res.iota = iota;
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
        sys::log() << "overallriskproblem::pseudogradientdescent" << std::endl;

        unsigned horizon = this->requiredhorizon(accuracy / 2);

        sys::log() << "Required horizon: " << horizon << std::endl;
        auto bestv = this->evaluatecrit(s0ind, ps, accuracy,params);

        for(unsigned i=0; i<params.fheuristicmaxiters; i++)
        {
            std::vector<finitepolicy> bestps = ps;
            sys::log() << "iteration " << i << " value=" << bestv.x << " ";
            for(unsigned x=0; x<ps.size(); x++)
            {
                for(unsigned y=0; y<ps[x].size(); y++)
                {
                    sys::log() << ps[x][y] << " ";
                }
                sys::log() << "- ";
            }
            sys::log() << std::endl;
            for(unsigned j=0; j<ps.size(); j++)
            {
                sys::log() << "time=" << j << std::endl;
                for(unsigned k=0; k<ps[j].size(); k++)
                {
                    sys::log() << "k=" << k;
                    for(unsigned n = 0; n <= 1; n++ )
                    {
                        std::vector<finitepolicy> p = ps;
                        if(n == 0 && p[j][k]>0)
                            p[j][k]--;
                        else if(n == 1 && p[j][k] < this->fconstraint.num()-1)
                            p[j][k]++;
                        else
                            continue;
                        for(unsigned x=0; x<ps.size(); x++)
                        {
                            sys::log()<< " -";
                            for(unsigned y=0; y<ps[x].size(); y++)
                            {
                                sys::log() << " " << p[x][y];
                            }
                        }
                        if(!this->fconstraint.isfeasible(p[j][k],k))
                        {
                            sys::log() << " infeasible" << std::endl;
                            continue;
                        }
                        auto v = this->evaluatecrit(s0ind, p, accuracy,params);

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


};





} // namespace

#endif // OVERALLRISKDP_HPP
