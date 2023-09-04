#include "orpp/dp.hpp"
//#include "orpp/overallriskdp.hpp"
#include "orpp/boostdist.hpp"
#include <boost/math/tools/roots.hpp>
#include <fstream>

using namespace orpp;



constexpr int overlinew = 5;

class testactionspace :
        public integerspace, public constrainedspace<unsigned int,unsigned int>
{
public:
    testactionspace() : integerspace(0,1) {}
    virtual bool isfeasible(const unsigned int& e, const unsigned int& c) const
    {
        if(c == 0 && e==1)
            return false;
        else
            return true;
    }
};

class teststatespace : public integerspace
{
public:
    teststatespace() : integerspace(0,overlinew) {}
};

class testreward : public dpreward<teststatespace, testactionspace>
{
public:
    double operator() (const dpcondition<unsigned int, unsigned int>& x) const
    { return x.a; }
};

class testtransition: public finitetransition<teststatespace,testactionspace>
{
public:
    testtransition(probability ap) : fp(ap) {}
private:
    virtual unsigned natoms_is(const dpcondition<unsigned int,unsigned int>&) const
    { return overlinew+1; }
    virtual atom<unsigned int> atom_is(unsigned int i, const dpcondition<unsigned int,unsigned int>& c) const
    {
//std::cout << "atom " << i << " = ";
        assert(!(c.s == 0 && c.a==1));

        auto s = c.s - c.a;
        if(s==0)
        {
            if(i==0)
               return {0,1};
            else
               return {i,0};
        }
        if(i < s)
            return {i, 0};
        int bincensorship = overlinew - s;
        assert(bincensorship>=0);
        boost::math::binomial d(s,fp);
        unsigned int bini = i - s;
        if(bini<bincensorship)
        {
            assert(i<=overlinew);
            if(bini <= s)
                return { i,  boost::math::pdf( d, bini ) };
            else
                return { i, 0 };
        }
        else if(bini==bincensorship)
        {
            double p = 0;
            for(int j=bincensorship; j<=s; j++)
                p +=  boost::math::pdf( d, j);
//std::cout << ""
            assert(i<=overlinew);
            return { i, p };
//std::cout << "i=" << i << " p=" << p << std::endl;
        }
        else
            return { i ,0};
    }
    virtual unsigned int natoms_is(const dpcondition<int,int>&) const { return overlinew +1; };
    virtual bool is_sorted() const { return true; }
    probability fp;
};

using testcrit = CVaR<ldistribution<double>,true>;


class testproblem : public finitehomodpproblem<testcrit, teststatespace,
        testactionspace, testtransition, testreward>
{
public:
    testproblem(probability alpha, probability pincrease, double gamma) :
        finitehomodpproblem<testcrit, teststatespace,
                    testactionspace, testtransition, testreward>
          (testcrit(alpha),teststatespace(), testactionspace(),
           testtransition(pincrease), testreward(), gamma,1) {}
    void setriskaversion(probability alpha)
    {
        fcrit = testcrit(alpha);
    }
};

void enumerate(const std::vector<finitepolicy>& ps,
                      const testproblem& prob,
                      unsigned horizon,
                      double accuracy,
                      unsigned depth,
                      std::vector<finitepolicy>& bestp,
                      statcounter& bestv)
{
    if(ps.size()==depth)
    {
        statcounter c = prob.evaluateraw(1, ps, horizon, accuracy);
        if(bestp.size()==0 || c.average() > bestv.average())
        {
            bestp = ps;
            bestv = c;
        }
    }
    else
    {
        std::vector<testproblem::policy> p = ps;
        p.push_back(testproblem::policy(prob,1));
        for(unsigned i=1; i< overlinew; i++)
        {
            for(unsigned j=0; j<=i; j++)
                p[p.size()-1][j] = 0;
            if(p.size()<=2)
                std::cout << p.size();
            enumerate(p,prob,horizon,accuracy,depth, bestp, bestv);
            if(p.size()==1)
                std::cout << std::endl;
        }
    }
}

void testhomogeneity(std::ostream& out)
{
    double accuracy = 0.001;
    testproblem problem(0,0.7,0.80);

    out << "Value iteration" << std::endl;

    testproblem::value initV(problem,0.5);

    auto vires=problem.valueiteration(initV,0.01);

    out << "Homogeneous enumeration." << std::endl;

    int horizon = problem.requiredhorizon(accuracy / 2.0 );

    out << "required horizon = " << horizon << std::endl;

    testproblem::policy besthomo(problem);
    realestimate besthv;

    for(unsigned i=1; i< pow(2,besthomo.size()); i++)
    {
        testproblem::policy policy(problem,0);
        auto decimal = i;
        unsigned k=0;
        bool feasible = true;
        while (decimal != 0) {
            unsigned int p = decimal % 2;
            if(!problem.constraint().feasible(p,k))
                feasible = false;
            assert(k < policy.size());
            policy[k++] = p;
            decimal = decimal / 2;
          }
        out << i << ": ";
        for(unsigned j=0; j<=overlinew; j++)
            out << policy[j];
        if(!feasible)
        {
            out << " infeasible" << std::endl;
            continue;
        }

        auto sc = problem.evaluateraw(1,{policy}, accuracy / 2);
        auto crit =  sc.estaverage();
        if(crit.x > besthv.x)
        {
            besthv = crit;
            besthomo = policy;
        }

        out  << " " << crit.x << " (" << crit.sd << ")";
        out << std::endl;
    }


    out << "Evaluating the value function by evaluateraw" << std::endl;

    testproblem::value Vraw(problem);
    for(unsigned i=0; i<= overlinew; i++)
    {
        auto sc = problem.evaluateraw(i,{besthomo}, accuracy / 2);
        Vraw[i] = sc.average();
    }


    out << "Evaluating the value function by method evaluate" << std::endl;

    testproblem::value Ve(problem, 0.5);
    Ve = problem.evaluate(Ve,besthomo,accuracy).x;


    out << std::endl << "Value functions: " << std::endl;

    out << "value iteraton: ";
    for(unsigned i=0; i<vires.v.size(); i++)
        out << vires.v[i] << ",";
    out << std::endl;

    out << "evaluate: ";
    for(unsigned i=0; i<Ve.size(); i++)
        out << Ve[i] << ",";
    out << std::endl;

    out << "evaluateraw: ";
    for(unsigned i=0; i<= overlinew; i++)
    {
//        auto sc = problem.evaluateraw(i,{besthomo}, horizon, accuracy / 2);
//        Vraw[i] = sc.average();
        out << Vraw[i] << ",";
    }
    out  << std::endl;
    out  << std::endl;
    out  << "Errors:" << std::endl;

    out << "accuracy: " << accuracy << std::endl;
    out << "VI - Evaluate: " << testproblem::dist(Ve,vires.v) << std::endl;
    out << "VI - Vraw: " << testproblem::dist(Vraw,vires.v) << std::endl;

    out  << std::endl;
    out << "Policies" << std::endl;
    out << "value iteration: " ;

    for(unsigned j=0; j<vires.p.size(); j++)
        std::cout << vires.p[j];
    out << std::endl;

    out << "homogeneous enumeration: ";
    for(unsigned j=0; j<=overlinew; j++)
        out << besthomo[j];
    out << std::endl;
}

valuewitherror<double> pseudogradientdescent(
                      orpp::index s0ind,
                      double kappa,
                      std::vector<testproblem::policy>& ps,
                      const testproblem& prob,
                      double accuracy,
                      unsigned numiters)
{
    unsigned horizon = prob.requiredhorizon(accuracy / 2);
    valuewitherror<double> bestv={0,0};
    for(unsigned i=0; i<numiters; i++)
    {
        std::vector<testproblem::policy> bestps = ps;
        for(unsigned j=0; j<ps.size(); j++)
        {
std::cout << "j=" << j << std::endl;
            for(unsigned k=1; k<ps[j].size(); k++)
            {
std::cout << "k=" << k << std::endl;
                std::vector<testproblem::policy> p = ps;
                p[j][k]=1-p[j][k];
                statcounter sc = prob.evaluateraw(s0ind, p, accuracy);

                MeanCVaR<empiricaldistribution,true> c(kappa,1.0);

                auto v = c(sc.dist);
                if(v.x > bestv.x + bestv.sd)
                {
                    bestv = v;
                    bestps = p;
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

        for(unsigned j=0; j<ps.size(); j++)
        {
            for(unsigned k=0; k<ps[j].size(); k++)
            {
                std::cout << ps[j][k];
            }
            std::cout << " ";
        }
        std::cout << bestv.x << std::endl;
    }
    return bestv;
}


void testcvar()
{
/*
    std::vector<double> a = { 1,2,5,7};

    empiricaldistribution d(a);


    for(double alpha = 0.0; alpha <= 0.99; alpha += 0.1)
    {
        CVaR<empiricaldistribution> c(alpha);
        MeanCVaR<empiricaldistribution> cc(alpha,0.3);
        MeanCVaR<empiricaldistribur> ccc(alpha,0.3);
        auto r = ccc(d);

        std::cout << alpha << " " << c(d,nothing()) <<
                    (alpha == 0.0 ? "=" : ">") << cc(d,nothing()) << "=" << r.x
                  << "(" << r.sd << ")" << std::endl;
    }
*/
}

struct rmtermination  {
  bool operator() (double min, double max)  {
    return abs(min - max) <= faccuracy / 2;
  }

  rmtermination(double accuracy) : faccuracy(accuracy) {}
  double accuracy() const { return faccuracy; }
private:
    double faccuracy;

};

template <typename Problem>
struct rmdifference  {
  double operator() (double iota)
  {
      fproblem.setriskaversion(iota);
      auto r = fproblem.evaluate(finitvalue,fpolicy,faccuracy / 2);
      if(r.sd > faccuracy / 2)
          throw exception("rmdifference: failed to achieve desired accuracy");
      return r.x[finitindex] - fconstant;
  }
  rmdifference(const Problem& problem, const typename Problem::policy& policy,
               double constant,
               const typename Problem::value& initvalue,
               orpp::index initindex,
               const rmtermination& term):
      fproblem(problem), fpolicy(policy), fconstant(constant),
      finitvalue(initvalue), finitindex(initindex), faccuracy(term.accuracy()) {}
private:
  Problem fproblem;
  typename Problem::policy fpolicy;
  double fconstant;
  typename Problem::value finitvalue;
  orpp::index finitindex;
  double faccuracy;

};


void testheuristicodl(std::ostream& out)
{
    double accuracy = 0.001;
    double kappa = 0.6;
    testproblem problem(kappa,0.7,0.85);

    unsigned numiters = 3;
    testproblem::value initV(problem,0.5);

    orpp::index s0ind = 1;

    double iota = kappa;
    double constant = 0;
    out << "Starting the heuristic algorithm" << std::endl;
    testproblem::policy bestp(problem);
    for(unsigned i=0; i<numiters; i++)
    {
        problem.setriskaversion(iota);

        auto vires=problem.valueiteration(initV,accuracy);
        initV = vires.v;
        bestp = vires.p;

        statcounter cs = problem.evaluateraw( s0ind, {vires.p}, accuracy / 2.0 );

        MeanCVaR<empiricaldistribution,true> cv(kappa,1.0);

        constant = cv(cs.dist).x;

        out << "Policy: ";
        for(unsigned j=0; j<vires.p.size(); j++)
            out << vires.p[j];
        out << " Iota = " << iota;
        out << " with overall Cvar vaůie = " << constant << std::endl;

        rmtermination term(accuracy);

        rmdifference<testproblem> dif(problem,vires.p,constant,vires.v, s0ind ,term);

        using boost::math::tools::bisect;
        double from = 0;  // The solution must lie in the interval [from, to], additionally f(from) <= 0 && f(to) >= 0
        double to = kappa;


        std::pair<double, double> result =
           bisect(dif, from, to, term);
        iota = (result.first + result.second) / 2;  // = 0.381966...

    }

    out << "Evaluating the value function by method evaluate" << std::endl;

    problem.setriskaversion(iota);
    testproblem::value Ve = problem.evaluate(initV,bestp,accuracy).x;
    out << "Result is " << Ve[s0ind] << std::endl;
    out << "comparing to " << constant << std::endl;

    out << "Solving the original problem by pseudogradientdescent" << std::endl;

    std::vector<testproblem::policy> ps(2,bestp);


    auto bestv = pseudogradientdescent(s0ind,kappa,ps,problem,accuracy,3);
    out << "With result " << bestv.x << " compared to " << constant << std::endl;
}


int main_old()
{
    std::ofstream out("testdpprotocol.txt");
    if(!out)
        throw "cannot open protocol";
    testheuristicodl(std::cout);
//    testhomogeneity(std::cout);
    return 0;
}
