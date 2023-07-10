#include "orpp/dp.hpp"
#include "orpp/random.hpp"
#include "orpp/boostdist.hpp"
#include "orpp/simulations.hpp"
#include <boost/math/tools/roots.hpp>


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


class testproblem : public finitedpproblem<testcrit, teststatespace,
        testactionspace, testtransition, testreward>
{
public:
    testproblem(probability alpha, probability pincrease, double gamma) :
        finitedpproblem<testcrit, teststatespace,
                    testactionspace, testtransition, testreward>
          (testcrit(alpha),teststatespace(), testactionspace(),
           testtransition(pincrease), testreward(), gamma) {}
    void setriskaversion(probability alpha)
    {
        fcrit = testcrit(alpha);
    }
};

void enumerate(const std::vector<testproblem::policy>& ps,
                      const testproblem& prob,
                      unsigned horizon,
                      double accuracy,
                      unsigned depth,
                      std::vector<testproblem::policy>& bestp,
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

bool testhomogeneity(double accuracy = 0.001)
{
    testproblem problem(0,0.7,0.80);

    std::cout << "Value iteration" << std::endl;

    testproblem::value initV(problem,0.5);

    auto vires=problem.valueiteration(initV,0.01);

    std::cout << "Homogeneous enumeration." << std::endl;

    int horizon = problem.horizon(accuracy / 2.0, 1.0 / (1-problem.gamma()) );

    std::cout << "horizon = " << horizon << std::endl;

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
        std::cout << i << ": ";
        for(unsigned j=0; j<=overlinew; j++)
            std::cout << policy[j];
        if(!feasible)
        {
            std::cout << " infeasible" << std::endl;
            continue;
        }

        auto sc = problem.evaluateraw(1,{policy}, horizon, accuracy / 2);
        auto crit =  sc.estaverage();
        if(crit.x > besthv.x)
        {
            besthv = crit;
            besthomo = policy;
        }

        std::cout  << " " << crit.x << " (" << crit.sd << ")";
        std::cout << std::endl;
    }


    std::cout << "Evaluating the value function by evaluateraw" << std::endl;

    testproblem::value Vraw(problem);
    for(unsigned i=0; i<= overlinew; i++)
    {
        auto sc = problem.evaluateraw(i,{besthomo}, horizon, accuracy / 2);
        Vraw[i] = sc.average();
    }


    std::cout << "Evaluating the value function by method evaluate" << std::endl;

    testproblem::value Ve(problem, 0.5);
    Ve = problem.evaluate(Ve,besthomo,accuracy).x;

/*
    unsigned stages = 3;

    std::cout << stages << " - step enumeration." << std::endl;

    std::vector<testproblem::policy> bestp;
    statcounter bestv;


    std::vector<testproblem::policy> initialp;
    enumerate(initialp,problem,horizon, accuracy / 2,
              stages,
              bestp,bestv);

    std::cout << "Optimal enumeration: " << bestv.average() << " (" << bestv.averagestdev() << ")" << std::endl;

    auto dif = bestv.estaverage() - besthv;
    std::cout << "Dif: " << dif.x << stars(dif.twosidedsignificance()) << std::endl;

    std::cout << "Same optimal value test ";
    if(dif.twosidedsignificance() == orpp::enotsignificant)
        std::cout << "passed" << std::endl;
    else
        std::cout << "failed" << std::endl;

    std::cout << "Same policy test " << std::endl;

    std::cout << "Homo policy" << std::endl;
    for(unsigned j=0; j<=overlinew; j++)
        std::cout << besthomo[j];
    std::cout << std::endl;


    std::cout << "Hetero policy" << std::endl;
    bool failed = false;
    for(unsigned i=0; i<bestp.size(); i++)
    {
        for(unsigned j=0; j<=overlinew; j++)
        {
            std::cout << bestp[i][j];
            if(i>0 && bestp[i][j] != besthomo[j]) // first policy can differ as we start from a single state
                failed = true;
        }
        std::cout << std::endl;
    }
    std::cout << (failed ? "failed" : "passed") << std::endl;
*/

    std::cout << "Value functions: ";
    std::cout << std::endl;

    std::cout << "Value iteraton: ";
    for(unsigned i=0; i<vires.v.size(); i++)
        std::cout << vires.v[i] << ",";
    std::cout << std::endl;

    std::cout << "Evaluate: ";
    for(unsigned i=0; i<Ve.size(); i++)
        std::cout << Ve[i] << ",";
    std::cout << std::endl;

    std::cout << "Evaluateraw: ";
    for(unsigned i=0; i<= overlinew; i++)
    {
//        auto sc = problem.evaluateraw(i,{besthomo}, horizon, accuracy / 2);
//        Vraw[i] = sc.average();
        std::cout << Vraw[i] << ",";
    }
    std::cout  << std::endl;

    std::cout << "Accuracy: " << accuracy << std::endl;
    std::cout << "Ve - Evaluate: " << testproblem::dist(Ve,vires.v) << std::endl;
    std::cout << "Ve - Vraw: " << testproblem::dist(Vraw,vires.v) << std::endl;

    std::cout  << std::endl;
    std::cout << "Policies" << std::endl;
    std::cout << "Value iteration: " ;

    for(unsigned j=0; j<vires.p.size(); j++)
        std::cout << vires.p[j];
    std::cout << std::endl;

    std::cout << "Homogeneous enumeration: ";
    for(unsigned j=0; j<=overlinew; j++)
        std::cout << besthomo[j];
    std::cout << std::endl;
    return true;
}


void testcvar()
{
    std::vector<double> a = { 1,2,5,7};

    empiricaldistribution d(a);


    for(double alpha = 0.0; alpha <= 0.99; alpha += 0.1)
    {
        CVaR<empiricaldistribution> c(alpha);
        MeanCVaR<empiricaldistribution> cc(alpha,0.3);
        empiricalMeanCVaR ccc(alpha,0.3);
        auto r = ccc(d);

        std::cout << alpha << " " << c(d,nothing()) <<
                    (alpha == 0.0 ? "=" : ">") << cc(d,nothing()) << "=" << r.x
                  << "(" << r.sd << ")" << std::endl;
    }
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


int main()
{
    double accuracy = 0.001;
    double kappa = 0.6;
    testproblem problem(kappa,0.7,0.85);

    unsigned numiters = 10;
    testproblem::value initV(problem,0.5);

    int horizon = problem.horizon(accuracy / 2.0, 1.0 / (1-problem.gamma()) );
    std::cout << "Horizon=" << horizon << std::endl;
    orpp::index s0ind = 1;

    double iota = kappa;
    double constant = 0;
    testproblem::policy bestp(problem);
    for(unsigned i=0; i<numiters; i++)
    {
        problem.setriskaversion(iota);

        testproblem::viresult vires=problem.valueiteration(initV,accuracy);
        initV = vires.v;
        bestp = vires.p;

        statcounter cs = problem.evaluateraw( s0ind, {vires.p}, horizon, accuracy / 2.0 );

        empiricalMeanCVaR<true> cv(kappa,1.0);

        constant = cv(cs.dist).x;

        std::cout << "Policy: ";
        for(unsigned j=0; j<vires.p.size(); j++)
            std::cout << vires.p[j];
        std::cout << " Iota = " << iota;
        std::cout << " with Cvar " << constant << std::endl;

        rmtermination term(accuracy);

        rmdifference<testproblem> dif(problem,vires.p,constant,vires.v, s0ind ,term);

        using boost::math::tools::bisect;
        double from = 0;  // The solution must lie in the interval [from, to], additionally f(from) <= 0 && f(to) >= 0
        double to = kappa;


        std::pair<double, double> result =
           bisect(dif, from, to, term);
        iota = (result.first + result.second) / 2;  // = 0.381966...





    }

    std::cout << "Evaluating the value function by method evaluate" << std::endl;

    problem.setriskaversion(iota);
    testproblem::value Ve = problem.evaluate(initV,bestp,accuracy).x;
    std::cout << "Result is " << Ve[s0ind] << std::endl;
    std::cout << "comparing to " << constant << std::endl;

    return 0;

}
