#include "orpp/dp.hpp"
#include "orpp/random.hpp"
#include "orpp/boostdist.hpp"
#include "orpp/simulations.hpp"

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
    std::cout << "Result:" << std::endl;
    std::cout << "V=";
    for(unsigned i=0; i<vires.v.size(); i++)
        std::cout << vires.v[i] << ",";
    std::cout << std::endl;
    for(unsigned j=0; j<=vires.p.size(); j++)
        std::cout << vires.p[j];
    std::cout << std::endl;

    std::cout << "Homogeneous enumeration." << std::endl;

    int horizon = log(accuracy / 2.0 * (1-problem.gamma())) / log( problem.gamma());

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

    std::cout << "Optimal homo: ";
    for(unsigned j=0; j<=overlinew; j++)
        std::cout << besthomo[j];

    std::cout << " " << besthv.x << "(" << besthv.sd << ")" << std::endl;

    std::cout << "Evaluating the value function by evaluateraw" << std::endl;


    testproblem::value Vraw(problem);
    for(unsigned i=0; i<= overlinew; i++)
    {
        auto sc = problem.evaluateraw(i,{besthomo}, horizon, accuracy / 2);
        Vraw[i] = sc.average();
        std::cout << i << " ~ " << Vraw[i] << std::endl;
    }


    std::cout << "Evaluating the value function by method evaluate " << std::endl;

    testproblem::value Ve(problem, 0.5);
    Ve = problem.evaluate(Ve,besthomo,accuracy).x;
    std::cout << "V=";
    for(unsigned i=0; i<Ve.size(); i++)
        std::cout << Ve[i] << ",";
    std::cout << std::endl;


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
    return !failed;
}



int main()
{
    testhomogeneity(0.005);
    return 0;


    std::vector<double> a = { 1,2,5,7};

    empiricaldistribution d(a);

/*
    for(double alpha = 0.0; alpha <= 0.99; alpha += 0.1)
    {
        CVaR<empiricaldistribution> c(alpha);
        MeanCVaR<empiricaldistribution> cc(alpha,0.3);
        empiricalMeanCVaR ccc(alpha,0.3);
        auto r = ccc(d);

        std::cout << alpha << " " << c(d,nothing()) <<
                    "=" << cc(d,nothing()) << "=" << r.x
                  << "(" << r.sd << ")" << std::endl;
    }

    return 0;

*/

    double alpha = 0;
    testproblem problem(alpha,0.7,0.85);


    double accuracy = 0.005;

    int horizon = log(accuracy / 2.0 * (1-problem.gamma())) / log( problem.gamma());

    std::cout << "horizon = " << horizon << std::endl;

    testproblem::policy besth(problem,1);
    realestimate besthv;

    for(unsigned i=1; i< overlinew; i++)
    {
        testproblem::policy policy(problem,1);
        for(unsigned j=0; j<i; j++)
            policy[j] = 0;
        auto sc = problem.evaluateraw(1,{policy}, horizon, accuracy / 2);
        empiricalMeanCVaR c(alpha,1);
        auto crit =  c(sc.dist);
        if(-crit.x > -besthv.x)
        {
            besthv = crit;
            besth = policy;
        }

        std::cout << i << ": " << crit.x << " (" << crit.sd << ") ";
        for(unsigned j=0; j<=overlinew; j++)
            std::cout << policy[j];
        std::cout << std::endl;
    }

    for(unsigned j=0; j<=overlinew; j++)
        std::cout << besth[j];
    std::cout << std::endl;

    std::cout << "Value iteration" << std::endl;

    testproblem::value V(problem,0.5);

    auto res=problem.valueiteration(V,0.01);
    std::cout << "Results" << std::endl;
    std::cout << "V=";
    for(unsigned i=0; i<V.size(); i++)
        std::cout << V[i] << ",";
    std::cout << std::endl;
    for(unsigned j=0; j<=overlinew; j++)
        std::cout << res.p[j];
    std::cout << std::endl;

    return 0;

/*    std::vector<testproblem::policy> bestp;
    statcounter bestv;
    enumerate(std::vector<testproblem::policy>(0),problem,horizon, accuracy / 2,
              4,
              bestp,bestv);

    std::cout << "Recursion: " << bestv.average() << " (" << bestv.averagestdev() << ")" << std::endl;

    auto dif = bestv.estaverage() - besthv;
    std::cout << "Dif: " << dif.x << stars(dif.twosidedsignificance());

    for(unsigned i=0; i<bestp.size(); i++)
    {
        for(unsigned j=0; j<=overlinew; j++)
            std::cout << bestp[i][j];
        std::cout << std::endl;
    }
    */
}
