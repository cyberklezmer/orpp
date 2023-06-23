#include "orpp/dp.hpp"
#include "orpp/random.hpp"
#include "orpp/boostdist.hpp"
#include "orpp/simulations.hpp"

using namespace orpp;



constexpr int overlinep = 10;
constexpr int overlinew = 5;

class testactionspace : public integeriteratedspace<int>
{
public:
    testactionspace() : integeriteratedspace(0, overlinew) {}
    virtual bool isfeasible(const int& e, const int& c) const
    {
        if(c == 0 && e==1)
            return false;
        else
            return true;
    }
};

class teststatespace : public integeriteratedspace<nothing>
{
public:
    teststatespace() : integeriteratedspace<nothing>(0,overlinew) {}
};

class testreward : public reward<teststatespace, testactionspace>
{
public:
    double operator() (const dpcondition<int, int>& x) const
    { return x.a; }
};

class testtransition: public fdistribution<unsigned int,dpcondition<int,int>>
{
public:
    testtransition(probability ap) : fp(ap) {}
private:
    virtual atom<unsigned int> atom_is(unsigned int i, const dpcondition<int,int>& c) const
    {
//std::cout << "atom " << i << " = ";
        assert(!(c.s == 0 && c.a==1));
        auto s = c.s - c.a;
        if(i < s)
            return {i, 0};
        int bincensorship = overlinew - s;
        assert(bincensorship>=0);
        boost::math::binomial d(s,fp);
        unsigned int bini = i - s;
        if(bini<bincensorship)
        {
            assert(i<=overlinew);
            return { i,  boost::math::pdf( d, bini ) };
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
    virtual bool is_sorted() const { return true; }
    probability fp;
};

using testcrit = CVaR<fdistribution<double,nothing>>;


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
        statcounter c = prob.evaluate(1, ps, horizon, accuracy);
        if(bestp.size()==0 || c.average() > bestv.average())
        {
            bestp = ps;
            bestv = c;
        }
    }
    else
    {
        std::vector<finitepolicy> p = ps;
        p.push_back(finitepolicy(overlinew+1,1));
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

int main()
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
                    "=" << cc(d,nothing()) << "=" << r.x
                  << "(" << r.sd << ")" << std::endl;
    }

    return 0;

    testproblem problem(0.95,0.7,0.85);


    double accuracy = 0.01;

    int horizon = log(accuracy / 2.0 * (1-problem.gamma())) / log( problem.gamma());

    std::cout << "horizon = " << horizon << std::endl;



    finitepolicy besth(overlinew+1,1);
    double besthv = 0;

    for(unsigned i=1; i< overlinew; i++)
    {
        finitepolicy policy(overlinew+1,1);
        for(unsigned j=0; j<=i; j++)
            policy[j] = 0;
        auto sc = problem.evaluate(1,{policy}, horizon, accuracy / 2);
        if(sc.average() > besthv)
        {
            besthv = sc.average();
            besth = policy;
        }

        std::cout << i << ": " << sc.average() << " (" << sc.averagestdev() << ")" << std::endl;
    }

    for(unsigned j=0; j<=overlinew; j++)
        std::cout << besth[j];
    std::cout << std::endl;


    std::vector<finitepolicy> bestp;
    statcounter bestv;
    enumerate(std::vector<finitepolicy>(0),problem,horizon, accuracy / 2,
              4,
              bestp,bestv);

    std::cout << "Recursion: " << bestv.average() << " (" << bestv.averagestdev() << ")" << std::endl;
    for(unsigned i=0; i<bestp.size(); i++)
    {
        for(unsigned j=0; j<=overlinew; j++)
            std::cout << bestp[i][j];
        std::cout << std::endl;
    }
}
