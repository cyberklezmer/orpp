#include "orpp/dp.hpp"
#include "orpp/random.hpp"
#include "orpp/boostdist.hpp"

using namespace orpp;



constexpr int overlinep = 10;
constexpr int overlinew = 10;

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

int main()
{
    std::vector<atom<int>> a = { {0,0.5}, {1,0.5} };

    ldistribution<int> d(a);
//    d.atoms<true>(a,nothing());

    CVaR<ldistribution<int>> c(0.95);
//    c(d,nothing());

/*    for(double alpha = 0.0; alpha <= 0.99; alpha += 0.1)
    {
        CVaR<ldistribution<int>> c(alpha);

       std::cout << alpha << " " << c(d,nothing()) << std::endl;;
    }*/

    testproblem problem(0.95,0.7,0.90);

//tbd test var.

    double accuracy = 0.01;

    int horizon = log(accuracy / 2.0 * (1-problem.gamma())) / log( problem.gamma());

    std::cout << "horizon = " << horizon << std::endl;



    for(unsigned i=1; i< overlinew; i++)
    {
        finitepolicy policy(overlinew+1,1);
        for(unsigned j=0; j<=i; j++)
            policy[j] = 0;
        auto sc = problem.evaluate(1,policy, horizon, accuracy / 2);
        std::cout << i << ": " << sc.average() << " (" << sc.averagestdev() << ")" << std::endl;
    }

}
