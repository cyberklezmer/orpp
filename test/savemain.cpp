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
        int maxval = overlinew-c.s;
        assert(maxval>=0);
        boost::math::binomial d(c.s,fp);
        if(i<maxval)
            return { i,  boost::math::pdf( d, i ) };
        else if(i==maxval)
        {
            double s = 0;
            for(unsigned j=maxval; j<=c.s; j++)
                s +=  boost::math::pdf( d, j );
            return { i, s };
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
   testproblem problem(0.95,0.1,0.9);
   finitepolicy policy(overlinew+1,0);
//   auto sc = problem.evaluate(overlinew / 2,policy, 20,0.01);
//   std::cout << sc.average() << " (" << sc.averagestdev() << ")" << std::endl;


/*    boost::math::binomial nd(1,0.5);

    std::cout << boost::math::pdf( nd, 0 ) << std::endl;
    std::cout << boost::math::pdf( nd, 1.0 ) << std::endl;
    std::cout << boost::math::pdf( nd, 0.5 ) << std::endl;*/
}

