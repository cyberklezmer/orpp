#include "orpp/overallriskdp.hpp"
#include <fstream>
#include "orpp/boostdist.hpp"
#include "orpp/test/testdp.hpp"

using namespace orpp;

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
    static constexpr int nstates = 10;
    teststatespace() : integerspace(0,nstates-1) {}
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
    { return teststatespace::nstates; }
    virtual atom<unsigned int> atom_is(unsigned int i, const dpcondition<unsigned int,unsigned int>& c) const
    {
        assert(i < teststatespace::nstates);
        assert(c.s < teststatespace::nstates);
//std::cout << "atom " << i << " = ";
        assert(!(c.s == 0 && c.a==1));
        assert(c.a <= 1);
        assert(i < teststatespace::nstates);
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
        if(i > 2 * s)
            return {i, 0};

//        int bincensorship = teststatespace::nstates - 1 - s;
//        assert(bincensorship>=0);
        boost::math::binomial d(s,fp);
        if(i < teststatespace::nstates - 1)
            return { i,  boost::math::pdf( d, i - s ) };
        else
        {
            assert(i == teststatespace::nstates - 1);
            double p = 0;
            for(int j=i-s; j<=s; j++)
                p +=  boost::math::pdf( d, j);
            return { i, p };
        }
    }
    virtual bool is_sorted() const { return true; }
    probability fp;
};

using testcrit = CVaR<ldistribution<double>,true>;

class testproblem : public overallriskproblem<testcrit,
        teststatespace, testactionspace, testtransition, testreward>
{
public:
    testproblem(probability alpha, probability pincrease, double gamma) :
        overallriskproblem<testcrit, teststatespace,
                    testactionspace, testtransition, testreward>
          (testcrit(alpha),teststatespace(), testactionspace(),
           testtransition(pincrease), testreward(), gamma,1,2.0 / (1-alpha)) {}
};


accuracytestresult testoverall(const testproblem& problem,
                               const std::vector<finitepolicy>& ps,
                               orpp::index s0ind, double accuracy,
                               unsigned testiters,
                               const typename testproblem::computationparams& params)
{
    sys::log() << "testoverall" << std::endl;
    return testevaluate<false>(problem, s0ind, ps ,accuracy, params, testiters);
}

class testhomoproblem: public testproblem::nestedproblem
{
public:
    testhomoproblem(double iota, double pincrease, double gamma) :
//        overallriskproblem<testcrit, teststatespace,
//                    testactionspace, testtransition, testreward>
          testproblem::nestedproblem(testcrit(iota),teststatespace(), testactionspace(),
           testtransition(pincrease), testreward(), gamma,1)
    {

    }
};

accuracytestresult testhomo(const testhomoproblem& problem,
                            double accuracy, orpp::index s0ind,
                            unsigned testiters,
                            const testhomoproblem::computationparams& cp
                            )
{
    finitevaluefunction initV(problem,0);
    auto res = problem.valueiteration(initV,accuracy,cp);

    return testevaluatehomo(problem, s0ind, res.p ,accuracy, cp, testiters);
}


void testhomotime(const testhomoproblem& problem,
                            double accuracy, orpp::index s0ind,
                            unsigned testiters,
                            const testhomoproblem::computationparams& cp
                            )
{
    finitevaluefunction initV(problem,0);
    for(unsigned i=0; i<10; i++)
       problem.valueiteration(initV,accuracy,cp);
}


template <bool test>
void proceed(double kappa, double pincrease, double gamma,
             orpp::index s0ind, double accuracy,
             unsigned testiters,
             const testproblem::computationparams& params)
{
   testproblem problem(kappa,pincrease,gamma);

// testproblem::heuristicresult res = problem.heuristic<false>(s0ind,accuracy,params);
std::vector<orpp::index> fv = { 0,0,0,1,1,1,1,1,1,1};
finitepolicy foo(fv);
    if constexpr(test)
    {
//        testoverall(problem,{res.p},s0ind,accuracy,testiters,params);
// testoverall(problem,{foo},s0ind,accuracy,testiters,params);


        //testhomoproblem hp(res.iota, pincrease,gamma);
//        testhomo(hp,accuracy,s0ind,testiters,params.fnestedparams);

    testhomoproblem hp(0.4, pincrease,gamma);
        testhomotime(hp,accuracy,s0ind,testiters,params.fnestedparams);
//std::vector<finitepolicy> ps(2,foo);
//    std::vector<finitepolicy> ps(2,res.p);

    //auto resp = problem.pseudogradientdescent(s0ind, ps, accuracy, params);
     

    }
std::vector<finitepolicy> ps(2,foo);
//    std::vector<finitepolicy> ps(2,res.p);

  auto resp = problem.pseudogradientdescent(s0ind, ps, accuracy, params);
//    sys::log() << "result heuristic = " << res.v << std::endl;
//    sys::log() << "result pseudo = " << resp.x << std::endl;
}


void measure(int threads)
{
    time_t start, end;

    time(&start);

    double pincrease = 0.7;
    double gamma = 0.85;
    double kappa = 0.6;// 0.6;
    double accuracy = 0.03;
    orpp::index s0ind = 1;
    unsigned testiters = 10;

    testproblem::computationparams pars;
    pars.fthreadstouse = pars.fnestedonedparams.fthreadstouse
            = pars.fnestedparams.fthreadstouse = threads;
    pars.fthreadbatch = pars.fnestedonedparams.fthreadbatch
            = pars.fnestedparams.fthreadbatch = 3000;
    pars.fmaxevaliterations = 2000000;

    proceed<true>(kappa, pincrease, gamma, s0ind, accuracy, testiters, pars);

    time(&end);

    // Calculating total time taken by the program.
    double time_taken = double(end - start);
    std::cout << "Time taken by program is : " << std::fixed
        << time_taken << std::setprecision(5);
    std::cout << " sec " << std::endl;
}


int main()
{
    sys::setlog(std::cout);
    sys::setloglevel(0);
    // measure(40);
    measure(0);
    return 0;
}

