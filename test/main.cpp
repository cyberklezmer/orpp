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
    static constexpr int nstates = 5;
    teststatespace() : integerspace(0,nstates) {}
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
    { return teststatespace::nstates+1; }
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
        int bincensorship = teststatespace::nstates - s;
        assert(bincensorship>=0);
        boost::math::binomial d(s,fp);
        unsigned int bini = i - s;
        if(bini<bincensorship)
        {
            assert(i<=teststatespace::nstates);
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
            assert(i<=teststatespace::nstates);
            return { i, p };
//std::cout << "i=" << i << " p=" << p << std::endl;
        }
        else
            return { i ,0};
    }
    virtual unsigned int natoms_is(const dpcondition<int,int>&) const { return teststatespace::nstates +1; };
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

using testhomoproblem = typename testproblem::nestedproblem;

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

template <bool test>
void proceed(double kappa, double pincrease, double gamma,
             orpp::index s0ind, double accuracy,
             unsigned testiters,
             const testproblem::computationparams& params)
{
    testproblem problem(kappa,pincrease,gamma);

    //testproblem::heuristicresult res = problem.heuristic<true>(s0ind,accuracy,params);
std::vector<orpp::index> fv = { 0,0,0,1,1,1};
finitepolicy foo(fv);
    if constexpr(test)
    {
//        testoverall(problem,{res.p},s0ind,accuracy,testiters,params);
testoverall(problem,{foo},s0ind,accuracy,testiters,params);


//        testhomoproblem hp(res.iota,pincrease,gamma);
//        testhomo(hp,accuracy,s0ind,testiters,params.fnestedparams);


    }
//std::vector<finitepolicy> ps(2,foo);
//  std::vector<finitepolicy> ps(2,res.p);

//    auto resp = problem.pseudogradientdescent(s0ind, ps, accuracy, params);
    //sys::log() << "result heuristic = " << res.v << std::endl;
//    sys::log() << "result pseudo = " << resp.x << std::endl;
}


void measure(int threads)
{
    time_t start, end;

    time(&start);

    sys::setlog(std::cout);
    sys::setloglevel(2);

    double pincrease = 0.7;
    double gamma = 0.85;
    double kappa = 0.6;// 0.6;
    double accuracy = 0.003;
    orpp::index s0ind = 1;
    unsigned testiters = 1;

    testproblem::computationparams pars;
    pars.fthreadstouse = threads;
    pars.fthreadbatch = 10000;
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
    measure(10);
    measure(0);
    return 0;
}
