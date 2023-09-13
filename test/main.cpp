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


void test(double kappa, double pincrease, double gamma,
             orpp::index s0ind, double accuracy,
             unsigned testiters,
             const testproblem::computationparams& params)
{
   testproblem problem(kappa,pincrease,gamma);

   testproblem::heuristicresult res = problem.heuristic<false>(s0ind,accuracy,params);
   testoverall(problem,{res.p},s0ind,accuracy,testiters,params);

   testhomoproblem hp(res.iota, pincrease,gamma);
   testhomo(hp,accuracy,s0ind,testiters,params.fnestedparams);
}


void examine(double kappa, double gamma,  double accuracy, std::ostream& report)
{

    double pincrease = 0.7;
    orpp::index s0ind = 1;
    unsigned testiters = 10;

    testproblem::computationparams pars;
    pars.fthreadstouse = pars.fnestedonedparams.fthreadstouse
            = pars.fnestedparams.fthreadstouse = 48;
    pars.fthreadbatch = pars.fnestedonedparams.fthreadbatch
            = pars.fnestedparams.fthreadbatch = 3000;
    pars.fmaxevaliterations = 2000000;

    testproblem problem(kappa,pincrease,gamma);

    report << kappa << "," << gamma << "," << accuracy << ",";

    unsigned tstart = sys::timems();
    testproblem::heuristicresult gdres = problem.heuristic<true>(s0ind,accuracy,pars);
    unsigned gdend = sys::timems();
    report << gdres.p << ","
           << gdres.v << "," << gdres.iota << ","
           << gdend - tstart << ",";

    testproblem::heuristicresult hres = problem.heuristic<false>(s0ind,accuracy,pars);
    unsigned hend = sys::timems();
    report << hres.p << ","
           << hres.v << "," << hres.iota << ","
           << hend - gdend << ",";



    testhomoproblem hp(0, pincrease,gamma);
    finitevaluefunction initV(hp,0);
    testhomoproblem::viresult vires = hp.valueiteration(initV,accuracy,pars.fnestedparams);

    unsigned eend = sys::timems();


    report << vires.p << ","
           << vires.v[s0ind] << ","
           << eend - hend << ",";

    std::vector<finitepolicy> ps(2,gdres.p);

    testproblem::pgdresult respg = problem.pseudogradientdescent(s0ind, ps, accuracy, pars);

    unsigned pgend = sys::timems();

    for(unsigned k=0;; k++)
    {
        report << respg.ps[k];
        if(k== ps.size()-1)
            break;
        report << "-";
    }
    report << respg.v.x << "," << pgend - eend; ;


    // Calculating total time taken by the program.
    double time_taken = (pgend - tstart) / 1000.0;
    sys::logline() << "Time taken by program is : " << std::fixed
        << time_taken << std::setprecision(5);
    sys::log() << " sec " << std::endl;

    report << std::endl;
}


int main()
{
    sys::setlog(std::cout);
    sys::setloglevel(0);
//    std::ofstream report("report.csv");
//    if(!report)
//        throw exception("Cannot open report.csv");
   std::ostringstream report;
   report << "kappa,gamma,"
           << "gdpolicy,gdcrit,accuracygdlambda,gdtime,"
           << "hpolicy,hcrit,hlambda,htime,"
           << "epolicy,ecrit,etime,"
           << "pgpolicy,pgcrit,pgtime" << std::endl;
    report << std::setprecision(5);
     std::vector<double> kappas = { 0.6, 0.75, 0.9 };
    std::vector<double> gammas = { 0.85, 0.9, 0.95 };

    examine(kappas[2],gammas[1], 0.05, report);
    std::cout << report.str() << std::endl;
    return 0;
}
