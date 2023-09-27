#include <fstream>
#include "orpp/boostdist.hpp"
#include "orpp/overallriskdp.hpp"
#include "orpp/test/testdp.hpp"
#include "orpp/test/testoverallriskdp.hpp"

using namespace orpp;

class testactionspace :
        public integerspace, public constrainedspace<unsigned int,unsigned int>
{
public:
    testactionspace(unsigned amaxcons) : integerspace(0,amaxcons) {}
    unsigned maxcons() const { return num()-1; }
    virtual bool isfeasible(const unsigned int& e, const unsigned int& c) const
    {
        return e<=maxcons() && e <= c;
    }
};

class teststatespace : public integerspace
{
public:
    int nstates() const { return num(); }
    teststatespace(unsigned anstates) : integerspace(0,anstates-1)  {}
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
    testtransition(probability ap, unsigned nstates)
        : fp(ap), fnstates(nstates) {}
private:
    virtual unsigned natoms_is(const dpcondition<unsigned int,unsigned int>&) const
    { return fnstates; }
    virtual atom<unsigned int> atom_is(unsigned int i, const dpcondition<unsigned int,unsigned int>& c) const
    {
        assert(i < fnstates);
        assert(c.s < fnstates);
//std::cout << "atom " << i << " = ";
        assert(c.a <= c.s);
        assert(i < fnstates);
        auto s = static_cast<int>(c.s) - static_cast<int>(c.a);
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
        if(i < fnstates - 1)
            return { i,  boost::math::pdf( d, i - s ) };
        else
        {
            assert(i == fnstates - 1);
            double p = 0;
            for(int j=i-s; j<=s; j++)
                p +=  boost::math::pdf( d, j);
            return { i, p };
        }
    }
    virtual bool is_sorted() const { return true; }
    probability fp;
    unsigned fnstates;
};

using testcrit = CVaR<ldistribution<double>,true>;

class testproblem : public overallriskproblem<testcrit,
        teststatespace, testactionspace, testtransition, testreward>
{
public:
    testproblem(unsigned nstates, unsigned maxcons, probability alpha, probability pincrease, double gamma) :
        overallriskproblem<testcrit, teststatespace,
                    testactionspace, testtransition, testreward>
          (testcrit(alpha),teststatespace(nstates), testactionspace(maxcons),
           testtransition(pincrease,nstates), testreward(), gamma,1,2.0 / (1-alpha)) {}
};


accuracytestresult testoverall(const testproblem& problem,
                               orpp::index p0,
                               const std::vector<finitepolicy>& ps,
                               orpp::index s0ind, double accuracy,
                               unsigned testiters,
                               const typename testproblem::computationparams& params)
{
    sys::log() << "testoverall" << std::endl;
    return testevaluate<false>(problem, s0ind, p0, ps, accuracy, params, testiters);
}

class testhomoproblem: public testproblem::nestedproblem
{
public:
    testhomoproblem(unsigned nstates, unsigned maxcons,
                    double iota, double pincrease, double gamma) :
//        overallriskproblem<testcrit, teststatespace,
//                    testactionspace, testtransition, testreward>
          testproblem::nestedproblem(testcrit(iota),teststatespace(nstates),
                                     testactionspace(maxcons),
           testtransition(pincrease,nstates), testreward(), gamma, 1)
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


void test(unsigned nstates, unsigned maxcons,
            double kappa, double pincrease, double gamma,
             orpp::index s0ind, double accuracy,
             unsigned testiters,
             const testproblem::computationparams& params)
{
   testproblem problem(nstates, maxcons, kappa,pincrease,gamma);

   testproblem::heuristicresult res = problem.heuristic(s0ind,accuracy,params);
   testoverall(problem,res.p[s0ind],{res.p},s0ind,accuracy,testiters,params);

   testhomoproblem hp(nstates, maxcons, res.iota, pincrease,gamma);
   testhomo(hp,accuracy,s0ind,testiters,params.fnestedparams);
}

/* tbd
 *
 *     if constexpr(test)
    {
        if(!testlipshitzproperty(problem,s0ind,foopp,accuracy,pars,2))
           throw exception("nonlipschitz or nonomonotonous problem");
    }
*/

struct examineprogram
{
    unsigned nstates;
    unsigned maxcons;
    double kappa;
    double gamma;
    double accuracy;
    double pincrease = 0.7;
    orpp::index s0ind = 1;
//    unsigned testiters = 10;
    bool heuristic = false;
    bool taylorheuristic = false;
    bool pseudogradienthomo = false;
    bool riskneutral = false;
//    bool pesudogradient = false;
    bool pseudogradienthetero = false;
    bool enumerate = false;
    testproblem::computationparams pars;
    unsigned fmaxstatestoenum;
};

void examine(examineprogram p, std::ostream& report)
{
    report << p.nstates << "," << p.maxcons << "," << p.kappa << "," << p.gamma << "," << p.accuracy << ",";

    testproblem problem(p.nstates, p.maxcons, p.kappa, p.pincrease, p.gamma);

    finitepolicy startingp(problem,0);

    // enumeration
    unsigned tstart = sys::gettimems();

    if(p.riskneutral)
    {
         testhomoproblem hp(p.nstates, p.maxcons, 0, p.pincrease,p.gamma);
         finitevaluefunction initV(hp,0);
         testhomoproblem::viresult vires = hp.valueiteration(initV,p.accuracy,p.pars.fnestedparams);
         report << vires.p << "," << vires.v[p.s0ind] << ",";
    }
    else
        report << ",,";
    unsigned rnend = sys::gettimems();

    report << rnend - tstart << ",";


    if(p.heuristic)
    {
        try
        {
            testproblem::heuristicresult hres = problem.heuristic(p.s0ind,p.accuracy,p.pars);
            report << hres.p << "," << hres.v << "," << hres.iota << ",";
        }
        catch(const timelimitexception& e)
        {
            report << ",outoftime,,";
        }

    }
    else
        report << ",,,";
    unsigned hend = sys::gettimems();
    report  << hend - rnend << ",";

    if(p.taylorheuristic)
    {
        try
        {
            testproblem::heuristicresult hres = problem.taylorheuristic(p.s0ind,p.accuracy,startingp,p.pars);
            report << hres.p << ","
                   << hres.v << "," << hres.iota << ",";
        }
        catch(const timelimitexception& e)
        {
            report << ",outoftime,,";
        }

    }
    else
        report << ",,,";
    unsigned tend = sys::gettimems();
    report << tend - hend << ",";

    if(p.pseudogradienthomo)
    {
        try
        {
            testproblem::pgdhomoresult respg = problem.pseudogradientdescent(p.s0ind, startingp, p.accuracy, p.pars);
            report << respg.p;
            report << "," <<respg.v.x << ",";
        }
        catch(const timelimitexception& e)
        {
            report << ",outoftime,";
        }

    }
    else
        report << ",,";

    unsigned cdend = sys::gettimems();
    report << cdend - rnend << ",";

    if(p.enumerate) // tbd still only to log
    {
        if(pow(p.maxcons+1,p.nstates) > p.fmaxstatestoenum)
        {
            sys::logline() << "Too much states to enumerate" << std::endl;
            report << ",toomuchstates,";
        }
        else try
        {
            testproblem::enumresult res = problem.enumeratehomo(p.s0ind, p.accuracy, p.pars);

            sys::log() << "best of enumerate:" << res.p << " "
                       << res.v.x  << std::endl;
            report << res.p << "," << res.v.x << ",";
        }
        catch(const timelimitexception& e)
        {
            report << ",outoftime,";
        }
    }
    else
        report << ",,";

    unsigned eend = sys::gettimems();

    report << eend-cdend<< ",";


    if(p.pseudogradienthetero)
    {
        try
        {
            testproblem::heteropolicy heterop = { startingp[p.s0ind], {startingp, startingp, startingp} };

            testproblem::pgdheteroresult respg = problem.pseudogradientdescent(p.s0ind, heterop, p.accuracy, p.pars);
            report << respg.p.p0;
            for(unsigned k=0;k < respg.p.ps.size(); k++)
            {
                report << "-";
                report << respg.p.ps[k];
            }
            report << "," <<respg.v.x << ",";
        }
        catch(const timelimitexception& e)
        {
            report << ",outoftime,";
        }

    }
    else
        report << ",,";

    unsigned pgend = sys::gettimems();
    report << pgend - eend << ",";


    report << std::endl;

    // Calculating total time taken by the program.
    double time_taken = (pgend - tstart) / 1000.0;
    sys::logline() << "Time taken by program is : " << std::fixed
        << time_taken << std::setprecision(5);
    sys::log() << " sec " << std::endl;

}


int main(int argc, char *argv[])
{
    unsigned nthreads = 1;
    if(argc>1)
    {
        try
        {
            nthreads = std::stoul(argv[1]);
        }
        catch(...)
        {
            std::cout << "Error converting everyn string " << argv[1] << std:: endl;
            throw;
        }
    }

    sys::setlog(std::cout);
    sys::logline() << "Using " << nthreads << " threads." << std::endl;

    sys::setloglevel(0);
     std::ostringstream report;


     report << "nstates,maxcons,kappa,gamma,accuracy,"
           << "rnpolicy,rncrit,rntime,"
           << "hpolicy,hcrit,hlambda,htime,"
           << "tpolicy,tcrit,tlambda,ttime,"
           << "pghomopolicy,pghomocrit,pghomotime,"
           << "enumgpolicy,enumgcrit,enumegtime,"
           << "pgheteropolicy,pgheterocrit,pgheterotime,"
           << std::endl;

     report << std::setprecision(5);
    std::vector<double> kappas = { 0.6, 0.75, 0.9 };
    std::vector<double> gammas = { 0.85, 0.9, 0.95 };

    examineprogram p;
    testproblem::computationparams pars;

    p.nstates = 10;
    p.maxcons = 4;



    pars.fopttimelimit = pars.fpseudogradienttimelimit = pars.fenumtimelimit = 100000;

    pars.fthreadstouse = pars.fnestedtaylorparams.fthreadstouse = pars.fnestedonedparams.fthreadstouse
             = pars.fnestedparams.fthreadstouse = nthreads;
    pars.fthreadbatch = pars.fnestedtaylorparams.fthreadbatch = pars.fnestedonedparams.fthreadbatch
             = pars.fnestedparams.fthreadbatch = 3000;
    pars.fmaxevaliterations = 2000000;



    p.fmaxstatestoenum = 10000;
    p.pars = pars;
    p.pincrease = 0.7;
    p.s0ind = 1;
    p.riskneutral = true;
    p.heuristic = true;
    p.taylorheuristic = true;
    p.pseudogradienthomo = true;
    p.enumerate = true;
    p.pseudogradienthetero = true;
    p.kappa = kappas[2];
    p.gamma = gammas[1];
    p.accuracy = 0.05;

    examine(p, report); // tbd
    std::cout << report.str() << std::endl;
    return 0;
}
