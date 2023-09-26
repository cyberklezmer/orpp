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

   testproblem::heuristicresult res = problem.heuristic<false>(s0ind,accuracy,params);
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
    bool coordinatedescent = false;
    bool riskneutral = false;
//    bool pesudogradient = false;
    bool pseudogradienthetero = false;
    bool enumerate = false;
    testproblem::computationparams pars;

};

void examine(examineprogram p, std::ostream& report)
{
    report << p.nstates << "," << p.maxcons << "," << p.kappa << "," << p.gamma << "," << p.accuracy << ",";

    testproblem problem(p.nstates, p.maxcons, p.kappa, p.pincrease, p.gamma);

    finitepolicy startingp(problem,0);

    // enumeration
    unsigned tstart = sys::timems();

    if(p.riskneutral)
    {
         testhomoproblem hp(p.nstates, p.maxcons, 0, p.pincrease,p.gamma);
         finitevaluefunction initV(hp,0);
         testhomoproblem::viresult vires = hp.valueiteration(initV,p.accuracy,p.pars.fnestedparams);
         report << vires.p << "," << vires.v[p.s0ind] << ",";
    }
    else
        report << ",,";
    unsigned rnend = sys::timems();

    report << rnend - tstart << ",";


    if(p.heuristic)
    {
        testproblem::heuristicresult hres = problem.heuristic<false>(p.s0ind,p.accuracy,p.pars);
        report << hres.p << "," << hres.v << "," << hres.iota << ",";
    }
    else
        report << ",,,";
    unsigned hend = sys::timems();
    report  << hend - rnend << ",";

    if(p.taylorheuristic)
    {
        testproblem::heuristicresult hres = problem.taylorheuristic(p.s0ind,p.accuracy,startingp,p.pars);
        report << hres.p << ","
               << hres.v << "," << hres.iota << ",";
    }
    else
        report << ",,,";
    unsigned tend = sys::timems();
    report << tend - hend << ",";


    if(p.coordinatedescent)
    {
        testproblem::heuristicresult hres = problem.heuristic<true>(p.s0ind,p.accuracy,p.pars);
        report << hres.p << "," << hres.v << "," << hres.iota << ",";
    }
    else
        report << ",,,";
    unsigned cdend = sys::timems();
    report << cdend - rnend << ",";

    if(p.enumerate) // tbd still only to log
    {
        sys::logline() << "enumerate" << std::endl;
        finitepolicy besthomo(problem);
        double besthomov = 0;
        for(unsigned i=1; i< pow(p.maxcons+1,p.nstates); i++)
        {
            finitepolicy policy(problem,0);
            auto decimal = i;
            unsigned k=0;
            bool feasible = true;
            while (decimal != 0) {
                unsigned int z = decimal % (p.maxcons+1);
                if(!problem.constraint().feasible(z,k))
                    feasible = false;
                assert(k < policy.size());
                policy[k++] = z;
                decimal = decimal / (p.maxcons+1);
              }
            if(!feasible)
            {
                sys::logline() << i << ": "  << policy << " infeasible" << std::endl;
                continue;
            }

            auto res = problem.evaluatecrit(p.s0ind,policy, p.accuracy / 2, p.pars);
            sys::logline() << i << ": "  << policy << " " <<res.x << "(" << res.sd << ")";
            if(res.x > besthomov)
            {
                besthomov = res.x;
                besthomo = policy;
                sys::log() << "*";
            }
            sys::log() << std::endl;
        }
        sys::log() << "best of enumerate:" << besthomo << " "
                   << besthomov << std::endl;
        report << besthomo << "," << besthomov << ","; ;
    }
    else
        report << ",,";

    unsigned eend = sys::timems();

    report << eend-cdend<< ",";


    if(p.pseudogradienthetero)
    {
        testproblem::heteropolicy heterop = { startingp[p.s0ind], {startingp, startingp, startingp} };

        testproblem::pgdresult respg = problem.pseudogradientdescent(p.s0ind, heterop, p.accuracy, p.pars);
        report << respg.p.p0;
        for(unsigned k=0;k < respg.p.ps.size(); k++)
        {
            report << "-";
            report << respg.p.ps[k];
        }
        report << "," <<respg.v.x << ",";
    }
    else
        report << ",,";

    unsigned pgend = sys::timems();
    report << pgend - eend << ",";


    report << std::endl;

    // Calculating total time taken by the program.
    double time_taken = (pgend - tstart) / 1000.0;
    sys::logline() << "Time taken by program is : " << std::fixed
        << time_taken << std::setprecision(5);
    sys::log() << " sec " << std::endl;

}


int main()
{
    sys::setlog(std::cout);
    sys::setloglevel(0);
     std::ostringstream report;

    testproblem::computationparams pars;

    pars.fthreadstouse = pars.fnestedtaylorparams.fthreadstouse = pars.fnestedonedparams.fthreadstouse
             = pars.fnestedparams.fthreadstouse = 8;
    pars.fthreadbatch = pars.fnestedtaylorparams.fthreadbatch = pars.fnestedonedparams.fthreadbatch
             = pars.fnestedparams.fthreadbatch = 3000;
    pars.fmaxevaliterations = 2000000;

    examineprogram p;
    p.pars = pars;


    p.pincrease = 0.7;
    p.s0ind = 1;
//    unsigned testiters = 10;
    p.riskneutral = true;
    p.heuristic = true;
    p.taylorheuristic = true;
    p.coordinatedescent = true;
    p.enumerate = false;
    p.pseudogradienthetero = true;


     report << "rnpolicy,rncrit,rntime,"
           << "nsatates,maxcons,kappa,gamma,accuracy"
           << "hpolicy,hcrit,hlambda,htime,"
           << "tpolicy,tcrit,tlambda,ttime,"
           << "cdolicy,cdcrit,cdlambda,cdtime,"
           << "pgpolicy,pgcrit,pgtime"
           << "egpolicy,egcrit,egtime"
           << std::endl;

     report << std::setprecision(5);
    std::vector<double> kappas = { 0.6, 0.75, 0.9 };
    std::vector<double> gammas = { 0.85, 0.9, 0.95 };

    p.nstates = 8;
    p.maxcons = 3;
    p.kappa = kappas[2];
    p.gamma = gammas[1];
    p.accuracy = 0.002;

    examine(p, report); // tbd
    std::cout << report.str() << std::endl;
    return 0;
}
