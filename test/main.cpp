#include <fstream>
#include "orpp/boostdist.hpp"
#include "orpp/overallriskdp.hpp"
#include "orpp/test/testdp.hpp"
// #include "orpp/test/testoverallriskdp.hpp"

using namespace orpp;

class invactionspace :
        public integerspace, public constrainedspace<unsigned int,unsigned int>
{
public:
    invactionspace(unsigned maxinv, unsigned lot) :
        integerspace(0,maxinv / lot ), flot(lot), fmaxinv(maxinv)
    {
    }
    virtual bool isfeasible(const unsigned int& a, const unsigned int& s) const
    {
        return flot * a <= s && flot * a <= fmaxinv;
    }
    unsigned lot() const { return flot; }
private:
    unsigned flot;
    unsigned fmaxinv;
};

class invstatespace : public integerspace
{
public:
    int nstates() const { return num(); }
    invstatespace(unsigned maxinv) : integerspace(0,2*maxinv), fmaxinv(maxinv)  {}
    unsigned maxinv() const { return fmaxinv; }
private:
    unsigned fmaxinv;
};

class invreward : public dpreward<invstatespace, invactionspace>
{
public:
    invreward(unsigned amaxinv) : fmaxinv(amaxinv) {}
    double operator() (const dpcondition<unsigned int, unsigned int>& x) const
    {
        return std::max(0, static_cast<int>(x.s) - static_cast<int>(fmaxinv));
    }
private:
    unsigned fmaxinv;
};

class invtransition: public finitetransition<invstatespace,invactionspace>
{
public:
    invtransition(unsigned maxinv, unsigned lot, probability ap, probability apcrash)
        : fp(ap), fpcrash(apcrash), fmaxinv(maxinv), fnstates(maxinv*2+1), flot(lot)
    {}
private:
    virtual unsigned natoms_is(const dpcondition<unsigned int,unsigned int>&) const
       { return fnstates; }
    virtual atom<unsigned int> atom_is(unsigned int i, const dpcondition<unsigned int,unsigned int>& c) const
    {
        assert(i < fnstates);
        assert(c.s < fnstates);
        unsigned toconsume = std::max(0, static_cast<int>(c.s) - static_cast<int>(fmaxinv));
        unsigned toinvest = c.s - toconsume;
        assert(toinvest <= fmaxinv);
        unsigned rflots = c.a;
        unsigned rfinv = c.a * flot;
        assert(rfinv <= fmaxinv);
        unsigned riskinv = toinvest - rfinv;
        unsigned rfandcresult = rfinv + rflots;

        if(riskinv == 0)
        {
            if(i==rfandcresult)
                return {i,1};
            else
                return {i,0};
        }
        if(i < rfandcresult)
            return {i,0};
        if(i == rfandcresult )
            return {i,fpcrash};
        if(i < rfandcresult+riskinv)
            return {i,0};
        if(i > rfandcresult+2 * riskinv )
            return {i,0};
        boost::math::binomial d(riskinv,fp);
        double p = (1-fpcrash) * boost::math::pdf( d, i - (rfandcresult+riskinv) );
        return {i,p};
    }
    virtual bool is_sorted() const { return true; }
    probability fp;
    probability fpcrash;
    unsigned fmaxinv;
    unsigned fnstates;
    unsigned flot;
};

using critcvar = CVaR<ldistribution<double>,true>;
class critmcv: public MeanCVaR<ldistribution<double>,true>
{
public:
    critmcv(double alpha) : MeanCVaR(0.95, alpha) {}
};

template <typename Crit>
class invproblem : public overallriskproblem<Crit,
        invstatespace, invactionspace, invtransition, invreward>
{
public:
    invproblem(unsigned maxinv, unsigned lot, probability alpha,
               probability pincrease, double gamma, probability pcrash) :
        overallriskproblem<Crit, invstatespace,
                    invactionspace, invtransition, invreward>
          (Crit(alpha),invstatespace(maxinv), invactionspace(maxinv,lot),
           invtransition(maxinv,lot,pincrease,pcrash), invreward(maxinv), gamma,maxinv/lot,2.0 / (1-alpha)) {}
    void doenumerateex(const finitepolicy& ap,
                           orpp::index s0ind,
                           double accuracy,
                           typename finitedpproblem<Crit,invstatespace, invactionspace, invtransition, invreward>::enumresult& thebest,
                           const typename finitedpproblem<Crit,invstatespace, invactionspace, invtransition, invreward>::computationparams& params,
                           timems st) const
    {
        orpp::index si = ap.size();
        if(si < this->statespace().num())
        {
            unsigned s = this->statespace()[si];

            finitepolicy p = ap;
            p.push_back(0); // all the same what we push

            auto mi = this->statespace().maxinv();
            bool justcopying = si > mi;

            orpp::index ai;
            if(justcopying)
                ai = ap[mi];
            else
                ai = 0;
            unsigned a;
            if(!this->constraint().first(a))
                throw exception("At least one action has to exist");
            bool atleastonefeasible = false;
            for(; ; ai++)
            {
                if(this->constraint().feasible(a,s))
                {
                    atleastonefeasible = true;
                    p[si] = ai;
                    doenumerateex(p,s0ind,accuracy, thebest,params,st);
                }
                if(justcopying || !this->constraint().next(a))
                   break;
            }
            if(!atleastonefeasible)
                throw exception("At least one action has to be feasible");
            if(!justcopying)
               assert(ai == this->constraint().num()-1);
        }
        else
        {
           auto t = sys::gettimems();
           if(t - st > params.fenumtimelimit)
           {
                sys::logline() << "Exiting for time reasons" << std::endl;
                throw timelimitexception(params.fenumtimelimit);
           }
           auto res = this->evaluatecrit(s0ind,ap, accuracy / 2, params);
           sys::logline() << ap << ": " <<res.x << "(" << res.sd << ")";
           if(res.x > thebest.v.x)
           {
                thebest.v = {res.x, accuracy / 2};
                thebest.p = ap;
                sys::log() << "*";
            }
            sys::log() << std::endl;
        }

    }
    typename finitedpproblem<Crit,invstatespace, invactionspace, invtransition, invreward>::enumresult enumeratehomoex(orpp::index s0ind,
                         double accuracy,
                         const typename finitedpproblem<Crit,invstatespace, invactionspace, invtransition, invreward>::computationparams& params)
    {
        sys::logline() << "finitedpproblem::enumerate" << std::endl;
        auto st = sys::gettimems();

       typename finitedpproblem<Crit,invstatespace, invactionspace, invtransition, invreward>::enumresult
                thebest =  { finitepolicy(*this), {0, 0} };

        finitepolicy p(0U);

        doenumerateex(p,s0ind,accuracy,thebest,params, st);


        sys::logline() << "finitedpproblem::enumerate ended" << std::endl;
        return thebest;
    }



};

template <typename Crit>
class invhomoproblem: public invproblem<Crit>::nestedproblem
{
public:
    invhomoproblem(unsigned maxinv, unsigned lot, probability alpha,
                   probability pincrease, double gamma, probability pcrash)  :
          invproblem<Crit>::nestedproblem(Crit(alpha),invstatespace(maxinv), invactionspace(maxinv,lot),
                                              invtransition(maxinv,lot,pincrease,pcrash), invreward(maxinv), gamma,
                                    static_cast<double>(maxinv))
    {
    }
};




class testactionspace :
        public integerspace, public constrainedspace<unsigned int,unsigned int>
{
public:
    testactionspace(unsigned amaxcons, unsigned numstates) : integerspace(0,amaxcons),
      fnumstates(numstates)
    {
        if(amaxcons < (numstates-1) - maxfinalstate())
            throw exception("Too little maxcons");
    }
    unsigned maxcons() const { return num()-1; }
    virtual bool isfeasible(const unsigned int& a, const unsigned int& s) const
    {
        if(s > maxfinalstate())
            return a == s-maxfinalstate();
        else
           return a<=maxcons() && a <= s;
    }
private:
    unsigned maxfinalstate() const { return (fnumstates-1) / 2; }

   unsigned fnumstates;
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
    testtransition(probability ap, unsigned nstates, probability apcrash)
        : fp(ap), fnstates(nstates), fpcrash(apcrash) {}
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
        if(i == 0)
            return {i, fpcrash};
        // now i>0, s>0
        if(i < s)
            return {i, 0};
        if(i > 2 * s)
            return {i, 0};

        boost::math::binomial d(s,fp);
        if(i <= fnstates - 1)
        {
            double p = (1-fpcrash) * boost::math::pdf( d, i - s );
            return { i,  p};
        }
        else
        {
            throw exception("should be obsolete");
            assert(i == fnstates - 1);
            double p = 0;
            for(int j=i-s; j<=s; j++)
                p += (1-fpcrash) * boost::math::pdf( d, j);
            return { i, p };
        }
    }
    virtual bool is_sorted() const { return true; }
    probability fp;
    unsigned fnstates;
    probability fpcrash;
};


template <typename Crit>
class testproblem : public overallriskproblem<Crit,
        teststatespace, testactionspace, testtransition, testreward>
{
public:
    testproblem(unsigned nstates, unsigned maxcons, probability alpha, probability pincrease, double gamma, probability pcrash) :
        overallriskproblem<Crit, teststatespace,
                    testactionspace, testtransition, testreward>
          (Crit(alpha),teststatespace(nstates), testactionspace(maxcons,nstates),
           testtransition(pincrease,nstates, pcrash), testreward(), gamma,maxcons,2.0 / (1-alpha)) {}
};


template <typename Crit>
accuracytestresult testoverall(const testproblem<Crit>& problem,
                               orpp::index p0,
                               const std::vector<finitepolicy>& ps,
                               orpp::index s0ind, double accuracy,
                               unsigned testiters,
                               const typename testproblem<Crit>::computationparams& params)
{
    sys::log() << "testoverall" << std::endl;
    return testevaluate<false>(problem, s0ind, p0, ps, accuracy, params, testiters);
}

template <typename Crit>
class testhomoproblem: public testproblem<Crit>::nestedproblem
{
public:
    testhomoproblem(unsigned nstates, unsigned maxcons,
                    double iota, double pincrease, double gamma, double pcrash) :
//        overallriskproblem<testcrit, teststatespace,
//                    testactionspace, testtransition, testreward>
          testproblem<Crit>::nestedproblem(Crit(iota),teststatespace(nstates),
                                     testactionspace(maxcons, nstates),
           testtransition(pincrease,nstates, pcrash), testreward(), gamma, maxcons)
    {
    }
};

template <typename Crit>
accuracytestresult testhomo(const testhomoproblem<Crit>& problem,
                            double accuracy, orpp::index s0ind,
                            unsigned testiters,
                            const typename testhomoproblem<Crit>::computationparams& cp
                            )
{
    finitevaluefunction initV(problem,0);
    auto res = problem.valueiteration(initV,accuracy,cp);

    return testevaluatehomo(problem, s0ind, res.p ,accuracy, cp, testiters);
}

template <typename Crit>
void testhomotime(const testhomoproblem<Crit>& problem,
                            double accuracy, orpp::index s0ind,
                            unsigned testiters,
                            const typename testhomoproblem<Crit>::computationparams& cp
                            )
{
    finitevaluefunction initV(problem,0);
    for(unsigned i=0; i<10; i++)
       problem.valueiteration(initV,accuracy,cp);
}

template <typename Crit>
void test(unsigned nstates, unsigned maxcons,
            double kappa, double pincrease, double gamma,
          double pcrash,
             orpp::index s0ind, double accuracy,
             unsigned testiters,
             const typename testproblem<Crit>::computationparams& params)
{
   testproblem<Crit> problem(nstates, maxcons, kappa,pincrease,gamma,pcrash);

   typename testproblem<Crit>::heuristicresult res = problem.heuristic(s0ind,accuracy,params);
   testoverall(problem,res.p[s0ind],{res.p},s0ind,accuracy,testiters,params);

   testhomoproblem<Crit> hp(nstates, maxcons, res.iota, pincrease,gamma, pcrash);
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
    double kappa;
    double gamma;
    double accuracy;
    double evalaccuracy;
    orpp::index s0ind = 1;
    double eastep = 0.2;
//    unsigned testiters = 10;
    bool heuristic = false;
    bool taylorheuristic = false;
    bool heuristicplus = false;
    bool pseudogradienthomo = false;
    bool riskneutral = false;
//    bool pesudogradient = false;
    bool pseudogradienthetero = false;
    bool enumerate = false;
    unsigned fmaxstatestoenum;

};

template <typename Crit>
struct testexamineprogram: public examineprogram
{
    double pincrease = 0.7;
    unsigned nstates;
    unsigned maxcons;
    double pcrash;
    typename testproblem<Crit>::computationparams pars;
};

template <typename Crit>
struct invexamineprogram: public examineprogram
{
    double pincrease = 0.7;
    unsigned maxinv;
    unsigned lot;
    double pcrash;
    typename invproblem<Crit>::computationparams pars;
};


template <bool enumeratex, typename P, typename HP, typename R>
void examineproblem(P& problem, HP& hp, const R& p, std::ostream& report,
                    bool donotenum, finitepolicy& lastep,
                    double& lastex
                    )
{
    timems tverystart;
    timems tstart = tverystart = sys::gettimems();

    finitepolicy probablybest(0U);
    if(p.riskneutral)
    {
         sys::logline() << "riskneutral" << std::endl;
         auto vires = problem.riskaversesolution(p.accuracy,p.pars);
         timems tend = sys::gettimems();
         auto res = problem.evaluatecrit(p.s0ind,vires.p,p.evalaccuracy,p.pars);
         sys::logline() << vires.p << ": " << vires.v[p.s0ind] << " crit= " << res.x << " (" << res.sd << ")," << std::endl;
         report << vires.p << "," << vires.v[p.s0ind] << "," << res.x << "," << tend - tstart << ",";
         probablybest = vires.p;
    }
    else
        report << ",,,,";

    tstart = sys::gettimems();

    if(p.heuristic)
    {
        try
        {
            typename P::heuristicresult hres = problem.template heuristic<true>(p.s0ind,p.accuracy,p.pars);
            timems tend = sys::gettimems();
            auto res = problem.evaluatecrit(p.s0ind,hres.p,p.evalaccuracy,p.pars);

            report << hres.p << "," << res.x << "," << hres.iota << "," << tend - tstart << ",";
        }
        catch(const timelimitexception& e)
        {
            report << ",outoftime,,,";
        }

        tstart = sys::gettimems();

        try
        {
            typename P::heuristicresult hres = problem.heuristic(p.s0ind,p.accuracy,p.pars);
            timems tend = sys::gettimems();
            auto res = problem.evaluatecrit(p.s0ind,hres.p,p.evalaccuracy,p.pars);

            report << hres.p << "," << res.x << "," << hres.iota << "," << tend - tstart << ",";
            probablybest = hres.p;

        }
        catch(const timelimitexception& e)
        {
            report << ",outoftime,,,";
        }

    }
    else
        report << ",,,," << ",,,,";

    tstart = sys::gettimems();

    if(p.taylorheuristic)
    {
        try
        {
            typename P::heuristicresult hres = problem.template taylorheuristic<true>(p.s0ind,p.accuracy,p.pars);
            timems tend = sys::gettimems();

            auto res = problem.evaluatecrit(p.s0ind,hres.p,p.evalaccuracy,p.pars);

            report << hres.p << ","
                   << res.x << "," << hres.iota << ","  << tend - tstart << ",";;;
        }
        catch(const timelimitexception& e)
        {
            report << ",outoftime,,,";
        }

        tstart = sys::gettimems();


        try
        {
            typename P::heuristicresult hres = problem.taylorheuristic(p.s0ind,p.accuracy,p.pars);
            timems tend = sys::gettimems();

            auto res = problem.evaluatecrit(p.s0ind,hres.p,p.evalaccuracy,p.pars);

            report << hres.p << ","
                   << res.x << "," << hres.iota << ","  << tend - tstart << ",";;;
            probablybest = hres.p;

        }
        catch(const timelimitexception& e)
        {
            report << ",outoftime,,,";
        }

    }
    else
        report << ",,,," << ",,,,";

    tstart = sys::gettimems();

//
    if(p.heuristicplus)
    {
        try
        {
            typename P::heuristicplusresult hpres =
            problem.heuristicplus(p.s0ind,p.accuracy,p.pars);
            timems tend = sys::gettimems();

            auto res = problem.evaluatecrit(p.s0ind,hpres.pgres.p ,p.evalaccuracy,p.pars);

            report << hpres.hres.p << ","
                   << hpres.hres.v << "," << hpres.hres.iota << ","
                   << hpres.pgres.p << "," << res.x << ","  << tend - tstart << ",";;;
            probablybest = hpres.pgres.p;
        }
        catch(const timelimitexception& e)
        {
            report << ",,,,outoftime,,";
        }

    }
    else
        report << ",,,,,,";

    tstart = sys::gettimems();

    if(p.pseudogradienthomo)
    {
        try
        {
            auto vires = problem.riskaversesolution(p.accuracy,p.pars);
            typename P::pgdhomoresult respg = problem.pseudogradientdescenthomo(p.s0ind, vires.p, p.accuracy, p.pars);
            timems tend = sys::gettimems();

            auto res = problem.evaluatecrit(p.s0ind,respg.p,p.evalaccuracy,p.pars);

            report << respg.p;
            report << "," <<res.x << ","  << tend - tstart << ",";;;
        }
        catch(const timelimitexception& e)
        {
            report << ",outoftime,,";
        }

    }
    else
        report << ",,,";

    tstart = sys::gettimems();

    if(p.enumerate) // tbd still only to log
    {
        if(donotenum)
        {
            if(lastep.size() == 0)
                report << ",outoftime,,";
            else
                report << lastep << "," << lastex << ","  <<  ",";
        }
        else
        {
            auto numpolicies = problem.numpolicyvalues();
            sys::logline() << "# of policy values: " << numpolicies << std::endl;
            if constexpr(!enumeratex)
            {
                if(numpolicies >= p.fmaxstatestoenum)
                {
                    sys::logline() << "Too much states to enumerate" << std::endl;
                    report << ",toomuchstates,";
                }
                else try
                {
                    typename P::enumresult res = problem.enumeratehomo(p.s0ind, p.evalaccuracy, p.pars);
                    timems tend = sys::gettimems();

                    sys::log() << "best of enumerate:" << res.p << " "
                               << res.v.x  << std::endl;
                    report << res.p << "," << res.v.x << ","  << tend - tstart << ",";

                    lastep = res.p;
                    lastex = res.v.x;
                }
                catch(const timelimitexception& e)
                {
                    lastep = finitepolicy(0U);
                    report << ",outoftime,,";
                }
            }
            else
            {
                try
                {
                    typename P::enumresult res = problem.enumeratehomoex(p.s0ind, p.evalaccuracy, p.pars);
                    timems tend = sys::gettimems();


                    sys::log() << "best of enumerate:" << res.p << " "
                               << res.v.x  << std::endl;
                    report << res.p << "," << res.v.x << "," << tend - tstart << ",";
                    lastep = res.p;
                    lastex = res.v.x;
                }
                catch(const timelimitexception& e)
                {
                    lastep = finitepolicy(0U);
                    report << ",outoftime,,";
                }
            }
        }
    }
    else
        report << ",,,";

    tstart = sys::gettimems();

    if(p.pseudogradienthetero)
    {
        if(probablybest.size() == 0)
            report << ",nostartpset,,";
        else try
        {
            typename P::heteropolicy heterop = { probablybest[p.s0ind], {probablybest, probablybest} };

            typename P::pgdheteroresult respg = problem.pseudogradientdescent(p.s0ind, heterop, p.accuracy, p.pars);
            timems tend = sys::gettimems();

            auto res = problem.evaluatecrit(p.s0ind,respg.p,p.evalaccuracy,p.pars);

            report << respg.p.p0;
            for(unsigned k=0;k < respg.p.ps.size(); k++)
            {
                report << "-";
                report << respg.p.ps[k];
            }
            report << "," <<respg.v.x << "," << tend - tstart << ",";
        }
        catch(const timelimitexception& e)
        {
            report << ",outoftime,sys::gettimems()-tstart,";
        }

    }
    else
        report << ",,,";



    // Calculating total time taken by the program.
    double time_taken = (sys::gettimems() - tverystart) / 1000.0;
    sys::logline() << "Time taken by program is : " << std::fixed
        << time_taken << std::setprecision(5);

    report << time_taken << std::endl;

    sys::log() << " sec " << std::endl;

}

enum eanalysis {egrid, elarge, eapprox, enumanalyses};

template <typename P, typename R, typename C>
void domain(unsigned nthreads, std::string repontname, eanalysis e)
{
    sys::setlog(std::cout);
    sys::logline() << "Using " << nthreads << " threads." << std::endl;

    sys::setloglevel(0);

    R p;



    typename P::computationparams pars;

    p.s0ind = 1;

    p.eastep = 0.2;


    pars.fthreadstouse = pars.fnestedtaylorparams.fthreadstouse
             = pars.fnestedparams.fthreadstouse = nthreads;
    pars.fthreadbatch = pars.fnestedtaylorparams.fthreadbatch
             = pars.fnestedparams.fthreadbatch = 5000;

    p.evalaccuracy = 0.0025;
    p.fmaxstatestoenum = 10000;
    p.pars = pars;

    p.riskneutral = true;
    p.pseudogradienthomo = true;
if(0)
{
    p.heuristicplus = true;
    p.heuristic = true;
    p.taylorheuristic = true;
    p.heuristic = true;
}
    if(e==egrid)
    {
        p.enumerate = true;
        p.pseudogradienthetero = true;
    }


    std::ofstream report(repontname);
    if(!report)
    {
        throw exception("cannot open rep");
    }

    report << "problem,crit,";
    report << "nstates/maxinv,maxcons/lot,kappa,gamma,pcrash,accuracy,evalaccuracy,";
    std::string id;
    if constexpr(std::is_same<testexamineprogram<C>,R>::value)
    {
            id = "cons";
    }

    if constexpr(std::is_same<invexamineprogram<C>,R>::value)
    {
            id = "inv";
    }

    std::string cid;
    if constexpr(std::is_same<critcvar,C>::value)
            cid = "cvar";
    if constexpr(std::is_same<critmcv,C>::value)
            cid = "MCV";


    if(e == eapprox)
    {
       report << "ref,";
       for(double q = 0; q < 1.01; q+=p.eastep)
           report << "a" << q << ",";
       for(double q = 0; q < 1.01; q+=p.eastep)
           report << "e" << q << ",";
       for(double q = 0; q < 1.01; q+=p.eastep)
            report << "p" << q << ",";
       report << std::endl;
    }
    else
       report << "rnpolicy,rnexp,rncrit,rntime,"
          << "qhpolicy,qhcrit,qhlambda,qhtime,"
          << "hpolicy,hcrit,hlambda,htime,"
          << "qtpolicy,qtcrit,qtlambda,qttime,"
          << "tpolicy,tcrit,tlambda,ttime,"
          << "hppolicyh,hpcrith,hplambdah,hppolicyp,hpcritp,hptime,"
          << "pghomopolicy,pghomocrit,pghomotime,"
          << "enumgpolicy,enumgcrit,enumegtime,"
          << "pgheteropolicy,pgheterocrit,pgheterotime,"
          << "timetaken"
          << std::endl;
    report << std::setprecision(5);

    double lastex = 0;
    finitepolicy lastep(0U);
    bool donotenum = false;


    if(e==egrid)
        for(double kappa = 0.1; kappa < 0.91; kappa += 0.2)
            for(double pcrash = 0.025; pcrash < 0.126; pcrash += 0.025)
            {
                pars.fopttimelimit = pars.fpseudogradienttimelimit
                        = pars.fenumtimelimit = 7200000;


                donotenum = false;
                for(double accuracy = 0.01; accuracy < 0.2; accuracy *= 4)
        //for(double kappa = 0.7; kappa < 0.71; kappa += 0.2)
        //for(double pcrash = 0.125; pcrash < 0.126; pcrash += 0.025)
                {
                    p.accuracy = accuracy;

                    p.pcrash = pcrash;
                    p.gamma = 0.8;
                    p.kappa = kappa;

                    report << id << "," << cid << ",";
                    if constexpr(std::is_same<testexamineprogram<C>,R>::value)
                    {
                        p.nstates = 10;
                        p.maxcons = 5;
                        p.pincrease = 0.7;
                        report << p.nstates << "," << p.maxcons << "," << p.kappa << "," << p.gamma << ","
                               << p.pcrash << "," << p.accuracy << "," << p.evalaccuracy << ",";

                        sys::logline() << "kappa, gamma, pcrash = "
                                       << p.kappa << ", " << p.gamma << ", "
                                       << p.pcrash << std::endl;
                        testproblem<C> problem(p.nstates, p.maxcons, p.kappa, p.pincrease, p.gamma, p.pcrash);
                        testhomoproblem<C> hp(p.nstates, p.maxcons, 0, p.pincrease,p.gamma, p.pcrash);

                        examineproblem<false>(problem,hp,p,report,donotenum,lastep,lastex);
                        sys::logline() << std::endl;

                    }
                    if constexpr(std::is_same<invexamineprogram<C>,R>::value)
                    {
                        p.maxinv = 4;
                        p.lot = 2;
                        p.pincrease = 0.7;
                        report << p.maxinv << "," << p.lot << "," << p.kappa << "," << p.gamma << ","
                               << p.pcrash << "," << p.accuracy << "," << p.evalaccuracy << ",";

                        sys::logline() << "kappa, gamma, pcrash = "
                                       << p.kappa << ", " << p.gamma << ", "
                                       << p.pcrash << std::endl;
                        invproblem<C> problem(p.maxinv,p.lot, p.kappa, p.pincrease, p.gamma, p.pcrash);
                        invhomoproblem<C> hp(p.maxinv,p.lot, p.kappa, p.pincrease, p.gamma, p.pcrash);

                        examineproblem<true>(problem,hp,p,report,donotenum,lastep,lastex);
                        sys::logline() << std::endl;

                    }
                    donotenum = true;
                }
            }

    if(e==elarge)
    {
        p.accuracy = 0.01;
        p.pcrash = 0.033;
        p.gamma = 0.8;
        p.kappa = 0.7;


        for(unsigned i=
            //4
8; i<=10; i+=2)
        {
            report << id << "," << cid << ",";
            pars.fopttimelimit = pars.fpseudogradienttimelimit
                    = pars.fenumtimelimit = 3*7200*1000;

            if constexpr(std::is_same<testexamineprogram<C>,R>::value)
            {
                p.nstates = i;
                p.maxcons = i / 2;
                p.pincrease = 0.7;
                report << p.nstates << "," << p.maxcons << "," << p.kappa << "," << p.gamma << ","
                       << p.pcrash << "," << p.accuracy << "," << p.evalaccuracy << ",";;

                sys::logline() << "kappa, gamma, pcrash = "
                               << p.kappa << ", " << p.gamma << ", "
                               << p.pcrash << std::endl;
                testproblem<C> problem(p.nstates, p.maxcons, p.kappa, p.pincrease, p.gamma, p.pcrash);
                testhomoproblem<C> hp(p.nstates, p.maxcons, 0, p.pincrease,p.gamma, p.pcrash);

                examineproblem<false>(problem,hp,p,report,donotenum,lastep,lastex);
                sys::logline() << std::endl;

            }
            if constexpr(std::is_same<invexamineprogram<C>,R>::value)
            {
                p.maxinv = i;
                p.lot = 2;
                p.pincrease = 0.7;
                report << p.maxinv << "," << p.lot << "," << p.kappa << "," << p.gamma << ","
                       << p.pcrash << "," << p.accuracy << "," << p.evalaccuracy << ",";

                sys::logline() << "kappa, gamma, pcrash = "
                               << p.kappa << ", " << p.gamma << ", "
                               << p.pcrash << std::endl;
                invproblem<C> problem(p.maxinv,p.lot, p.kappa, p.pincrease, p.gamma, p.pcrash);
                invhomoproblem<C> hp(p.maxinv,p.lot, p.kappa, p.pincrease, p.gamma, p.pcrash);

                examineproblem<true>(problem,hp,p,report,donotenum,lastep,lastex);
                sys::logline() << std::endl;

            }
        }
    }
    if(e==eapprox)
    {
        if constexpr(std::is_same<invexamineprogram<C>,R>::value)
        {
            p.evalaccuracy = 0.01;
            std::vector<unsigned> pb = { 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 6,6,6,6,6,6,6,6,6,6,6 };

            p.maxinv = 12;
            double q = 0.5;
            std::vector<unsigned> r;
            for(unsigned i=0; i<pb.size(); i++)
            {
                double x = (i % 2) ? 0.51 : 0.49;
                r.push_back( static_cast<unsigned>(q*pb[i] + x));
            }
            finitepolicy pr(r);
            sys::logline() << "referemce:" << pr << std::endl;


            for(double kappa = 0.1; kappa < 0.91; kappa += 0.2)
                for(double pcrash = 0.025; pcrash < 0.126; pcrash += 0.025)
                {
                    p.pcrash = pcrash;
                    p.gamma = 0.8;
                    p.kappa = kappa;
                    p.lot = 2;

                    invproblem<C> problem(p.maxinv,p.lot, p.kappa, p.pincrease, p.gamma, p.pcrash);

                    auto res = problem.evaluatecrit(p.s0ind,pr,p.evalaccuracy,p.pars);


                    double iota = problem.findiota(pr,res.x,finitevaluefunction(problem,0),
                                     p.s0ind,
                                     p.evalaccuracy,
                                     p.pars.fnestedparams
                                     );

                    invhomoproblem<C> hp(p.maxinv,p.lot, iota, p.pincrease, p.gamma, p.pcrash);

                    sys::logline() << "kappa=" << kappa
                                      << " pcrash=" << pcrash
                                      << " iota=" << iota
                                      << std::endl ;
                    report << id << "," << cid << ",";
                    report << p.maxinv << "," << p.lot << "," << p.kappa << "," << p.gamma << ","
                           << p.pcrash << ",," << p.evalaccuracy << ","
                            << pr << ",";

                    std::vector<finitepolicy> ps;
                    std::vector<double> es;
                    for(double q=0; q<=1.0001; q+=p.eastep)
                    {
                        std::vector<unsigned> r;
                        for(unsigned i=0; i<pb.size(); i++)
                        {
                            double x = (i % 2) ? 0.51 : 0.49;
                            r.push_back( static_cast<unsigned>(q*pb[i] + x));
                        }
                        finitepolicy cp(r);
                        auto hres = hp.evaluate(finitevaluefunction(problem,0),cp,
                                       p.evalaccuracy,p.pars.fnestedparams);
                        double av = hres.x[p.s0ind];
                        auto res = problem.evaluatecrit(p.s0ind,cp,p.evalaccuracy,p.pars);
                        double ev = res.x;

                        sys::logline() << "policy:" << cp
                                       << " (" << ev << "-" << av << ")="
                                       << ev-av << std::endl;
                        report << av << ",";
                        es.push_back(ev);
                        ps.push_back(cp);
                    }
                    for(unsigned k=0; k<es.size(); k++)
                        report << es[k] << ",";
                    for(unsigned k=0; k<ps.size(); k++)
                        report << ps[k] << ",";
                    report << std::endl;
                }
        }
        else
            throw exception("approx not implemented for cons");
    }

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

    eanalysis e = egrid;
    if(argc > 2)
        switch(argv[2][0])
        {
        case 'G':
            e = egrid;
            break;
        case 'L':
            e = elarge;
            break;
        case 'A':
            e = eapprox;
            break;
        default:
            throw exception("Unknown option in the second argument");
        }

    bool test = false;
    if(argc>3)
    {
        if(argv[3][0] == 'I')
            test = false;
        if(argv[3][0] == 'C')
            test = true;
    }

    bool cvar = true;
    if(argc>4)
    {
        if(argv[4][0] == 'M')
            cvar = false;
        if(argv[4][0] == 'C')
            cvar = true;
    }

    std::string rn;
    if(e == egrid)
        rn = "grid";
    else if(e == elarge)
        rn = "large";
    else if(e == eapprox)
        rn = "approx";

    if(test)
    {
        if(cvar)
        {
            using crit = critcvar;
            domain<testproblem<crit>,testexamineprogram<crit>,crit>(nthreads, rn + "CC.csv",e);
        }
        else
        {
            using mcrit = critmcv;
            domain<testproblem<mcrit>,testexamineprogram<mcrit>,mcrit>(nthreads, rn +"CM.csv",e);
        }
    }
    else
    {
        if(cvar)
        {
            using crit = critcvar;
            domain<invproblem<crit>,invexamineprogram<crit>,crit>(nthreads, rn + "IC.csv",e );
        }
        else
        {
            using mcrit = critmcv;
            domain<invproblem<mcrit>,invexamineprogram<mcrit>,mcrit>(nthreads, rn + "IM.csv",e);
        }
    }
    return 0;
}




