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
    orpp::index s0ind = 1;
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
void examineproblem(P& problem, HP& hp, const R& p, std::ostream& report)
{
    finitepolicy startingp(problem,0);

    unsigned tstart = sys::gettimems();


    timems viend = tstart;

    if(p.riskneutral)
    {
         sys::logline() << "riskneutral" << std::endl;
         auto vires = problem.riskaversesolution(p.accuracy,p.pars);
         viend = sys::gettimems();
         auto res = problem.evaluatecrit(p.s0ind,vires.p,p.accuracy,p.pars);
         sys::logline() << vires.p << ": " << vires.v[p.s0ind] << " crit= " << res.x << " (" << res.sd << ")," << std::endl;
         report << vires.p << "," << vires.v[p.s0ind] << "," << res.x << ",";
         startingp = vires.p;
    }
    else
        report << ",,,";
    unsigned rnend = sys::gettimems();

    report << viend - tstart << ",";


    if(p.heuristic)
    {
        try
        {
            typename P::heuristicresult hres = problem.heuristic(p.s0ind,p.accuracy,p.pars);
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
        if(!p.riskneutral)
            throw exception("Initial policy not set (p.riskneutral==false)");
        try
        {
            typename P::heuristicresult hres = problem.taylorheuristic(p.s0ind,p.accuracy,startingp,p.pars);
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

//
    if(p.heuristicplus)
    {
        try
        {
            typename P::heuristicplusresult hpres =
            problem.heuristicplus(p.s0ind,p.accuracy,p.pars);
            report << hpres.hres.p << ","
                   << hpres.hres.v << "," << hpres.hres.iota << ","
                   << hpres.pgres.p << "," << hpres.pgres.v.x << "," ;
            startingp = hpres.pgres.p;
        }
        catch(const timelimitexception& e)
        {
            report << ",,,,outoftime,";
        }

    }
    else
        report << ",,,,,";
    unsigned plusend = sys::gettimems();
    report << plusend - tend << ",";

    if(p.pseudogradienthomo)
    {
        if(!p.riskneutral)
            throw exception("Initial policy not set (p.riskneutral==false)");
        try
        {
            typename P::pgdhomoresult respg = problem.pseudogradientdescenthomo(p.s0ind, startingp, p.accuracy, p.pars);
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
    report << cdend - plusend << ",";

    if(p.enumerate) // tbd still only to log
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
                typename P::enumresult res = problem.enumeratehomo(p.s0ind, p.accuracy, p.pars);

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
        {
            try
            {
                typename P::enumresult res = problem.enumeratehomoex(p.s0ind, p.accuracy, p.pars);

                sys::log() << "best of enumerate:" << res.p << " "
                           << res.v.x  << std::endl;
                report << res.p << "," << res.v.x << ",";
            }
            catch(const timelimitexception& e)
            {
                report << ",outoftime,";
            }
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
            typename P::heteropolicy heterop = { startingp[p.s0ind], {startingp, startingp} };

            typename P::pgdheteroresult respg = problem.pseudogradientdescent(p.s0ind, heterop, p.accuracy, p.pars);
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

enum eanalysis {egrid, elarge, enumanalyses};

template <typename P, typename R, typename C>
void domain(unsigned nthreads, std::string repontname, eanalysis e)
{
    sys::setlog(std::cout);
    sys::logline() << "Using " << nthreads << " threads." << std::endl;

    sys::setloglevel(0);

    R p;

    typename P::computationparams pars;

    p.s0ind = 1;

    pars.fopttimelimit = pars.fpseudogradienttimelimit
            = pars.fenumtimelimit = 7200000;

    pars.fthreadstouse = pars.fnestedtaylorparams.fthreadstouse = pars.fnestedonedparams.fthreadstouse
             = pars.fnestedparams.fthreadstouse = nthreads;
    pars.fthreadbatch = pars.fnestedtaylorparams.fthreadbatch = pars.fnestedonedparams.fthreadbatch
             = pars.fnestedparams.fthreadbatch = 5000;

    p.fmaxstatestoenum = 10000;
    p.pars = pars;
    
    p.riskneutral = true;
    p.pseudogradienthomo = true;
    p.heuristicplus = true;
    p.heuristic = true;
    p.taylorheuristic = true;
    p.heuristic = true;

    if(e==egrid)
    {
        p.enumerate = true;
        p.pseudogradienthetero = true;
    }

//p.accuracy = 0.0015;

    std::ofstream report(repontname);
    if(!report)
    {
        throw exception("cannot open rep");
    }

    report << "problem,crit,";
    std::string id;
    if constexpr(std::is_same<testexamineprogram<C>,R>::value)
    {
            id = "cons";
            report << "nstates,maxcons,kappa,gamma,pcrash,accuracy,";
    }

    if constexpr(std::is_same<invexamineprogram<C>,R>::value)
    {
            id = "inv";
            report << "maxinv,lot,kappa,gamma,pcrash,accuracy,";
    }

    std::string cid;
    if constexpr(std::is_same<critcvar,C>::value)
            cid = "cvar";
    if constexpr(std::is_same<critmcv,C>::value)
            cid = "MCV";


    report << "rnpolicy,rnexp,rncrit,rntime,"
          << "hpolicy,hcrit,hlambda,htime,"
          << "tpolicy,tcrit,tlambda,ttime,"
          << "hppolicyh,hpcrith,hplambdah,hppolicyp,hpcritp,hptime,"
          << "pghomopolicy,pghomocrit,pghomotime,"
          << "enumgpolicy,enumgcrit,enumegtime,"
          << "pgheteropolicy,pgheterocrit,pgheterotime,"
          << std::endl;
    report << std::setprecision(5);

    if(e==egrid)
       for(double accuracy = 0.1; accuracy > 0.02; accuracy /= 2)
        for(double kappa = 0.1; kappa < 0.91; kappa += 0.2)
            for(double pcrash = 0; pcrash < 0.126; pcrash += 0.025)
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
                           << p.pcrash << "," << p.accuracy << ",";

                    sys::logline() << "kappa, gamma, pcrash = "
                                   << p.kappa << ", " << p.gamma << ", "
                                   << p.pcrash << std::endl;
                    testproblem<C> problem(p.nstates, p.maxcons, p.kappa, p.pincrease, p.gamma, p.pcrash);
                    testhomoproblem<C> hp(p.nstates, p.maxcons, 0, p.pincrease,p.gamma, p.pcrash);

                    examineproblem<false>(problem,hp,p,report);
                    sys::logline() << std::endl;

                }
                if constexpr(std::is_same<invexamineprogram<C>,R>::value)
                {
                    p.maxinv = 4;
                    p.lot = 2;
                    p.pincrease = 0.7;
                    report << p.maxinv << "," << p.lot << "," << p.kappa << "," << p.gamma << ","
                           << p.pcrash << "," << p.accuracy << ",";

                    sys::logline() << "kappa, gamma, pcrash = "
                                   << p.kappa << ", " << p.gamma << ", "
                                   << p.pcrash << std::endl;
                    invproblem<C> problem(p.maxinv,p.lot, p.kappa, p.pincrease, p.gamma, p.pcrash);
                    invhomoproblem<C> hp(p.maxinv,p.lot, p.kappa, p.pincrease, p.gamma, p.pcrash);

                    examineproblem<true>(problem,hp,p,report);
                    sys::logline() << std::endl;

                }
            }

    if(e==elarge)
    {
        p.accuracy = 0.01;
        p.pcrash = 0.033;
        p.gamma = 0.8;
        p.kappa = 0.7;

        report << id << "," << cid << ",";
        if constexpr(std::is_same<testexamineprogram<C>,R>::value)
        {
            p.nstates = 14;
            p.maxcons = 7;
            p.pincrease = 0.7;
            report << p.nstates << "," << p.maxcons << "," << p.kappa << "," << p.gamma << ","
                   << p.pcrash << "," << p.accuracy << ",";

            sys::logline() << "kappa, gamma, pcrash = "
                           << p.kappa << ", " << p.gamma << ", "
                           << p.pcrash << std::endl;
            testproblem<C> problem(p.nstates, p.maxcons, p.kappa, p.pincrease, p.gamma, p.pcrash);
            testhomoproblem<C> hp(p.nstates, p.maxcons, 0, p.pincrease,p.gamma, p.pcrash);

            examineproblem<false>(problem,hp,p,report);
            sys::logline() << std::endl;

        }
        if constexpr(std::is_same<invexamineprogram<C>,R>::value)
        {
            p.maxinv = 15;
            p.lot = 2;
            p.pincrease = 0.7;
            report << p.maxinv << "," << p.lot << "," << p.kappa << "," << p.gamma << ","
                   << p.pcrash << "," << p.accuracy << ",";

            sys::logline() << "kappa, gamma, pcrash = "
                           << p.kappa << ", " << p.gamma << ", "
                           << p.pcrash << std::endl;
            invproblem<C> problem(p.maxinv,p.lot, p.kappa, p.pincrease, p.gamma, p.pcrash);
            invhomoproblem<C> hp(p.maxinv,p.lot, p.kappa, p.pincrease, p.gamma, p.pcrash);

            examineproblem<true>(problem,hp,p,report);
            sys::logline() << std::endl;

        }
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
    else
        rn = "large";

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

