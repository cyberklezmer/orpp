#include <vector>
#include <iostream>
#include <assert.h>
#include <cassert>
#include "orpp/boostdist.hpp"
#include "radp.hpp"
#include "invproblem.hpp"
#include "sellproblem.hpp"
#include "testproblem.hpp"


enum eanalysis {eexp1, eexp2, eapprox, enumanalyses};

using namespace orpp;


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
inline void examineproblem(P& problem, HP& hp, const R& p, std::ostream& report,
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


template <typename P, typename R, typename C>
inline void domain(unsigned nthreads, std::string repontname, eanalysis e)
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

    if(e==eexp1)
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


    if(e==eexp1)
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

    if(e==eexp2)
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

void generatehomostrats(int pos, unsigned endowment,
              std::vector<unsigned>& v,
              std::vector<std::vector<unsigned>>& results) {
    if (pos > endowment) {
        results.push_back(v);  // save a copy of the vector
        return;
    }

    unsigned llimit, hlimit;
    if(pos == 0)
    {
        llimit = hlimit = 0;
    }
    else
    {
        llimit = 1;
        hlimit = pos;
    }

    for (unsigned val = llimit; val <= hlimit; ++val) {
        v[pos] = val;
        generatehomostrats(pos + 1, endowment, v, results);
    }
}


ldistribution<double, true> overalldistrecourse( fullstrategy<int>& strat,
                                    unsigned inventory,
                                  probability bparam, double gamma,
                                   std::vector<unsigned> index)
{
    unsigned level = index.size();
    bool laststage = level == strat.T();
    binomialdistribution b(strat.k()-1,bparam);

    std::vector<probability> weights(strat.k());
    std::vector<ldistribution<double,true>> dists;
    //    std::vector<scaledfdistribution<binomialdistribution>> bdists;
    for(unsigned k=0; k < strat.k(); k++)
    {
        weights[k] = b(k).p;
        auto newindex = index;
        newindex.push_back(k);
        unsigned inv = strat[index];
        unsigned newinv = strat[newindex];
        assert(inv >= newinv);
        unsigned selling = inv - newinv;
        scaledfdistribution<binomialdistribution> thisstage(b,selling);

        if(laststage)
        {
            if(newinv > 0)
            {
                scaledfdistribution<binomialdistribution> sb(b,gamma * newinv);
                ldistribution<double,true> cv = convolution(thisstage,sb);
                dists.push_back(cv);
            }
            else
            {
                diracdistribution<double> dd(0.0);
                // this is done because push_back wont accept different type
                ldistribution<double,true> cv = convolution(thisstage,dd);
                dists.push_back(cv);
            }
        }
        else
        {
           auto d =  overalldistrecourse(strat,
                                      newinv,
                                      bparam, gamma,
                                      newindex);
           scaledfdistribution<ldistribution<double,true>> sd(d,gamma);
           auto topush = convolution(sd,thisstage);
           dists.push_back(topush);
        }
    }
    auto result = mixture(dists,weights);
    return result;
}



ldistribution<double, true> overalldist( fullstrategy<int>& strat,
                                        unsigned endowment,
                                        unsigned initprice,
                                        probability bparam, double gamma)
{
    std::vector<unsigned> index(1,0);
    int newinv = strat[index];
    unsigned selling = endowment - newinv;
    ldistribution<double,true> first({static_cast<double>(selling*initprice)},true);
    auto next = overalldistrecourse(strat,newinv, bparam, gamma, index);
    scaledfdistribution<ldistribution<double,true>> snext(next,gamma);
    auto result = convolution(first,snext);
    return result;
}

double findbeststragegy(const std::vector<fullstrategy<int>>&  list,
                                  const unsigned& endowment,
                                  const unsigned initprice,
                                  const probability bparam,
                                  const double gamma,
                                  double kappa,
                                  fullstrategy<int>& bestp)
{
    bool firsttime = true;
    double bestcrit;

    for (auto s : list)
    {
        ldistribution<double, true> dist = overalldist(s, endowment, initprice, bparam, gamma);

        CVaR<ldistribution<double,true>,true> crit(kappa);
        double c = crit(dist);
        bool star = false;
        if( firsttime || c > bestcrit )
        {
            bestcrit = c;
            bestp = s;
            star = true;
        }
        firsttime = false;
        if (sys::loglevel() >= 2)
        {
            sys::log() << "Strategy " << std::endl;
            if (sys::loglevel() >= 2)
                s.output(sys::log());

            if (s.T() == 1)
            {
                sys::log() << "Check: " << s[{0}] << " - ";
                for (unsigned i = 0; i < s.k(); ++i)
                    sys::log() << s[{0, i}] << " ";
                sys::log() << std::endl;
            }

            if (s.T() == 2)
            {
                sys::log() << "Check: " << s[{0}] << " - ";
                for (unsigned i = 0; i < s.k(); ++i)
                {
                    sys::log() << s[{0, i}] << "(";
                    for (unsigned j = 0; j < s.k(); ++j)
                        sys::log() << s[{0, i, j}] << " ";
                    sys::log() << ") ";
                }
                sys::log() << std::endl;
            }

            if (sys::loglevel() >= 2)
            {
                sys::log() << "Distribution " << std::endl;

                std::vector<atom<double>> atoms;
                dist.atoms(atoms);

                sys::out() << std::fixed << std::setprecision(6);
                for (const auto& a : atoms)
                {
                    sys::out() << "x = " << a.x << ", p = " << a.p << "\n";
                }
                sys::log().flush();
            }
            sys::logline() << "crit = " << c;
            if(star)
                sys::log() << "*";
            sys::log() << std::endl;
        }

    }
    return bestcrit;
}

template <bool homo=true>
bool is_homo(fullstrategy<int>& strat, std::vector<unsigned> index,
             unsigned endowment,
             std::vector<std::vector<std::vector<int>>>& notebook )
{
    assert(index.size()>=1);

    if(index.size()==strat.T()+1)
    {
        if(strat[index] > 1)
            throw "Strategy too short";
    }
    unsigned k = index[index.size()-1];
    unsigned prevstate;
    if(index.size() == 1)
        prevstate = endowment;
    else
    {
        auto index2 = index;
        index2.pop_back();
        prevstate = strat[index2];
    }
    unsigned hi = homo ? 0 : index.size()-1;
    int& ntb = notebook[hi][prevstate][k];
    if(ntb == -1)
        ntb = strat[index];
    else if(ntb != strat[index])
        return false;
    if(index.size()<=strat.T())
    {
        for(unsigned i=0; i<strat.k(); i++)
        {
            auto newindex = index;
            newindex.push_back(i);
            // endowment for nothing in nested levels, but lazy to create different funciton
            if(!is_homo<homo>(strat,newindex,endowment,notebook))
                return false;
        }
    }
    return true;
}

template <bool homo>
bool is_homo(fullstrategy<int>& strat, unsigned endowment)
{
    std::vector<std::vector<std::vector<int>>>
        notebook(strat.T()+1,std::vector<std::vector<int>>(endowment+1,std::vector<int>(strat.k(),-1)));


       return is_homo<homo>(strat,{0}, endowment, notebook);
}

template <typename C>
void dosell(unsigned nthreads, std::string repontname)
{

    sys::setlog(std::cout);
    sys::logline() << "Using " << nthreads << " threads." << std::endl;

    sys::setloglevel(1);

    typename sellproblem<C>::computationparams pars;

    pars.fthreadstouse = pars.fnestedtaylorparams.fthreadstouse
        = pars.fnestedparams.fthreadstouse = nthreads;
    pars.fthreadbatch = pars.fnestedtaylorparams.fthreadbatch
        = pars.fnestedparams.fthreadbatch = 5000;

    double evalaccuracy = 0.0025
                          *10

        ;
    unsigned pfmaxstatestoenum = 10000;

    std::ofstream report(repontname);
    if(!report)
    {
        throw exception("cannot open rep");
    }

    report << "problem,crit,";
    report << "nstates/maxinv,maxcons/lot,kappa,gamma,pcrash,accuracy,evalaccuracy,";
    std::string id;

    id = "sell";

    std::string cid;
    if constexpr(std::is_same<critcvar,C>::value)
        cid = "cvar";
    if constexpr(std::is_same<critmcv,C>::value)
        cid = "MCV";

    report << std::endl;
    report << std::setprecision(5);

    double gamma = 0.948;
    double kappa = 0.7;
    probability bparam = 0.76;
    unsigned endowment = 5;
    unsigned nprices = 3;
    unsigned initprice = 1;
    {
        sellproblem<C> problem(endowment, nprices, bparam, kappa, gamma);

        std::vector<fullstrategy<int>> list;
        createpolicies(list, nprices, endowment);
        sys::log() << "Number of policies: " << list.size() << std::endl;

        std::vector<fullstrategy<int>> homo;
        homo.reserve(list.size()); // optional
        for(auto s: list)
        {
            if(is_homo<true>(s,endowment))
                homo.push_back(s);
        }
        sys::log() << "Number of homogenous policies: " << homo.size() << std::endl;

        std::vector<fullstrategy<int>> hetero;
        hetero.reserve(list.size()); // optional
        hetero.reserve(list.size()); // optional
        for(auto s: list)
        {
            if(is_homo<false>(s,endowment))
                hetero.push_back(s);
        }
        sys::log() << "Number of homogenous policies: " << hetero.size() << std::endl;


        throw;

        fullstrategy<int> bests(list[0].T(),list[0].k());
        double bestenum;
        if constexpr(std::is_same<critcvar,C>::value)
        {
                double bestenum = findbeststragegy(list,
                                   endowment,
                                   initprice,
                                   bparam,
                                   gamma,
                                   kappa,
                                   bests);
            sys::log() << "Enumeration: opt. crit=" << bestenum << std::endl;
            bests.output(sys::log());
        }
        else
        {
            throw "not implemented";
        }


        sys::log().flush();
        throw;

        auto initstate = sellstatespace::state(endowment,initprice,nprices);

        report << id << "," << cid << ",";

        pars.fenumtimelimit = 3*7200*1000;


        report << kappa << "," << gamma << ","
           << evalaccuracy << ",";

        sys::logline() << "kappa, gamma = " << kappa << "," << gamma << std::endl;
//    invhomoproblem<C> hp(p.maxinv,p.lot, p.kappa, p.pincrease, p.gamma, p.pcrash);

       if(sys::loglevel() >= 2)
       {
           finitepolicy tp(problem);
           sys::log() << "Testing transition given policy " << tp << std::endl;

           selltransition st(nprices, bparam, endowment);
           for(unsigned i=0; i<nprices * (endowment+1); i++)
           {
               sys::log() << "state " << i << " atoms " << st.natoms({tp[i],i}) << std::endl;
               for(unsigned j=0; j < nprices * (endowment+1); j++)
               {
                   auto a = st(j,{tp[i],i});
                   sys::log() << a.p << " ";
               }
               sys::log() << std::endl;
           }
       }


        timems tverystart;
        timems tstart = tverystart = sys::gettimems();

        unsigned endowment = 3; // example
        std::vector<unsigned> v(endowment + 1, 0);
        std::vector<std::vector<unsigned>> results;

        sys::logline() << "Solution by nested problem " << std::endl;
        generatehomostrats(0, endowment, v, results);
        std::vector<finitepolicy> pols;
        // Print results
        sys::logline() << "Generating feasible policies" << std::endl;
        for (const auto& vec : results)
        {
            finitepolicy policy(problem);
            unsigned k = 0;
            for(unsigned i=0; i<endowment+1; i++)
                for(unsigned j=0; j<nprices; j++)
                    policy[k++]=vec[i];
            if(sys::loglevel() >= 2)
            {
                sys::logline() << policy << std::endl;
            }
            pols.push_back(policy);
        }

        auto hres = problem.enumheuristic(initstate, evalaccuracy,
                      pols,
                      pars );

/*        auto res = problem.evaluatecrit(initstate,vires.p,evalaccuracy,pars);
        sys::logline() << vires.p << ": " << vires.v[initstate] << " crit= " << res.x << " (" << res.sd << ")," << std::endl;

    timems tend = sys::gettimems();
    report << vires.p << "," << vires.v[initstate] << "," << res.x << "," << tend - tstart << ",";
    probablybest = vires.p;

    //
    try
    {
        typename sellproblem<C>::heuristicplusresult hpres =
            problem.heuristicplus(initstate,accuracy,pars);
        timems tend = sys::gettimems();

        auto res = problem.evaluatecrit(initstate,hpres.pgres.p ,evalaccuracy,pars);

        report << hpres.hres.p << ","
               << hpres.hres.v << "," << hpres.hres.iota << ","
               << hpres.pgres.p << "," << res.x << ","  << tend - tstart << ",";;;
        probablybest = hpres.pgres.p;
    }
    catch(const timelimitexception& e)
    {
        report << ",,,,outoftime,,";
    }


    tstart = sys::gettimems();

    // Calculating total time taken by the program.
    double time_taken = (sys::gettimems() - tverystart) / 1000.0;
    sys::logline() << "Time taken by program is : " << std::fixed
                   << time_taken << std::setprecision(5);

    report << time_taken << std::endl;

    sys::log() << " sec " << std::endl;

    sys::logline() << std::endl; */
    }
}


int main(int argc, char *argv[])
{
    unsigned nthreads = 1;
    using crit = critcvar;

    dosell<crit>(nthreads,std::string("testreport"));
    return 1;


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

    eanalysis e = eexp1;
    if(argc > 2)
        switch(argv[2][0])
        {
        case 'G':
            e = eexp1;
            break;
        case 'L':
            e = eexp2;
            break;
        case 'A':
            e = eapprox;
            break;
        default:
            throw exception("Unknown option in the second argument");
        }

    enum etask {etest, einv, esell};
    etask et;
    if(argc>3)
    {
        if(argv[3][0] == 'I')
            et = einv;
        if(argv[3][0] == 'C')
            et = etest;
        if(argv[3][0] == 'S')
            et = esell;

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
    if(e == eexp1)
        rn = "ex1";
    else if(e == eexp2)
        rn = "ex2";
    else if(e == eapprox)
        rn = "approx";

    switch(et)
    {
    case etest:
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
        break;
    case einv:
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
        break;
    case esell:
        if(cvar)
        {
            using crit = critcvar;
            dosell<crit>(nthreads, rn + "SC.csv");
        }
        else
        {
            using mcrit = critmcv;
            dosell<crit>(nthreads, rn + "SM.csv" );
        }
        break;

    }
    return 0;
}







