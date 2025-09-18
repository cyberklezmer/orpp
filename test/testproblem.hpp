// testproblem.hpp



#ifndef TESTPROBLEM_HPP
#define TESTPROBLEM_HPP

#include "orpp/test/testdp.hpp"


namespace orpp {


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


}
#endif

