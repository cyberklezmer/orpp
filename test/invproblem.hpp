// invproblem.hpp



#ifndef INVPROBLEM_HPP
#define INVPROBLEM_HPP

#include "orpp/overallriskdp.hpp"
#include "orpp/test/testdp.hpp"
#include "radp.hpp"

namespace orpp {

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

}
#endif

