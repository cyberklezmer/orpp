// sellproblem.hpp

#ifndef SELLPROBLEM_HPP
#define SELLPROBLEM_HPP


#include "radp.hpp"
namespace orpp 
{

class sellstatespace : public integerspace
{
public:
    sellstatespace(unsigned endowment, unsigned nprices) :
        integerspace(0,(endowment+1)*nprices-1),
        fendowment(endowment), fnprices(nprices) {}
    static unsigned price(unsigned s, unsigned nprices)  { return s % nprices; }
    static unsigned inventory(unsigned s, unsigned nprices) { return s / nprices; }
    static unsigned state(unsigned inventory, unsigned price, unsigned nprices)
      { return inventory * nprices + price; }
private:
    unsigned fendowment;
    unsigned fnprices;
};

class sellactionspace :
                        public integerspace, public constrainedspace<unsigned int,unsigned int>
{
public:
    sellactionspace(unsigned endowment, unsigned nprices) :
        integerspace(0,endowment), fendowment(endowment), fnprices(nprices)
    {
    }
    virtual bool isfeasible(const unsigned int& a, const unsigned int& s) const
    {
        auto e = sellstatespace::inventory(s,fnprices);
        if(e == 0)
            return a == 0;
        else
            return a > 0 && a <= e;
    }
private:
    unsigned fendowment;
    unsigned fnprices;
};

class sellreward : public dpreward<sellstatespace, sellactionspace>
{
public:
    sellreward(unsigned nprices) : fnprices(nprices) {}
    double operator() (const dpcondition<unsigned int, unsigned int>& x) const
    {
        return x.a * sellstatespace::price(x.s,fnprices);
    }
private:
    unsigned fnprices;
};

class selltransition: public finitetransition<sellstatespace,sellactionspace>
{
public:
    selltransition(unsigned nprices, probability fp, unsigned endowment)
        : fnprices(nprices), fendowment(endowment), fD(nprices-1,fp)
    {}
private:
    virtual unsigned natoms_is(const dpcondition<unsigned int,unsigned int>&) const
    { return fnprices * (fendowment+1); }
    virtual atom<unsigned int> atom_is(unsigned int i, const dpcondition<unsigned int,unsigned int>& c) const
    {
        unsigned inv = sellstatespace::inventory(c.s,fnprices);
        assert(c.a <= inv);
        unsigned newinv = inv-c.a;
        unsigned iinv = sellstatespace::inventory(i,fnprices);
        if(iinv == newinv)
        {
            unsigned iprice = sellstatespace::price(i,fnprices);
            probability p =  fD(iprice).p;
            return {i,p};
        }
        else
            return {i,0};
    }
    virtual bool is_sorted() const { return true; }
    unsigned fnprices;
    unsigned fendowment;
    binomialdistribution fD;
};

template <typename Crit>
class sellproblem : public overallriskproblem<Crit,
     sellstatespace, sellactionspace, selltransition, sellreward>
{
public:
    sellproblem(unsigned endowment, unsigned nprices, probability param,
                probability alpha, double gamma) :
        overallriskproblem<Crit, sellstatespace,
         sellactionspace, selltransition, sellreward>
           (Crit(alpha),sellstatespace(endowment,nprices),
               sellactionspace(endowment, nprices),
               selltransition(nprices, param, endowment),
               sellreward(nprices), gamma, nprices-1,2.0 / (1-alpha)) {}


};

template<typename A, int NA = -1>
class fullstrategy
{
private:
    static std::vector<std::vector<A>> makedata(unsigned T, unsigned k)
    {
        std::vector<std::vector<A>> result(T+1);
        unsigned num = 1;
        for(unsigned i=0; i<=T; i++)
        {
            result[i] = std::vector<A>(num,NA);
            num *= k;
        }
        return result;
    }
public:
    fullstrategy(unsigned T, unsigned k): fT(T), fk(k), fData(makedata(T,k))
    {
    }

    A& operator[](const std::vector<unsigned>& args)
    {
        auto len = args.size();
        assert(len <= fT+1);
        unsigned multiplier = 1;
        unsigned index = 0;
        for(unsigned i=0;; i++)
        {
            assert(args[i] < fk);
            index += multiplier * args[len-i-1];
            if(i+1== len)
            {
                assert(index < fData[i].size());
                return fData[i][index];
            }
            multiplier *= fk;
        }
    }
    std::vector<A>& level(unsigned i)
    {
        assert(i<=fT);
        return fData[i];
    }
    // put this in your class:
    void output(std::ostream& o) const
    {
        // path[i] is the child index at depth i (0..fk-1)
        std::vector<unsigned> path;
        // DFS that prints value at each node, then recurses to children
        std::function<void(unsigned)> dfs = [&](unsigned depth)
        {
            // compute linear index in fData[depth] from the current path (base fk)
            size_t idx = 0;
            for (unsigned i = 0; i < depth; ++i)
                idx = idx * fk + path[i];

            // indent by depth and print the stored value
            for (unsigned i = 0; i < depth; ++i) o << "  ";
            o << fData[depth][idx] << '\n';

            // stop at leaves
            if (depth == fT) return;

            // recurse over children
            for (unsigned a = 0; a < fk; ++a) {
                path.push_back(a);
                dfs(depth + 1);
                path.pop_back();
            }
        };

        dfs(0); // start from root (depth 0, path empty ⇒ idx = 0)
    }
    unsigned T() const { return fT; }
    unsigned k() const { return fk; }
private:

    std::vector<std::vector<A>> fData;
    unsigned fT;
    unsigned fk;
};



template <class F>
void enumerate_under(const std::vector<int>& s, F&& visit) {
    const size_t n = s.size();
    if (n == 0) { visit(std::vector<int>{}); return; }

    std::vector<int> x(n, 0);
    // Force zeros where s[i]==0 (others start at 0 anyway).
    for (size_t i = 0; i < n; ++i) if (s[i] == 0) x[i] = 0;

    // Edge case: if any s[i] < 0, there are no valid vectors (adjust if needed).
    for (int v : s) if (v < 0) return;

    // Count of combinations could be huge; iterate lexicographically.
    while (true) {
        visit(x);

        // increment like a mixed-radix counter
        size_t i = 0;
        for (; i < n; ++i) {
            if (s[i] == 0) continue;            // fixed 0 digit
            if (x[i] + 1 < s[i]) {              // can increment this digit
                ++x[i];
                // reset all lower "digits" that are not fixed to 0
                for (size_t j = 0; j < i; ++j) if (s[j] > 0) x[j] = 0;
                break;
            }
        }
        if (i == n) break; // overflowed all positions: done
    }
}

void createpolicies( std::vector<fullstrategy<int>>& list,  unsigned k /* nprices */, unsigned n /*maxinv*/ )
{
    list.clear();

    auto T = n-2;
    std::vector<fullstrategy<int>> src;

    fullstrategy<int> ns(T,k);
    for(unsigned i=0; i<n; i++)
    {
        std::vector<unsigned> index(1,0);
        ns[index] = i;
        src.push_back(ns);
    }
    std::vector<fullstrategy<int>> dst;
    for(unsigned j=1; j<=T; j++)
    {
        for(unsigned i=0; i<src.size(); i++)
        {
            fullstrategy<int> cs = src[i];
            std::vector<int> bounds = cs.level(j-1);
            std::vector<int> boundswide;
            for(unsigned m=0; m<bounds.size(); m++)
                for(unsigned n=0; n<k; n++)
                    boundswide.push_back(bounds[m]);
            enumerate_under(boundswide, [cs,j,&dst](const std::vector<int>& x){
                fullstrategy<int> nd = cs;
                nd.level(j) = x;
                dst.push_back(nd);
            });
        }
        src = dst;
        dst.clear();
    }
    list = src;
}


}
#endif

