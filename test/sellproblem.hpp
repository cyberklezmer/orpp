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
    class fullstrategy
    {
    private:
        static std::vector<std::vector<int>> makedata(unsigned T, unsigned k)
        {
            std::vector<std::vector<int>> result(T+1);
            unsigned num = 1;
            for(unsigned i=0; i<=T; i++)
            {
                result[i] = std::vector<int>(num,-1);
                num *= k;
            }
            return result;
        }
    public:
        fullstrategy(unsigned T, unsigned k): fT(T), fk(k), fData(makedata(T,k))
        {
        }

        int& operator[](const std::vector<unsigned>& args)
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
        std::vector<int>& level(unsigned i)
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

        std::vector<std::vector<int>> fData;
        unsigned fT;
        unsigned fk;
    };


    static void createpolicies( std::vector<fullstrategy>& list,  unsigned k /* nprices */, unsigned n /*maxinv*/ )
    {
        list.clear();

        auto T = n-2;
        std::vector<fullstrategy> src;

        fullstrategy ns(T,k);
        for(unsigned i=0; i<n; i++)
        {
            std::vector<unsigned> index(1,0);
            ns[index] = i;
            src.push_back(ns);
        }
        std::vector<fullstrategy> dst;
        for(unsigned j=1; j<=T; j++)
        {
            for(unsigned i=0; i<src.size(); i++)
            {
                fullstrategy cs = src[i];
                std::vector<int> bounds = cs.level(j-1);
                std::vector<int> boundswide;
                for(unsigned m=0; m<bounds.size(); m++)
                    for(unsigned n=0; n<k; n++)
                        boundswide.push_back(bounds[m]);
                enumerate_under(boundswide, [cs,j,&dst](const std::vector<int>& x){
                    fullstrategy nd = cs;
                    nd.level(j) = x;
                    dst.push_back(nd);
                });
            }
            src = dst;
            dst.clear();
        }
        list = src;
    }
    template <class F>
    static void enumerate_under(const std::vector<int>& s, F&& visit) {
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


    static ldistribution<double, true> overalldistrecourse( fullstrategy& strat,
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
            diracdistribution<double> thisstage(selling*k);

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



    static ldistribution<double, true> overalldist( fullstrategy& strat,
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

    static double findbeststragegy(const std::vector<fullstrategy>&  list,
                            const unsigned& endowment,
                            const unsigned initprice,
                            const probability bparam,
                            const double gamma,
                            double kappa,
                            const MeanCVaR<ldistribution<double,true>,true>& crit,
                            fullstrategy& bestp)
    {
        bool firsttime = true;
        double bestcrit;

        for (auto s : list)
        {
            ldistribution<double, true> dist = overalldist(s, endowment, initprice, bparam, gamma);

 //           CVaR<ldistribution<double,true>,true> crit(kappa);
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
    static bool is_homo(fullstrategy& strat, std::vector<unsigned> index,
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
    static bool is_homo(fullstrategy& strat, unsigned endowment, std::vector<std::vector<std::vector<int>>>& notebook)
    {
        notebook = std::vector<std::vector<std::vector<int>>>
            (strat.T()+1,std::vector<std::vector<int>>(endowment+1,std::vector<int>(strat.k(),-1)));

        return is_homo<homo>(strat,{0}, endowment, notebook);
    }

};


}
#endif

