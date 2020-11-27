#ifndef OLDORPP_HPP_
#define OLDORPP_HPP_

class object
{
public:
    enum elogginglevel { nologging, basiclogging, extendedlogging };
    elogginglevel logging;
protected:
    ofstream nulstream;
    ostream& lout;
protected:
    object() : logging(nologging), lout(cout)
    {}
    virtual ~object() {}
public:
    void setlogging(elogginglevel al = basiclogging)
    {
        logging = al;
    }
    bool islogging() { return logging == nologging ? false : true; }
 };

/// work in progress, neodszkousene

class summary: public object
{
    vector<unsigned int> hist;
protected:
    virtual unsigned int getN() const = 0;
    virtual double X(unsigned int i) const = 0;
public:
    const vector<unsigned int>& histogram(double binw,
                    double lower, double upper)
    {
        unsigned int nbins = (upper-lower) / binw + 0.5 + 2;
        hist.resize(nbins);
        for(unsigned int i=0; i<nbins; i++)
            hist[i]=0;
        unsigned int n=getN();
        for(unsigned int i=0; i<n; i++)
        {
            int j=(X(i)-lower)/binw+1;
            if(j<0)
                j = 0;
            else if((unsigned)j>= nbins)
                j=nbins - 1;
            hist[j]++;
        }
        return hist;
    }
};

class estimator: public virtual object
{
    vector<paraminfo> params;
protected:
    const vector<paraminfo>& getparams()
    {
        return params;
    }
public:
    estimator(const vector<paraminfo>& aparams) :
        params(aparams)
        {}

    ///  \p initialpars.size()==0 means no ninitial parameter
    virtual double estimate(vector<paramresult>& result,
                   bool stds = false) = 0;

};


class nloptuser: public virtual object
{
    nlopt::algorithm alg;
    bool ismax;
    double maxtime;
    double xtolrel;
    double ftolrel;
    double xtolabs;
    double ftolabs;
    bool gradient;

    string what;

    static double eval(const vector<double> &x,
                     vector<double> &g, void* instance)
    {
        nloptuser* self = (nloptuser*) instance;

        double res = na;
        self->what = "";
        try
        {
            res = self->objf(x,g);
        }
        catch(exception& e)
        {
            self->what = e.what();
            throw;
        }
        return res;
    }
protected:
    nloptuser(nlopt::algorithm aalg, bool amax):
       alg(aalg), ismax(amax), maxtime(HUGE_VAL),
       xtolrel(na), ftolrel(na), xtolabs(na), ftolabs(na),
       gradient(true)
       {}
    virtual double objf(const vector<double> &x,
                                  vector<double> &g) = 0;

    double opt(const vector<paraminfo>& params,
                             vector<paramresult>& result  )
    {
        unsigned int q = params.size();

        vector<double> x(q);
        vector<double> lower(q);
        vector<double> upper(q);

        for(unsigned int i=0; i<q; i++)
        {
            double in = params[i].initial;
            x[i] = in==na ? 0 : in;
            lower[i] = params[i].lower;
            upper[i] = params[i].upper;
        }
        nlopt::opt op( // tbd use alg
          gradient ? nlopt::LN_BOBYQA : nlopt::LN_BOBYQA , q); //LD_LBFGS
        if(xtolrel != na)
            op.set_xtol_rel(xtolrel);
        if(ftolrel != na)
            op.set_ftol_rel(ftolrel);
        if(xtolabs != na)
            op.set_xtol_abs(xtolabs);
        if(ftolabs != na)
            op.set_ftol_abs(ftolabs);

        op.set_maxtime(maxtime);
        op.set_lower_bounds(lower);
        op.set_upper_bounds(upper);
        if(ismax)
            op.set_max_objective(eval, this);
        else
            op.set_min_objective(eval, this);

        double r;
        nlopt::result optr;

        try
        {
           optr = op.optimize(x, r);
clog << "nlopt res: " << optr << endl;
           result.clear();
           for(unsigned int i=0; i<q; i++)
           {
               result.push_back(paramresult(params[i]));
               result[i].value = x[i];
           }
        }
        catch(...) // added in order to bridge "swallowing" oo
                    // exceptions by nlopt
        {
            if(what.length() > 0)
            {
                cerr << "Exception occured: " << what << endl;
                throw runtime_error(what);
            }
            else
                throw;
        }
        const char *msg;
        switch(optr)
        {
            case NLOPT_SUCCESS:
                msg = "Generic success return value.";
                break;
            case NLOPT_STOPVAL_REACHED:
                msg = "Optimization stopped because stopval was reached.";
                break;
            case NLOPT_FTOL_REACHED:
                msg = "Optimization stopped because ftol_rel or ftol_abs was reached.";
                break;
            case NLOPT_XTOL_REACHED:
                msg = "Optimization stopped because xtol_rel or xtol_abs was reached.";
                break;
            case NLOPT_MAXEVAL_REACHED:
                msg = "Optimization stopped because maxeval (above) was reached.";
                break;
            case nlopt::MAXTIME_REACHED:
                 throw runtime_error("nlopt: axtime reached");
            case nlopt::FAILURE:
                 throw runtime_error("nlopt: failure");
            case nlopt::INVALID_ARGS:
                 throw runtime_error("nlopt: invalid args");
            case nlopt::OUT_OF_MEMORY:
                 throw runtime_error("nlopt: out of memory");
            case nlopt::ROUNDOFF_LIMITED:
                 throw runtime_error("nlopt: roundoff error");
            case nlopt::FORCED_STOP:
                 throw runtime_error("nlopt: forced stop");
            default:
                msg = "Unknown nlopt return value";
        }
        if(islogging())
                lout << msg << endl;

        return r;
    }
public:
    void setmaxtime(double m)
    {
        maxtime = m;
    }
    void setftolrel(double f)
    {
        ftolrel = f;
    }
    void setxtolrel(double f)
    {
        xtolrel = f;
    }
    void setftolabs(double f)
    {
        ftolabs = f;
    }
    void setxtolabs(double f)
    {
        xtolabs = f;
    }
    void setgradient(bool g)
    {
        gradient = g;
    }
};

//class sample
//{
//protected:
//    virtual unsigned int getn() = 0;
//    virtual unsigned int getdim() = 0;
//
//};


/// Empirical distribution on {0,1,...}

class empdistn
{
    vector<unsigned int> dist;
    unsigned int n;
public:
    empdistn() : n(0) {}
    void add(unsigned int obs)
    {
        if(obs > dist.size())
        {
            unsigned int olds = dist.size();
            dist.resize(olds+10);
            for(unsigned int i=olds; i<dist.size(); i++)
                dist[i] = 0;
            dist[obs]++;
            n++;
        }
    }
    void getdist(vector<double> d)
    {
        if(n==0)
        {
            d.resize(0);
        }
        else
        {
            d.resize(dist.size());
            for(unsigned int i=0; i<dist.size(); i++)
                d[i] = (double) dist[i] / (double) n;
        }
    }
};

/// Class computing MLE estinate


class mle : public estimator, public nloptuser
{
    // state variables
    boost::numeric::ublas::matrix<double> J;
    bool computeJ;

    double objf(const vector<double> &x, vector<double> &g)
    {

        unsigned int R = x.size();
        beforeloglikeval(x);
        for(unsigned int j=0; j<g.size(); j++)
            g[j] = 0;
        boost::numeric::ublas::matrix<double>
          I(boost::numeric::ublas::zero_matrix<double>(R,R));
        double s=0;

        unsigned int M = getN();
        for(unsigned int i=0; i<M; i++)
        {
            vector<double> grad(R);
            s += evallogdensity(i,x,grad);
            for(unsigned int j=0; j<g.size(); j++)
                g[j] += grad[j];
            if(computeJ)
                for(unsigned int j=0; j<R; j++)
                {
                    for(unsigned int k=0; k<R; k++)
                        I(j,k)+=grad[j]*grad[k] / (double) M;
                }
        }
        afterloglikeval(x, g, s);
        if(computeJ)
        {
            try
            {
                tools::invermatrix(I, J);
            }
            catch(...)
            {
                if(islogging())
                {
                    lout << "Cannot invert matrix" << endl;
                    for(unsigned int i=0; i<R; i++)
                    {
                        for(unsigned int j=0; j<R; j++)
                            lout << I(i,j) << ",";
                        lout << endl;
                    }
                }
                for(unsigned int i=0; i<R; i++)
                    for(unsigned int j=0; j<R; j++)
                        J(i,j) = na;
            }
        }
        if(logging >= extendedlogging)
        {
            lout << "#" << " LL=" << s << " args=(";
            for(unsigned int j=0;j<x.size();j++)
                lout << x[j] << " ";
            lout << ")";
            if(g.size())
            {
               lout << " grad=(";
               for(unsigned int j=0;j<g.size();j++)
                  lout << g[j] << " ";
               lout << ")";
            }
            lout << endl;
        }
//clog << "l(" << xx[0]<<","<< xx[1]<<"," << xx[2] << ") "
//        << "-l(" << x[0]<<","<< x[1]<<"," << x[2] << ") "
//        << ss << "-" << s << "..." << (ss-s) / 0.0001 << "=mlr=" << g[2] ;
//clog << " returning "<< s << endl;
        return s;
    }

protected:
    /// called once before repeated evaluation of evallogdensity (with the same params but different i)
    virtual int getN() = 0;
    virtual void beforeloglikeval(const vector<double> aparams) {}
    virtual double evallogdensity(int i, const vector<double>& aArgs,
                                                  vector<double>& g) = 0;
    virtual void afterloglikeval(const vector<double>& aparams,
                         const vector<double>& grad, double loglik) {}
public:

    mle(const vector<paraminfo>& aparams):
        estimator(aparams), nloptuser(nlopt::LD_LBFGS, true),
        J(aparams.size(),aparams.size()),
        computeJ(false)
      {}

    ///  \p initialpars.size()==0 means no ninitial parameter
    double estimate(vector<paramresult>& result,
                    bool stds = false)
    {
        unsigned int R=getparams().size();
        if(islogging())
        {
            lout << "MLE: Estimating " << R << " parameters using "
                << getN() << " observations." << endl;
        }

        computeJ = false;

        double loglik = opt(getparams(),result);

        if(stds)
        {
            computeJ = true;
            vector<double> x(R);
            vector<double> g(R);

            for(unsigned int i=0; i<R; i++)
                x[i] = result[i].value;

            loglik = objf(x,g);

            for(unsigned int i=0; i<R; i++)
            {
                double K = J(i,i);
                if(K!=na)
                    result[i].std=sqrt(K / (double) getN());
                result[i].grad = g[i];
            }
        }
        if(islogging())
        {
            lout << "Optimal params:" << endl;
            for(unsigned int j=0; j<result.size(); j++)
                lout << result[j] << endl;
            lout << endl;
        }

        return loglik;
    }

    vector<vector<double>> getJ()
    {
        unsigned int R=getparams().size();
        vector<vector<double>> res(R,vector<double>(R));
        for(unsigned i=0; i<R; i++)
            for(unsigned j=0; j<R; j++)
                res[i][j]=J(i,j);
        return res;
    }

/*	virtual string description()
    {
        string s = "MLE";
        return s;
    }*/
/*	void result(ostream& ofs, mle* refm = 0, bool tex = false)
    {
        using namespace std;
        if(!estimated)
            throw 1;

        ostream& t = ofs;

        if(!tex)
            ofs << endl << description() << " n=" << getN() << " r=" << getR() << endl << endl;

        const boost::numeric::ublas::matrix<double> Jm = J();

        for(unsigned int i=0; i<getR(); i++)
        {
            string s;
            double p = params[i];
            double sd = sqrt(Jm(i,i) / (double) getN());
            if(!tex)
                ofs << getlabel(i) << " " << p << "(" << sd << ") ";

            double z = fabs(p)/sd;

            if(!tex)
                ofs << z;

            string stars;
            if(z > 1.96)
                stars+= "*";
            if(z > 2.576)
                stars+= "*";
            if(z > 3.29)
                stars+= "*";

            if(!tex)
                ofs << stars << endl;
//			*mlog << m.getlabel(i) << "," << p << "," << sd
//					<< "," << z << "," << stars << endl;
            if(tex)
            {
                t << "$" << getlabel(i) << resetiosflags( ios::floatfield )
                  << "$ = $" << p << "$ ($" << sd << "$)" << stars << "\\\\" << endl;
            }
        }

        if(!tex)
            ofs << endl <<  "loglik=" << loglik() << endl;
    //    *mlog << endl <<  "loglik," << loglik << endl;


    //	*mlog << "rho," << rho << endl;

        double refll;
        int refdf;
        string rd;
        double rho;
        if(refm)
        {
            refll = refm->loglik();
            refdf = getR()-refm->getR();
            rd = refm->description();
            rho = 1.0 - loglik() / refm->loglik();
        }
        else
        {
            refll = loglik();
            refdf = getR();
            rd = "zero parameters";
            rho = 0;
        }
        if(!tex)
            ofs << "rho = " << rho << endl;

        if(!tex)
            ofs << "ratio = " << -2*(refll - loglik()) << " df=" << refdf << endl << endl;
        if(0 && tex)
        {
            t << "\\hline" << endl;
                        t << "\\end{tabular}" << endl << endl;
                        t << "\\vspace{5mm}" << endl << endl;

            t << "\\begin{tabular}{lrlr}" << endl;
            t << "$\\rho$ & " << rho << " & observations &" << getN() << "\\\\" << endl;
            t << "likelihood ratio & " << -2*(refll - loglik()) << " & d.f & "
                    << refdf << "\\\\" << endl;
            t << "\\end{tabular}" << endl;
        //    t << "%(" << m.description() << ") vs (" << rd << ")" << endl;
            t << "\\end{center}" << endl<< endl<< endl<< endl;
        }
    }*/

    double loglik(const vector<double> &x, vector<double> &g)
    {
        return objf(x,g);
    }

};

/// A class encapsulting calls of GRETL

class gretl
{
    string dir;
    string tmpfn()
    {
        return dir+"/_gretl_tmp";
    }
    string scriptresult;
public:
    gretl(const string& adir): dir(adir) {}
    string scdir()
    {
        return dir;
    }
    void runscript(const string& scfn)
    {
        using namespace std;
        string cmd = "gretlcli -b \"" + scfn + "\" > " + tmpfn();
        cout << cmd << endl;
        if(system(cmd.c_str()))
        {
            cerr << "Failed to run " << cmd << endl;
            throw 1;
        }
        ifstream p(tmpfn().c_str());
        if(!p)
        {
            cerr << "Error opening " << tmpfn() << endl;
            throw 1;
        }

        getline(p,scriptresult,'\0');
    }
    string getscriptoutput()
    {
        ifstream s(tmpfn().c_str());
        if(!s)
        {
                        cerr << "Error opening " << tmpfn().c_str() << endl;
            throw 1;
        }
        string str;
        getline(s,str);
        return str;
    }
    double findpar(const string& parname) // assumes that "? parname" appears in the output
    {
        string label = "? " + parname;
        unsigned int pos;
        if((pos = scriptresult.find(label))==string::npos)
        {
            cerr << "Cannot find parameter " << parname << " in gretl script output " << endl;
            throw 1;
        }
        while(scriptresult[pos++] != '\n')
            ;
        return atof(scriptresult.c_str()+pos);
    }
};

/// A class handling CSV input

class csvRow
{
    public:
        string const& operator[](size_t index) const
        {
            return m_data[index];
        }
        size_t size() const
        {
            return m_data.size();
        }
        void readNextRow(istream& str)
        {
            string         line;
            getline(str,line);

            stringstream   lineStream(line);
            string         cell;

            m_data.clear();
            while(getline(lineStream,cell,','))
            {
                m_data.push_back(cell);
            }
        }
    private:
        vector<string> m_data;
};

/// Operator reading a single csv row
/** TODO: handle end of line
*/

inline istream& operator>>(istream& str,csvRow& data)
{
    data.readNextRow(str);
    return str;
}


inline vector<double> operator+ (const vector<double>& x, const vector<double>& y)
{
    unsigned int q=x.size();
    if(y.size() != q)
        throw runtime_error("Different sizes in vector + oparation");
    vector<double> r(q);
    for(unsigned int s=0; s<q; s++)
        r[s] = x[s] + y[s];
    return r;
}



inline vector<double> operator- (const vector<double>& x, const vector<double>& y)
{
    unsigned int q=x.size();
    vector<double> r(q);
    for(unsigned int s=0; s<q; s++)
        r[s] = x[s] - y[s];
    return r;
}


inline vector<double> operator* (double c, const vector<double>& x)
{
    unsigned int q=x.size();
    vector<double> r(q);
    for(unsigned int s=0; s<q; s++)
        r[s] = c *x[s];
    return r;
}

inline vector<double> operator* (const vector<double>& y, double c)
{
    return c*y;
}

inline vector<double> operator/ (const vector<double>& y, double c)
{
    return y*(1.0/c);
}


inline vector<double> const_vector(double c, size_t q)
{
    vector<double> r(q);
    for(unsigned int s=0; s<q; s++)
        r[s] = c;
    return r;
}



inline vector<double> zero_vector(size_t q)
{
    return const_vector(0,q);
}


/*
class ols
{
    vector<double> **X;
    int r;
    const vector<double>& Y;
    bool c;
    bool estimated;
    boost::numeric::ublas::vector<double> bresult;
    vector<double> mhatY;
    vector<double> mhate;
    double mSSE;
public:
    ols(vector<double>** aX, int ar, vector<double>& aY, bool aConst = true):
        X(aX), r(ar), Y(aY), c(aConst), estimated(false),
        bresult(ar + (aConst?1:0)),
        mhatY(aY.size()), mhate(aY.size()) {}
    double b(int i) { if(estimated) return bresult[i]; else throw 1;}
    const vector<double>& hatY() { if(estimated) return mhatY; else throw 1;}
    const vector<double>& hate() { if(estimated) return mhate; else throw 1;}
    double SSE() { if(estimated) return mSSE; else throw 1;};
    void estimate()
    {
        int cshift = c ? 1 : 0;
        using namespace boost::numeric::ublas;
        int k=bresult.size();
        int n=Y.size();
        matrix<double> XtX(k,k);
        boost::numeric::ublas::vector<double> XtY(k);
        for(int i=0; i<k; i++)
        {
            double t=0;
            for(int h=0; h<n; h++)
            {
                int ii = i-cshift;
                double x = ii < 0 ? 1 : (*X[ii])[h];
                t +=x*Y[h];
            }
            XtY[i] = t;

            for(int j=0; j<k; j++)
            {
                double s=0;
                for(int h=0; h<n; h++)
                {
                    int ii = i-cshift;
                    double x1 = ii < 0 ? 1 : (*X[ii])[h];

                    int jj = j-cshift;
                    double x2 = jj < 0 ? 1 : (*X[jj])[h];

                    s += x1*x2;
                }
                XtX(i,j)=s;
            }
        }

        matrix<double> XtXinv(k,k);
        try
        {
            tools::InvertMatrix(XtX,XtXinv);
        }
        catch(...)
        {
            cerr << "unable to invert: " << endl;
            cerr << XtX << endl;
            throw;
        }
        bresult = prod(XtXinv, XtY);

        mSSE = 0;
        for(int i=0; i<n; i++)
        {
            double y = c ? bresult[0] : 0;
            for(int j=0; j<r; j++)
                y += bresult[j+cshift] * (*X[j])[i];
            mhatY[i] = y;
            double e = Y[i] - y;
            mhate[i] = e;
            mSSE += e * e;
        }
        estimated = true;
    }
};
*/

class uncertain
{
public:
    uncertain(const vector<double>& x, const matrix& var) : fx(x), fvar(var)
        {}
    uncertain(const vector<double>& x) : fx(x), fvar(matrix(x.size(),x.size(),0))
        {}
    uncertain() : fx(0), fvar(matrix(0,0,0))
        {}
    vector<double> x() const { return fx;}
    unsigned dim() const { return fx.size(); }
    matrix var() const { return fvar; }
private:
    vector<double> fx;
    matrix fvar;
};

#endif // OLDORPP_HPP_
