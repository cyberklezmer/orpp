#ifndef MATRICES_HPP
#define MATRICES_HPP
#include <Eigen/Eigenvalues>
#include "orpp.hpp"

namespace orpp
{



//inline std::ostream& operator<<(std::ostream& str,const dmatrix& r);

using dvector=Eigen::VectorXd;//<double, Eigen::Dynamic>;
using dmatrix=Eigen::Matrix<double, Eigen::Dynamic,Eigen::Dynamic>;



inline double radius(const dmatrix& m)
{
    assert(m.rows()==m.cols());
    using namespace Eigen;
    EigenSolver<dmatrix> s(m);
    Vector<std::complex<double>,Dynamic> lambdas = s.eigenvalues();
    double lambda = 0;
    for(unsigned i=0; i<m.rows(); i++)
        if(fabs(lambdas(i))>lambda)
            lambda = fabs(lambdas(i));
    return lambda;
}

inline dmatrix pseudoinverse(const dmatrix& m)
{
    //tbd https://stackoverflow.com/questions/44465197/eigen-library-pseudo-inverse-of-dmatrix-matlab-pinv

    using namespace Eigen;
    assert(m.rows()==m.cols());

    dmatrix B(m.rows(),m.cols());
    B = m.completeOrthogonalDecomposition().pseudoInverse();
    return B;
}



inline dvector stackv(const dmatrix& x, const dmatrix& y)
{
    assert(x.cols()==1);
    assert(x.cols()==1);

    dmatrix ret(x.rows()+y.rows(),x.cols());
    unsigned i=0;
    for(; i<x.rows(); i++)
        for(unsigned j=0; j<x.cols(); j++)
            ret(i,j) = x(i,j);
    for(unsigned k=0; k<y.rows(); k++,i++)
        for(unsigned j=0; j<y.cols(); j++)
            ret(i,j) = y(k,j);
    return ret;
}

inline dmatrix stackm(const dmatrix& x, const dmatrix& y)
{
    assert(x.cols()==y.cols());

    dmatrix ret(x.rows()+y.rows(),x.cols());
    unsigned i=0;
    for(; i<x.rows(); i++)
        for(unsigned j=0; j<x.cols(); j++)
            ret(i,j) = x(i,j);
    for(unsigned k=0; k<y.rows(); k++,i++)
        for(unsigned j=0; j<y.cols(); j++)
            ret(i,j) = y(k,j);
    return ret;
}

inline dmatrix block(const dmatrix& x, const dmatrix& y,
                    const dmatrix& z, const dmatrix& zz)
{
    assert(x.rows()==y.rows());
    assert(z.rows()==zz.rows());
    assert(x.cols()==z.cols());
    assert(y.cols()==zz.cols());

    dmatrix ret(x.rows()+z.rows(),x.cols()+y.cols());
    unsigned i=0;
    for(; i<x.rows(); i++)
    {
        unsigned j=0;
        for(j=0; j<x.cols(); j++)
            ret(i,j) = x(i,j);
        for(unsigned k=0; k<y.cols(); k++, j++)
            ret(i,j) = y(i,k);
    }
    for(unsigned k=0; k<z.rows(); k++,i++)
    {
        unsigned j=0;
        for(; j<z.cols(); j++)
            ret(i,j) = z(k,j);
        for(unsigned kk=0; kk<zz.cols(); kk++,j++)
            ret(i,j) = zz(k,kk);
    }
    return ret;
}

inline dmatrix diag(const dvector& v)
{
    dmatrix ret(v.size(), v.size());
    ret.setZero();
    for(unsigned i=0; i<v.size(); i++)
        ret(i,i) = v[i];
    return ret;
}


inline dvector dv(const std::vector<double> &d )
{
    dvector ret(d.size());
    for(unsigned i=0; i<d.size(); i++)
        ret[i] = d[i];
    return ret;
}

inline std::vector<double> vd(const dvector& dv)
{
    std::vector<double> ret(dv.size());
    for(unsigned i=0; i<dv.size(); i++)
        ret[i] = dv[i];
    return ret;
}

inline auto m2csv(const dmatrix& m)
{
    Eigen::IOFormat f(Eigen::StreamPrecision, Eigen::DontAlignCols,",");
    return m.format(f);
}












} // namespace

#endif // MATRICES_HPP

/// Routine computing dmatrix inversion

/**
 *  Taken from http://www.crystalclearsoftware.com/cgi-bin/boost_wiki/wiki.pl?LU_Matrix_Inversion
 *  Uses lu_factorize and lu_substitute in uBLAS to invert a dmatrix

template<class T>
bool invermatrix (const boost::numeric::ublas::dmatrix<T>& input, boost::numeric::ublas::dmatrix<T>& inverse) {
    using namespace boost::numeric::ublas;
    typedef boost::numeric::ublas::permutation_matrix<size_t> pmatrix;
    // create a working copy of the input
    boost::numeric::ublas::dmatrix<T> A(input);
    // create a permutation dmatrix for the LU-factorization
    pmatrix pm(A.size1());

    // perform LU-factorization
    int res = lu_factorize(A,pm);
       if( res != 0 ) return false;

    // create identity dmatrix of "inverse"
    inverse.assign(boost::numeric::ublas::identity_matrix<T>(A.size1()));

    // backsubstitute to get the inverse
    lu_substitute(A, pm, inverse);

    return true;
}

} // namespace

    dmatrix inverse()
    {
        assert(cols()==rows());
        unsigned n = cols();
        dmatrix ret(n,n);
        boost::numeric::ublas::dmatrix<double> A(n,n);
        for(unsigned i=0; i<n; i++)
            for(unsigned j=0; j<n; j++)
                A(i,j)=(*this)(i,j);
        using namespace boost::numeric::ublas;
        typedef permutation_matrix<size_t> pmatrix;
        // create a permutation dmatrix for the LU-factorization
        pmatrix pm(A.size1());

        // perform LU-factorization
        int res = lu_factorize(A,pm);
           if( res != 0 )
           {
               std::cerr << "cannot invert dmatrix" << std::endl;
               std::cerr << *this << std::endl;
               throw "cannot invert dmatrix";
           }

        boost::numeric::ublas::dmatrix<double> inverse(n,n);
        // create identity dmatrix of "inverse"
        inverse.assign(boost::numeric::ublas::identity_matrix<double>(A.size1()));

        // backsubstitute to get the inverse
        lu_substitute(A, pm, inverse);

        for(unsigned i=0; i<n; i++)
            for(unsigned j=0; j<n; j++)
                ret(i,j)=inverse(i,j);
        return ret;
    }


    bool iszero()
    {
        for(unsigned i=0; i<rows(); i++)
            for(unsigned j=0; j<cols(); j++)
                if((*this)(i,j) != 0)
                    return false;
        return true;
    }

inline double cdist(const dmatrix& a, const dmatrix& b)
{
    assert(a.rows()==b.rows());
    assert(a.cols()==b.cols());
    double sum = 0;
    for(unsigned i=0;i<a.rows();i++)
        for(unsigned j=0;j<a.cols();j++)
            sum += fabs(a(i,j)-b(i,j));
    return sum;
}

template <typename D>
inline D sum(const std::vector<std::vector<D>> c)
{
    D sum = 0;
    for(unsigned i=0; i<c.size(); i++)
        for(unsigned j=0; j<c[i].size(); j++)
            sum += c[i][j];
    return sum;
}

template <typename D>
inline D sum(const std::vector<D> c)
{
    D sum = 0;
    for(unsigned i=0; i<c.size(); i++)
            sum += c[i];
    return sum;
}


template<typename T>
inline void mcsvout(std::ostream& os, const std::vector<std::vector<T>>& s)
{
    for(unsigned i=0; i<s.size(); i++)
    {
        for(unsigned j=0; j<s[i].size(); j++)
            os << s[i][j] << ",";
        os << std::endl;
    }
    os << std::endl;
}

*/
