#ifndef MATRICES_HPP
#define MATRICES_HPP
#include <Eigen/Eigenvalues>
#include "orpp.hpp"

using namespace orpp;

inline double range(const matrix& mat)
{
    using namespace Eigen;
    assert(mat.r());
    assert(mat.r()==mat.c());

    Eigen::Matrix<double, Dynamic,Dynamic> A(mat.r(),mat.c());
    for(unsigned i=0; i<mat.r(); i++)
        for(unsigned j=0; j<mat.c(); j++)
            A(i,j)=mat(i,j);

   EigenSolver<Eigen::Matrix<double, Dynamic,Dynamic> > s(A);
   Vector<complex<double>,Dynamic> lambdas = s.eigenvalues();
   double lambda = 0;
   for(unsigned i=0; i<mat.r(); i++)
       if(fabs(lambdas(i))>lambda)
           lambda = fabs(lambdas(i));
   return lambda;
}

/// tbd to orpp
inline double det(const matrix& mat)
{
    using namespace Eigen;
    assert(mat.r()==mat.c());

    Eigen::Matrix<double, Dynamic,Dynamic> A(mat.r(),mat.c());
    for(unsigned i=0; i<mat.r(); i++)
        for(unsigned j=0; j<mat.c(); j++)
            A(i,j)=mat(i,j);

    return A.determinant();
}


inline matrix pseudoinverse(const matrix& mat)
{
//tbd https://stackoverflow.com/questions/44465197/eigen-library-pseudo-inverse-of-matrix-matlab-pinv

    using namespace Eigen;
    assert(mat.r()==mat.c());

    Eigen::Matrix<double, Dynamic,Dynamic> A(mat.r(),mat.c());
    for(unsigned i=0; i<mat.r(); i++)
        for(unsigned j=0; j<mat.c(); j++)
            A(i,j)=mat(i,j);

    Eigen::Matrix<double, Dynamic,Dynamic> B(mat.r(),mat.c());
    B = A.completeOrthogonalDecomposition().pseudoInverse();
    matrix ret(mat.r(),mat.c());
    for(unsigned i=0; i<mat.r(); i++)
        for(unsigned j=0; j<mat.c(); j++)
            ret(i,j)=B(i,j);
    return ret;
}


#endif // MATRICES_HPP
