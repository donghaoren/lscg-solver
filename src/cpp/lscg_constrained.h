#ifndef LSCG_SOLVER_LSCG_CONSTRAINED_H
#define LSCG_SOLVER_LSCG_CONSTRAINED_H

#include "third_party.h"

namespace Eigen {

class CLSCGIdentityPreconditioner {
  public:
    CLSCGIdentityPreconditioner() {}

    template <typename MatrixType> explicit CLSCGIdentityPreconditioner(const MatrixType &) {}

    template <typename MatrixType> CLSCGIdentityPreconditioner &analyzePattern(const MatrixType &) {
        return *this;
    }

    template <typename MatrixType> CLSCGIdentityPreconditioner &factorize(const MatrixType &) {
        return *this;
    }

    template <typename MatrixType> CLSCGIdentityPreconditioner &compute(const MatrixType &) {
        return *this;
    }

    template <typename Rhs> inline const Rhs &solve(const Rhs &b) const { return b; }

    ComputationInfo info() { return Success; }
};

template <typename MatrixType, typename Rhs, typename Dest, typename Preconditioner>
void constrained_least_square_conjugate_gradient(const MatrixType &matA, const MatrixType &matC,
                                                 const Rhs &rhsB, const Rhs &rhsD, Dest &x,
                                                 const Preconditioner &precond, Index &iters,
                                                 typename Dest::RealScalar &tol_error) {
    using std::abs;
    using std::sqrt;
    typedef typename Dest::RealScalar RealScalar;
    typedef typename Dest::Scalar Scalar;
    typedef Matrix<Scalar, Dynamic, 1> VectorType;

    RealScalar tol = tol_error;
    Index maxIters = iters;

    // Index m = mat.rows(), n = mat.cols();
    Index nX = matA.cols();
    Index nLambda = matC.rows();

    VectorType lambda(nLambda);
    lambda.setZero();

    // A'b - A'Ax - C'lambda
    VectorType residualX = matA.adjoint() * (rhsB - matA * x) - matC.adjoint() * lambda;
    // D - C'x
    VectorType residualLambda = rhsD - matC * x;

    RealScalar rhsNorm2 = (matA.adjoint() * rhsB).squaredNorm() + rhsD.squaredNorm();
    ;
    if (rhsNorm2 == 0) {
        x.setZero();
        iters = 0;
        tol_error = 0;
        return;
    }
    RealScalar threshold = tol * tol * rhsNorm2;
    RealScalar residualNorm2 = residualX.squaredNorm() + residualLambda.squaredNorm();
    if (residualNorm2 < threshold) {
        iters = 0;
        tol_error = sqrt(residualNorm2 / rhsNorm2);
        return;
    }

    VectorType pX(nX);
    VectorType pLambda(nLambda);

    // p = precond.solve(normal_residual); // initial search direction
    pX = residualX;
    pLambda = residualLambda;

    VectorType zX(nX), zLambda(nLambda), tmpA(matA.rows()), tmpC(matC.cols());
    // the square of the absolute value of r scaled by invM
    RealScalar absNew = numext::real(residualX.dot(pX) + residualLambda.dot(pLambda));
    Index i = 0;
    MatrixXd denseA = matA;
    MatrixXd denseC = matC;
    while (i < maxIters) {
        tmpA.noalias() = matA * pX;
        tmpC.noalias() = matC.adjoint() * pLambda;

        // the amount we travel on dir
        Scalar alpha = absNew / (tmpA.squaredNorm() + pX.dot(tmpC) * 2);

        // update solution
        x += alpha * pX;
        lambda += alpha * pLambda;

        // update residual
        residualX -= alpha * (matA.adjoint() * tmpA + tmpC);
        residualLambda -= alpha * (matC * pX);

        // update residualNorm2
        residualNorm2 = residualX.squaredNorm() + residualLambda.squaredNorm();
        if (residualNorm2 < threshold)
            break;

        // z = precond.solve(normal_residual); // approximately solve for "A'A z = normal_residual"
        zX = residualX;
        zLambda = residualLambda;

        RealScalar absOld = absNew;

        // update the absolute value of r
        absNew = numext::real(residualX.dot(zX) + residualLambda.dot(zLambda));

        // calculate the Gram-Schmidt value used to create the new search direction
        RealScalar beta = absNew / absOld;

        // update search direction
        pX = zX + beta * pX;
        pLambda = zLambda + beta * pLambda;
        i++;
    }
    tol_error = sqrt(residualNorm2 / rhsNorm2);
    iters = i;
}

} // namespace Eigen

#endif