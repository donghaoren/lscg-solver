#ifndef LSCG_SOLVER_LSCG_H
#define LSCG_SOLVER_LSCG_H

#include "third_party.h"

namespace Eigen {

template <typename _Scalar>
class RegularizedLeastSquareDiagonalPreconditioner : public DiagonalPreconditioner<_Scalar> {
    typedef _Scalar Scalar;
    typedef typename NumTraits<Scalar>::Real RealScalar;
    typedef DiagonalPreconditioner<_Scalar> Base;
    using Base::m_invdiag;

    Scalar m_k2;

  public:
    RegularizedLeastSquareDiagonalPreconditioner(Scalar k) : Base() { m_k2 = k * k; }

    template <typename MatType>
    explicit RegularizedLeastSquareDiagonalPreconditioner(const MatType &mat) : Base() {
        compute(mat);
    }

    template <typename MatType>
    RegularizedLeastSquareDiagonalPreconditioner &analyzePattern(const MatType &) {
        return *this;
    }

    template <typename MatType>
    RegularizedLeastSquareDiagonalPreconditioner &factorize(const MatType &mat) {
        // Compute the inverse squared-norm of each column of mat
        m_invdiag.resize(mat.cols());
        if (MatType::IsRowMajor) {
            m_invdiag.setZero();
            for (Index j = 0; j < mat.outerSize(); ++j) {
                for (typename MatType::InnerIterator it(mat, j); it; ++it)
                    m_invdiag(it.index()) += numext::abs2(it.value());
            }
            for (Index j = 0; j < mat.cols(); ++j) {
                Scalar sum = numext::real(m_invdiag(j));
                sum += m_k2;
                if (sum > RealScalar(0))
                    m_invdiag(j) = RealScalar(1) / sum;
            }
        } else {
            for (Index j = 0; j < mat.outerSize(); ++j) {
                RealScalar sum = mat.innerVector(j).squaredNorm();
                sum += m_k2;
                if (sum > RealScalar(0))
                    m_invdiag(j) = RealScalar(1) / sum;
                else
                    m_invdiag(j) = RealScalar(1);
            }
        }
        Base::m_isInitialized = true;
        return *this;
    }

    template <typename MatType>
    RegularizedLeastSquareDiagonalPreconditioner &compute(const MatType &mat) {
        return factorize(mat);
    }

    ComputationInfo info() { return Success; }

  protected:
};

// Minimize ||Ax - b||^2 + k ||x - x0||^2
template <typename MatrixType, typename Rhs, typename Dest, typename Preconditioner>
void regularized_least_square_conjugate_gradient(const MatrixType &mat, const Rhs &rhs, Dest &x,
                                                 Dest &x0, typename MatrixType::Scalar k,
                                                 const Preconditioner &precond, Index &iters,
                                                 typename Dest::RealScalar &tol_error) {
    using std::abs;
    using std::sqrt;
    typedef typename Dest::RealScalar RealScalar;
    typedef typename Dest::Scalar Scalar;
    typedef Matrix<Scalar, Dynamic, 1> VectorType;

    RealScalar tol = tol_error;
    Index maxIters = iters;

    Index m = mat.rows(), n = mat.cols();

    VectorType residual = rhs - mat * x;
    VectorType residual_r = (x0 - x) * k;
    VectorType normal_residual = mat.adjoint() * residual + residual_r * k;

    RealScalar rhsNorm2 = (mat.adjoint() * rhs + x0 * (k * k)).squaredNorm();
    if (rhsNorm2 == 0) {
        x.setZero();
        iters = 0;
        tol_error = 0;
        return;
    }
    RealScalar threshold = tol * tol * rhsNorm2;
    RealScalar residualNorm2 = normal_residual.squaredNorm();
    if (residualNorm2 < threshold) {
        iters = 0;
        tol_error = sqrt(residualNorm2 / rhsNorm2);
        return;
    }

    VectorType p(n);
    p = precond.solve(normal_residual); // initial search direction

    VectorType z(n), tmp(m), tmp_r(n);
    RealScalar absNew = numext::real(
        normal_residual.dot(p)); // the square of the absolute value of r scaled by invM
    Index i = 0;
    while (i < maxIters) {
        tmp.noalias() = mat * p;
        tmp_r.noalias() = k * p;

        Scalar alpha =
            absNew / (tmp.squaredNorm() + tmp_r.squaredNorm()); // the amount we travel on dir
        x += alpha * p;                                         // update solution
        residual -= alpha * tmp;                                // update residual
        residual_r -= alpha * tmp_r;
        normal_residual =
            mat.adjoint() * residual + residual_r * k; // update residual of the normal equation

        residualNorm2 = normal_residual.squaredNorm();
        if (residualNorm2 < threshold)
            break;

        z = precond.solve(normal_residual); // approximately solve for "A'A z = normal_residual"

        RealScalar absOld = absNew;
        absNew = numext::real(normal_residual.dot(z)); // update the absolute value of r
        RealScalar beta =
            absNew /
            absOld; // calculate the Gram-Schmidt value used to create the new search direction
        p = z + beta * p; // update search direction
        i++;
    }
    tol_error = sqrt(residualNorm2 / rhsNorm2);
    iters = i;
}

} // namespace Eigen

#endif