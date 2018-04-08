#include "api.h"
#include "common.h"
#include "third_party.h"

using namespace Eigen;

struct MatrixContainer {
    MatrixXd value;
};

EXPORT matrix_t linalg_matrix_create() { return (matrix_t) new MatrixContainer(); }

EXPORT void linalg_matrix_init(matrix_t matrix, length_t rows, length_t cols, number_t *values) {
    auto *m = (MatrixContainer *)matrix;
    m->value.setZero(rows, cols);
    if (values != nullptr) {
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                m->value(i, j) = *values++;
            }
        }
    }
}

EXPORT void linalg_matrix_destroy(matrix_t matrix) { delete (MatrixContainer *)matrix; }

EXPORT number_t *linalg_matrix_data(matrix_t matrix) {
    auto *m = (MatrixContainer *)matrix;
    return m->value.data();
}

EXPORT void linalg_matrix_fill(matrix_t matrix, number_t value) {
    auto *v = (MatrixContainer *)matrix;
    v->value.fill(value);
}

EXPORT length_t linalg_matrix_rows(matrix_t matrix) {
    auto *m = (MatrixContainer *)matrix;
    return m->value.rows();
}

EXPORT length_t linalg_matrix_cols(matrix_t matrix) {
    auto *m = (MatrixContainer *)matrix;
    return m->value.cols();
}

EXPORT length_t linalg_matrix_size(matrix_t matrix) {
    auto *m = (MatrixContainer *)matrix;
    return m->value.size();
}

EXPORT length_t linalg_matrix_row_stride(matrix_t matrix) {
    auto *m = (MatrixContainer *)matrix;
    return m->value.rowStride();
}

EXPORT length_t linalg_matrix_col_stride(matrix_t matrix) {
    auto *m = (MatrixContainer *)matrix;
    return m->value.colStride();
}

EXPORT number_t linalg_matrix_norm(matrix_t matrix) {
    auto *v = (MatrixContainer *)matrix;
    return v->value.norm();
}

EXPORT number_t linalg_matrix_l1_norm(matrix_t matrix, int p) {
    auto *v = (MatrixContainer *)matrix;
    return v->value.lpNorm<1>();
}

EXPORT void linalg_matrix_add(matrix_t dest, matrix_t a, matrix_t b) {
    auto *rd = (MatrixContainer *)dest;
    auto *ra = (MatrixContainer *)a;
    auto *rb = (MatrixContainer *)b;
    rd->value = ra->value + rb->value;
}

EXPORT void linalg_matrix_sub(matrix_t dest, matrix_t a, matrix_t b) {
    auto *rd = (MatrixContainer *)dest;
    auto *ra = (MatrixContainer *)a;
    auto *rb = (MatrixContainer *)b;
    rd->value = ra->value - rb->value;
}

EXPORT void linalg_matrix_emul(matrix_t dest, matrix_t a, matrix_t b) {
    auto *rd = (MatrixContainer *)dest;
    auto *ra = (MatrixContainer *)a;
    auto *rb = (MatrixContainer *)b;
    rd->value.array() = ra->value.array() * rb->value.array();
}

EXPORT void linalg_matrix_ediv(matrix_t dest, matrix_t a, matrix_t b) {
    auto *rd = (MatrixContainer *)dest;
    auto *ra = (MatrixContainer *)a;
    auto *rb = (MatrixContainer *)b;
    rd->value.array() = ra->value.array() / rb->value.array();
}

EXPORT void linalg_matrix_scale(matrix_t dest, matrix_t a, number_t scaler) {
    auto *rd = (MatrixContainer *)dest;
    auto *ra = (MatrixContainer *)a;
    rd->value = ra->value * scaler;
}

EXPORT void linalg_matrix_add_scale(matrix_t dest, matrix_t a, number_t a_scale, matrix_t b,
                                    number_t b_scale) {
    auto *rd = (MatrixContainer *)dest;
    auto *ra = (MatrixContainer *)a;
    auto *rb = (MatrixContainer *)b;
    rd->value = ra->value * a_scale + rb->value * b_scale;
}

EXPORT void linalg_matrix_mmul(matrix_t dest, matrix_t a, matrix_t b) {
    auto *rd = (MatrixContainer *)dest;
    auto *ra = (MatrixContainer *)a;
    auto *rb = (MatrixContainer *)b;
    rd->value = ra->value * rb->value;
}

EXPORT void linalg_solve_linear_system(matrix_t X, matrix_t ker, matrix_t A, matrix_t B) {
    auto *rX = (MatrixContainer *)X;
    auto *rker = (MatrixContainer *)ker;
    auto *rA = (MatrixContainer *)A;
    auto *rB = (MatrixContainer *)B;

    FullPivLU<MatrixXd> lu(rA->value);
    rker->value = lu.kernel();
    rX->value = lu.solve(rB->value);
}
