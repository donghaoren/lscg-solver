#ifndef WASMSOLVER_API_H
#define WASMSOLVER_API_H

typedef double number_t;
typedef int length_t;

extern "C" {

// Memory allocator

void *memory_alloc(length_t length);
void memory_free(float *pointer);

// Linear algebra

typedef struct matrix__ *matrix_t;

matrix_t linalg_matrix_create();
void linalg_matrix_init(matrix_t matrix, length_t rows, length_t cols, number_t *values);
void linalg_matrix_fill(matrix_t matrix, number_t value);
void linalg_matrix_destroy(matrix_t matrix);

number_t *linalg_matrix_data(matrix_t matrix);
length_t linalg_matrix_size(matrix_t matrix);
length_t linalg_matrix_rows(matrix_t matrix);
length_t linalg_matrix_cols(matrix_t matrix);
length_t linalg_matrix_col_stride(matrix_t matrix);
length_t linalg_matrix_row_stride(matrix_t matrix);

number_t linalg_matrix_norm(matrix_t matrix);
number_t linalg_matrix_l1_norm(matrix_t matrix, int p);

void linalg_matrix_add(matrix_t dest, matrix_t a, matrix_t b);
void linalg_matrix_sub(matrix_t dest, matrix_t a, matrix_t b);
void linalg_matrix_scale(matrix_t dest, matrix_t a, number_t scaler);
void linalg_matrix_add_scale(matrix_t dest, matrix_t a, number_t a_scale, matrix_t b,
                             number_t b_scale);
void linalg_matrix_emul(matrix_t dest, matrix_t a, matrix_t b);
void linalg_matrix_ediv(matrix_t dest, matrix_t a, matrix_t b);
void linalg_matrix_mmul(matrix_t dest, matrix_t a, matrix_t b);

void linalg_solve_linear_system(matrix_t X, matrix_t ker, matrix_t A, matrix_t B);

// Constraint solver

typedef struct solver__ *solver_t;

solver_t solver_create();
void solver_destroy(solver_t solver);

const int kSolverStrength_HARD = 0;
const int kSolverStrength_STRONG = 1;
const int kSolverStrength_MEDIUM = 2;
const int kSolverStrength_WEAK = 3;
const int kSolverStrength_WEAKER = 4;
const int kSolverStrength_DISABLED = 10;

const int kSolverAttribute_NUM_ITERATIONS = 1;
const int kSolverAttribute_TOLERANCE = 2;
const int kSolverAttribute_FLAGS = 3;

const int kSolverAttribute_NUM_VARIABLES = 10;
const int kSolverAttribute_NUM_CONSTRAINTS = 11;
const int kSolverAttribute_MAX_ITERATIONS = 12;
const int kSolverAttribute_ERROR = 13;

const int kSolverAttribute_HARD_LOSS = 20;
const int kSolverAttribute_SOFT_LOSS = 21;

const int kSolverAttribute_REGULARIZER_WEIGHT = 30;

const int kSolverFlag_DEFAULT = 0;
const int kSolverFlag_REDUCE = 1 << 1;
const int kSolverFlag_LAGRANGE = 1 << 2;

void solver_add_variable(solver_t solver, int variable_name, number_t value, bool edit);
void solver_make_constant(solver_t solver, int variable_name);
int solver_add_constraint(solver_t solver, int strength, number_t bias, int count,
                          int *variable_names, number_t *weights);
void solver_add_constraint_coefficient(solver_t solver, int constraint, int variable_name,
                                       number_t weight);
void solver_set_constraint_strength(solver_t solver, int constraint, int strength);
void solver_set_constraint_bias(solver_t solver, int constraint, number_t bias);
void solver_clear_constraint_coefficients(solver_t solver, int constraint);
void solver_set_value(solver_t solver, int variable_name, number_t value);
number_t solver_get_value(solver_t solver, int variable_name);
void solver_set_values(solver_t solver, int count, const int *variable_names,
                       const number_t *values);
void solver_get_values(solver_t solver, int count, const int *variable_names, number_t *output);

void solver_set_attribute_i(solver_t solver, int attribute, int value);
void solver_set_attribute_f(solver_t solver, int attribute, number_t value);
int solver_get_attribute_i(solver_t solver, int attribute);
number_t solver_get_attribute_f(solver_t solver, int attribute);

void solver_solve(solver_t solver);
void solver_compute_loss(solver_t solver);

} // extern "C"

#endif