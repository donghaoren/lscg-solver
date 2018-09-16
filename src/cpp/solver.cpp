#include "api.h"
#include "common.h"
#include "impl_solver.h"

EXPORT solver_t solver_create() { return (solver_t) new SolverImpl; }

EXPORT void solver_destroy(solver_t solver) { delete (SolverImpl *)solver; }

EXPORT void solver_add_variable(solver_t solver, int variable_name, number_t value, bool edit) {
    ((SolverImpl *)solver)->add_variable(variable_name, value, edit);
}

EXPORT void solver_make_constant(solver_t solver, int variable_name) {
    ((SolverImpl *)solver)->make_constant(variable_name);
}

EXPORT solver_constraint_t solver_add_constraint(solver_t solver, int strength, number_t bias,
                                                 int count, int *variable_names,
                                                 number_t *weights) {
    return (solver_constraint_t)(
        ((SolverImpl *)solver)->add_constraint(strength, bias, count, variable_names, weights));
}

EXPORT void solver_add_constraint_coefficient(solver_t solver, solver_constraint_t constraint,
                                              int variable_name, number_t weight) {
    Constraint *c = (Constraint *)constraint;
    c->weights.push_back(weight);
    c->variable_names.push_back(variable_name);
}

EXPORT void solver_set_constraint_strength(solver_t solver, solver_constraint_t constraint,
                                           int strength) {
    Constraint *c = (Constraint *)constraint;
    c->strength = strength;
}
EXPORT void solver_set_constraint_bias(solver_t solver, solver_constraint_t constraint,
                                       number_t bias) {
    Constraint *c = (Constraint *)constraint;
    c->bias = bias;
}
EXPORT void solver_clear_constraint_coefficients(solver_t solver, solver_constraint_t constraint) {
    Constraint *c = (Constraint *)constraint;
    c->weights.clear();
    c->variable_names.clear();
}

EXPORT void solver_set_value(solver_t solver, int variable_name, number_t value) {
    ((SolverImpl *)solver)->set_value(variable_name, value);
}

EXPORT number_t solver_get_value(solver_t solver, int variable_name) {
    return ((SolverImpl *)solver)->get_value(variable_name);
}

EXPORT void solver_set_values(solver_t solver, int count, const int *variable_names,
                              const number_t *values) {
    ((SolverImpl *)solver)->set_values(count, variable_names, values);
}

EXPORT void solver_get_values(solver_t solver, int count, const int *variable_names,
                              number_t *output) {
    ((SolverImpl *)solver)->get_values(count, variable_names, output);
}

EXPORT void solver_solve(solver_t solver) { ((SolverImpl *)solver)->solve(); }
EXPORT void solver_compute_loss(solver_t solver) { ((SolverImpl *)solver)->computeLoss(); }

EXPORT void solver_set_attribute_i(solver_t solver, int attribute, int value) {
    ((SolverImpl *)solver)->set_attribute_i(attribute, value);
}
EXPORT void solver_set_attribute_f(solver_t solver, int attribute, number_t value) {
    ((SolverImpl *)solver)->set_attribute_f(attribute, value);
}
EXPORT int solver_get_attribute_i(solver_t solver, int attribute) {
    return ((SolverImpl *)solver)->get_attribute_i(attribute);
}
EXPORT number_t solver_get_attribute_f(solver_t solver, int attribute) {
    return ((SolverImpl *)solver)->get_attribute_f(attribute);
}