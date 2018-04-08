#ifndef IMPL_SOLVER_H
#define IMPL_SOLVER_H

#include "api.h"
#include "common.h"

#include <unordered_map>
#include <unordered_set>
#include <vector>

struct Variable {
    int location;
    bool edit;
    number_t value;
};

struct Constraint {
    int strength;
    std::vector<int> variable_names;
    std::vector<number_t> weights;
    number_t bias;
};

struct AttributeValue {
    union {
        int int_value;
        number_t float_value;
    };
};

class SolverImpl {
  public:
    SolverImpl();

    void add_variable(int variable_name, number_t value, bool edit);
    void make_constant(int variable_name);
    Constraint *add_constraint(int strength, number_t bias, int count, int *variable_names,
                               number_t *weights);
    void set_values(int count, const int *variable_names, const number_t *values);
    void get_values(int count, const int *variable_names, number_t *output);
    void set_value(int variable_name, number_t value);
    number_t get_value(int variable_name);

    void solve();

    void solveNormal();
    void solveNormalLagrange();
    void solveWithReduce();

    void computeLoss();

    void set_attribute_i(int name, int value) { _attributes[name].int_value = value; }
    void set_attribute_f(int name, number_t value) { _attributes[name].float_value = value; }
    int get_attribute_i(int name) { return _attributes[name].int_value; }
    number_t get_attribute_f(int name) { return _attributes[name].float_value; }

  private:
    std::unordered_map<int, Variable> _variables;
    std::vector<Constraint> _constraints;
    std::unordered_map<int, AttributeValue> _attributes;
};

#endif