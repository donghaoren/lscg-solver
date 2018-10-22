#include "impl_solver.h"
#include "lscg.h"
#include "lscg_constrained.h"
#include "third_party.h"
#include <iostream>
#include <list>
#include <queue>
#include <set>

using namespace Eigen;

#define EPSILON 1e-8

template <typename MatrixType, typename Rhs, typename Dest>
void lscg_solve_main(const MatrixType &mat, const Rhs &rhs, Dest &x, Dest &x0,
                     typename MatrixType::Scalar k, Index &iters,
                     typename Dest::RealScalar &tol_error) {
    if (k == 0) {
        LeastSquaresConjugateGradient<MatrixType> lscg;
        lscg.setMaxIterations(iters);
        lscg.setTolerance(tol_error);
        lscg.compute(mat);
        x = lscg.solveWithGuess(rhs, x);
        tol_error = lscg.tolerance();
        iters = lscg.iterations();
    } else {
        RegularizedLeastSquareDiagonalPreconditioner<number_t> preconditioner(k);
        preconditioner.compute(mat);
        regularized_least_square_conjugate_gradient(mat, rhs, x, x0, k, preconditioner, iters,
                                                    tol_error);
    }
}

SolverImpl::SolverImpl() {
    _attributes[kSolverAttribute_REGULARIZER_WEIGHT].float_value = 0;
    _attributes[kSolverAttribute_FLAGS].int_value = 0;
}

void SolverImpl::add_variable(int name, number_t value, bool edit) {
    Variable var;
    var.edit = edit;
    var.value = value;
    _variables[name] = var;
}

void SolverImpl::make_constant(int variable_name) { _variables[variable_name].edit = false; }

int SolverImpl::add_constraint(int strength, number_t bias, int count, int *variable_names,
                               number_t *weights) {
    Constraint c;
    c.strength = strength;
    c.bias = bias;
    c.weights = std::vector<number_t>(weights, weights + count);
    c.variable_names = std::vector<int>(variable_names, variable_names + count);
    _constraints.push_back(c);
    return _constraints.size() - 1;
}

void SolverImpl::set_values(int count, const int *variable_names, const number_t *values) {
    for (int i = 0; i < count; i++) {
        _variables[variable_names[i]].value = values[i];
    }
}

void SolverImpl::get_values(int count, const int *variable_names, number_t *output) {
    for (int i = 0; i < count; i++) {
        output[i] = _variables[variable_names[i]].value;
    }
}

void SolverImpl::set_value(int variable_name, number_t value) {
    _variables[variable_name].value = value;
}

number_t SolverImpl::get_value(int variable_name) { return _variables[variable_name].value; }

void SolverImpl::solve() {
    if (_attributes[kSolverAttribute_FLAGS].int_value & kSolverFlag_REDUCE) {
        solveWithReduce();
    } else {
        if (_attributes[kSolverAttribute_FLAGS].int_value & kSolverFlag_LAGRANGE) {
            solveNormalLagrange();
        } else {
            solveNormal();
        }
    }
}

number_t strength_to_scale(int strength) {
    switch (strength) {
    case kSolverStrength_HARD:
        return 10;
    case kSolverStrength_STRONG:
        return 1;
    case kSolverStrength_MEDIUM:
        return 0.1;
    case kSolverStrength_WEAK:
        return 0.01;
    case kSolverStrength_WEAKER:
        return 0.001;
    case kSolverStrength_DISABLED:
        return 0;
    }
    return 0;
}

void SolverImpl::computeLoss() {
    int constraint_count = _constraints.size();
    number_t hard_loss = 0;
    number_t soft_loss = 0;
    for (int i = 0; i < constraint_count; ++i) {
        const Constraint &c = _constraints[i];
        number_t value = c.bias;
        for (int j = 0; j < c.variable_names.size(); j++) {
            const Variable &var = _variables[c.variable_names[j]];
            value += var.value * c.weights[j];
        }
        if (c.strength == kSolverStrength_HARD) {
            hard_loss += value * value;
        } else {
            soft_loss += value * value;
        }
    }
    _attributes[kSolverAttribute_HARD_LOSS].float_value = hard_loss;
    _attributes[kSolverAttribute_SOFT_LOSS].float_value = soft_loss;
}

void SolverImpl::solveNormal() {
    int variable_count = 0;
    for (auto &it : _variables) {
        if (it.second.edit) {
            it.second.location = variable_count++;
        }
    }
    int constraint_count = _constraints.size();

    SparseMatrix<double> A(constraint_count, variable_count);
    VectorXd x(variable_count);
    VectorXd x0(variable_count);
    VectorXd b(constraint_count);

    // Construct A
    std::vector<Triplet<double>> tripletList;
    tripletList.reserve(constraint_count * 5);
    for (int i = 0; i < constraint_count; ++i) {
        const Constraint &c = _constraints[i];
        number_t scaler = strength_to_scale(c.strength);
        number_t bias = c.bias * scaler;
        for (int j = 0; j < c.variable_names.size(); j++) {
            const Variable &var = _variables[c.variable_names[j]];
            if (var.edit) {
                tripletList.push_back(Triplet<double>(i, var.location, c.weights[j] * scaler));
            } else {
                bias += var.value * c.weights[j] * scaler;
            }
        }
        b[i] = -bias;
    }

    A.setFromTriplets(tripletList.begin(), tripletList.end());

    for (auto &it : _variables) {
        if (it.second.edit) {
            x0[it.second.location] = it.second.value;
        }
    }

    Eigen::Index iters = 2 * A.cols();
    double tol_error = NumTraits<number_t>::epsilon();

    if (_attributes.find(kSolverAttribute_TOLERANCE) != _attributes.end()) {
        tol_error = _attributes[kSolverAttribute_TOLERANCE].float_value;
    }
    if (_attributes.find(kSolverAttribute_MAX_ITERATIONS) != _attributes.end()) {
        iters = _attributes[kSolverAttribute_MAX_ITERATIONS].int_value;
    }

    x = x0;
    lscg_solve_main(A, b, x, x0, _attributes[kSolverAttribute_REGULARIZER_WEIGHT].float_value,
                    iters, tol_error);

    _attributes[kSolverAttribute_NUM_VARIABLES].int_value = variable_count;
    _attributes[kSolverAttribute_NUM_CONSTRAINTS].int_value = constraint_count;
    _attributes[kSolverAttribute_NUM_ITERATIONS].int_value = iters;
    _attributes[kSolverAttribute_ERROR].float_value = tol_error;

    for (auto &it : _variables) {
        if (it.second.edit) {
            it.second.value = x[it.second.location];
        }
    }
}

void SolverImpl::solveNormalLagrange() {
    int constraint_count = _constraints.size();
    if (constraint_count == 0) {
        return;
    }

    int soft_constraint_count = 0;
    int soft_triplets = 0;
    int hard_constraint_count = 0;
    int hard_triplets = 0;

    for (int i = 0; i < constraint_count; ++i) {
        const Constraint &c = _constraints[i];
        if (c.strength != kSolverStrength_DISABLED) {
            if (c.strength == kSolverStrength_HARD) {
                hard_constraint_count++;
                hard_triplets += c.variable_names.size();
            } else {
                soft_constraint_count++;
                soft_triplets += c.variable_names.size();
            }
        }
    }

    if (soft_constraint_count == 0 || hard_constraint_count) {
        solveNormal();
        return;
    }

    // Assign variable locations
    int variable_count = 0;
    for (auto &it : _variables) {
        if (it.second.edit) {
            it.second.location = variable_count++;
        }
    }

    // Problem: minimize || b - Ax ||^2, s.t., C x = d
    SparseMatrix<double> matA(soft_constraint_count, variable_count);
    VectorXd rhsB(soft_constraint_count);
    SparseMatrix<double> matC(hard_constraint_count, variable_count);
    VectorXd rhsD(hard_constraint_count);
    VectorXd x(variable_count);
    VectorXd x0(variable_count);

    // Construct A
    std::vector<Triplet<double>> tripletList;
    {
        tripletList.reserve(soft_triplets);
        int index = 0;
        for (const Constraint &c : _constraints) {
            if (c.strength == kSolverStrength_HARD || c.strength == kSolverStrength_DISABLED) {
                continue;
            }
            number_t scaler = strength_to_scale(c.strength);
            number_t bias = c.bias * scaler;
            for (int j = 0; j < c.variable_names.size(); j++) {
                const Variable &var = _variables[c.variable_names[j]];
                if (var.edit) {
                    tripletList.push_back(
                        Triplet<double>(index, var.location, c.weights[j] * scaler));
                } else {
                    bias += var.value * c.weights[j] * scaler;
                }
            }
            rhsB[index] = -bias;
            index += 1;
        }
        matA.setFromTriplets(tripletList.begin(), tripletList.end());
    }
    {
        tripletList.clear();
        tripletList.reserve(hard_triplets);
        int index = 0;
        for (const Constraint &c : _constraints) {
            if (c.strength != kSolverStrength_HARD || c.strength == kSolverStrength_DISABLED) {
                continue;
            }
            number_t scaler = strength_to_scale(c.strength);
            number_t bias = c.bias * scaler;
            for (int j = 0; j < c.variable_names.size(); j++) {
                const Variable &var = _variables[c.variable_names[j]];
                if (var.edit) {
                    tripletList.push_back(
                        Triplet<double>(index, var.location, c.weights[j] * scaler));
                } else {
                    bias += var.value * c.weights[j] * scaler;
                }
            }
            rhsD[index] = -bias;
            index += 1;
        }
        matC.setFromTriplets(tripletList.begin(), tripletList.end());
    }

    for (auto &it : _variables) {
        if (it.second.edit) {
            x0[it.second.location] = it.second.value;
        }
    }

    Index num_iterations = (variable_count + hard_constraint_count) * 2;
    double tol_error = NumTraits<number_t>::epsilon();

    if (_attributes.find(kSolverAttribute_TOLERANCE) != _attributes.end()) {
        tol_error = _attributes[kSolverAttribute_TOLERANCE].float_value;
    }
    if (_attributes.find(kSolverAttribute_MAX_ITERATIONS) != _attributes.end()) {
        num_iterations = _attributes[kSolverAttribute_MAX_ITERATIONS].int_value;
    }

    constrained_least_square_conjugate_gradient(
        matA, matC, rhsB, rhsD, x, CLSCGIdentityPreconditioner(), num_iterations, tol_error);

    _attributes[kSolverAttribute_NUM_VARIABLES].int_value = variable_count;
    _attributes[kSolverAttribute_NUM_CONSTRAINTS].int_value = constraint_count;
    _attributes[kSolverAttribute_NUM_ITERATIONS].int_value = num_iterations;
    _attributes[kSolverAttribute_ERROR].float_value = tol_error;

    for (auto &it : _variables) {
        if (it.second.edit) {
            it.second.value = x[it.second.location];
        }
    }
}

#define MAX_ELIMINATE_SIZE 5

class Snapshot {
    struct ConstraintItem {
        struct ConstraintItemRef {
            ConstraintItem *item;
            int index;

            int value() { return item->items.size(); }

            bool operator<(const ConstraintItemRef &rhs) const {
                return item->items.size() > rhs.item->items.size();
            }
        };

        number_t bias;
        int strength;
        std::map<int, number_t> items;
        int eliminate_variable;
        std::list<ConstraintItem::ConstraintItemRef>::iterator heap_iterator;
        ConstraintItemRef ref;
    };

    struct VariableItem {
        int variable_name;
        int index;
        int location;
        Variable *variable_ptr;
        std::set<int> connected_constraints;
        bool eliminated;
    };

    std::vector<VariableItem> _variables;
    std::unordered_map<int, int> _variable_name_to_index;
    std::vector<ConstraintItem> _constraints;
    // boost::heap::fibonacci_heap<ConstraintItem::ConstraintItemRef> _constraint_queue;
    std::list<ConstraintItem::ConstraintItemRef> _constraint_queue[MAX_ELIMINATE_SIZE + 1];
    std::vector<int> _eliminate_order;
    std::unordered_map<int, AttributeValue> &_attributes;

  public:
    Snapshot(std::unordered_map<int, AttributeValue> &attributes) : _attributes(attributes) {}
    void initialize(std::unordered_map<int, Variable> &variables,
                    const std::vector<Constraint> &constraints) {
        for (auto &variable : variables) {
            if (variable.second.edit) {
                VariableItem item;
                item.index = _variables.size();
                item.variable_name = variable.first;
                item.variable_ptr = &variable.second;
                item.eliminated = false;
                _variables.push_back(item);
                _variable_name_to_index[item.variable_name] = item.index;
            }
        }
        _constraints.reserve(constraints.size());
        int index = 0;
        for (auto &constraint : constraints) {
            if (constraint.strength == kSolverStrength_DISABLED)
                continue;
            _constraints.push_back(ConstraintItem());
            ConstraintItem &item = _constraints[index];
            item.bias = constraint.bias;
            item.strength = constraint.strength;
            item.eliminate_variable = -1;
            for (int j = 0; j < constraint.variable_names.size(); ++j) {
                int variable_name = constraint.variable_names[j];
                if (variables[variable_name].edit) {
                    item.items[_variable_name_to_index[variable_name]] += constraint.weights[j];
                } else {
                    item.bias += constraint.weights[j] * variables[variable_name].value;
                }
            }
            purne_items(item, index);
            for (auto &it : item.items) {
                _variables[it.first].connected_constraints.insert(index);
            }
            item.ref.item = &item;
            item.ref.index = index;
            if (item.strength == kSolverStrength_HARD && item.ref.value() <= MAX_ELIMINATE_SIZE) {
                _constraint_queue[item.ref.value()].push_front(item.ref);
                item.heap_iterator = _constraint_queue[item.ref.value()].begin();
            }
            index += 1;
        }
    }

    void purne_items(ConstraintItem &c, int index) {
        static std::vector<int> purne_list;
        purne_list.clear();
        for (auto iter = c.items.begin(); iter != c.items.end(); ++iter) {
            if (abs(iter->second) < EPSILON) {
                _variables[iter->first].connected_constraints.erase(index);
                purne_list.push_back(iter->first);
            }
        }
        for (auto i : purne_list) {
            c.items.erase(i);
        }
    }

    void eliminate_hard_constraint(int index) {
        ConstraintItem &constraint = _constraints[index];
        int min_connected_constraints_size = 0x7fffffff;
        for (auto &pair : constraint.items) {
            if (abs(pair.second) > EPSILON) {
                if (constraint.eliminate_variable == -1 ||
                    _variables[pair.first].connected_constraints.size() <
                        min_connected_constraints_size) {
                    constraint.eliminate_variable = pair.first;
                    min_connected_constraints_size =
                        _variables[pair.first].connected_constraints.size();
                }
            }
        }
        if (constraint.eliminate_variable == -1) {
            return;
        }
        VariableItem &var = _variables[constraint.eliminate_variable];
        if (var.connected_constraints.size() > 10) {
            constraint.eliminate_variable = -1;
            return;
        }

        _eliminate_order.push_back(index);
        // Eliminate the variable.
        var.eliminated = true;
        number_t elim_weight = constraint.items[constraint.eliminate_variable];
        // Look through all constraints with this variable
        for (auto cindex : var.connected_constraints) {
            if (cindex == index)
                continue;
            ConstraintItem &c = _constraints[cindex];
            auto iter = c.items.find(var.index);
            if (iter != c.items.end()) {
                int original_size = c.ref.value();
                c.items.erase(iter);
                number_t weight = iter->second / elim_weight;
                c.bias -= weight * constraint.bias;
                for (auto it : constraint.items) {
                    if (it.first != constraint.eliminate_variable) {
                        auto itt = c.items.find(it.first);
                        if (itt != c.items.end()) {
                            itt->second -= weight * it.second;
                        } else {
                            c.items.insert(std::make_pair(it.first, -weight * it.second));
                            _variables[it.first].connected_constraints.insert(cindex);
                        }
                    }
                }
                purne_items(c, cindex);
                int new_size = c.ref.value();
                if (c.strength == kSolverStrength_HARD && c.eliminate_variable == -1) {
                    if (original_size != new_size) {
                        if (original_size <= MAX_ELIMINATE_SIZE) {
                            _constraint_queue[original_size].erase(c.heap_iterator);
                        }
                        if (new_size <= MAX_ELIMINATE_SIZE) {
                            _constraint_queue[new_size].push_front(c.ref);
                            c.heap_iterator = _constraint_queue[new_size].begin();
                        }
                    }
                }
            }
        }
    }

    void eliminate() {
        int total = 0;
        while (true) {
            bool found = false;
            for (int i = 0; i <= MAX_ELIMINATE_SIZE; ++i) {
                if (!_constraint_queue[i].empty()) {
                    auto ref = _constraint_queue[i].back();
                    _constraint_queue[i].pop_back();
                    eliminate_hard_constraint(ref.index);
                    found = true;
                    break;
                }
            }
            if (!found)
                break;
        }
    }

    void solve() {
        int variable_count = 0;
        for (auto &it : _variables) {
            if (!it.eliminated) {
                it.location = variable_count;
                variable_count++;
            } else {
                it.location = -1;
            }
        }

        // Construct A
        std::vector<Triplet<double>> tripletList;
        int constraint_count = 0;
        for (auto &it : _constraints) {
            if (it.eliminate_variable != -1)
                continue;
            constraint_count++;
        }
        tripletList.reserve(constraint_count * 5);
        SparseMatrix<double> A(constraint_count, variable_count);
        VectorXd x(variable_count);
        VectorXd x0(variable_count);
        VectorXd b(constraint_count);

        int index = 0;
        for (auto &it : _constraints) {
            if (it.eliminate_variable != -1)
                continue;
            number_t scaler = strength_to_scale(it.strength);
            for (auto item : it.items) {
                int variable_index = item.first;
                number_t weight = item.second;
                tripletList.push_back(
                    Triplet<double>(index, _variables[item.first].location, weight * scaler));
            }
            b[index] = -it.bias * scaler;
            index++;
        }

        A.setFromTriplets(tripletList.begin(), tripletList.end());

        for (auto &it : _variables) {
            if (!it.eliminated) {
                x0[it.location] = it.variable_ptr->value;
            }
        }

        Eigen::Index iters = 2 * A.cols();
        double tol_error = NumTraits<number_t>::epsilon();

        if (_attributes.find(kSolverAttribute_TOLERANCE) != _attributes.end()) {
            tol_error = _attributes[kSolverAttribute_TOLERANCE].float_value;
        }
        if (_attributes.find(kSolverAttribute_MAX_ITERATIONS) != _attributes.end()) {
            iters = _attributes[kSolverAttribute_MAX_ITERATIONS].int_value;
        }

        x = x0;
        if (constraint_count > 0) {
            lscg_solve_main(A, b, x, x0,
                            _attributes[kSolverAttribute_REGULARIZER_WEIGHT].float_value, iters,
                            tol_error);
            _attributes[kSolverAttribute_NUM_ITERATIONS].int_value = iters;
            _attributes[kSolverAttribute_ERROR].float_value = tol_error;
        } else {
            _attributes[kSolverAttribute_NUM_ITERATIONS].int_value = 0;
            _attributes[kSolverAttribute_ERROR].float_value = 0;
        }

        _attributes[kSolverAttribute_NUM_VARIABLES].int_value = variable_count;
        _attributes[kSolverAttribute_NUM_CONSTRAINTS].int_value = constraint_count;

        for (auto &it : _variables) {
            if (!it.eliminated) {
                it.variable_ptr->value = x[it.location];
            }
        }
        for (int i = (int)_eliminate_order.size() - 1; i >= 0; --i) {
            auto &c = _constraints[_eliminate_order[i]];
            number_t other = c.bias;
            number_t scale = c.items[c.eliminate_variable];
            for (auto &it : c.items) {
                if (it.first != c.eliminate_variable) {
                    other += it.second * _variables[it.first].variable_ptr->value;
                }
            }
            _variables[c.eliminate_variable].variable_ptr->value = -other / scale;
        }
    }

    void solveLagrange() {
        // Construct A, B
        int constraint_count = 0;
        int soft_constraint_count = 0;
        int hard_constraint_count = 0;
        int soft_triplet_count = 0;
        int hard_triplet_count = 0;
        for (auto &it : _constraints) {
            if (it.eliminate_variable != -1 || it.items.size() == 0)
                continue;
            constraint_count++;
            if (it.strength == kSolverStrength_HARD) {
                hard_constraint_count++;
                hard_triplet_count += it.items.size();
            } else {
                soft_constraint_count++;
                soft_triplet_count += it.items.size();
            }
        }
        if (soft_constraint_count == 0 || hard_constraint_count == 0) {
            solve();
            return;
        }

        int variable_count = 0;
        for (auto &it : _variables) {
            if (!it.eliminated) {
                it.location = variable_count;
                variable_count++;
            } else {
                it.location = -1;
            }
        }

        SparseMatrix<double> matA(soft_constraint_count, variable_count);
        VectorXd rhsB(soft_constraint_count);
        SparseMatrix<double> matC(hard_constraint_count, variable_count);
        VectorXd rhsD(hard_constraint_count);
        VectorXd x(variable_count);
        VectorXd x0(variable_count);

        std::vector<Triplet<double>> tripletList;
        {
            tripletList.reserve(soft_triplet_count);
            int index = 0;
            for (auto &it : _constraints) {
                if (it.eliminate_variable != -1 || it.items.size() == 0) {
                    continue;
                }
                if (it.strength == kSolverStrength_HARD) {
                    continue;
                }
                number_t scaler = strength_to_scale(it.strength);
                for (auto item : it.items) {
                    int variable_index = item.first;
                    number_t weight = item.second;
                    tripletList.push_back(
                        Triplet<double>(index, _variables[item.first].location, weight * scaler));
                }
                rhsB[index] = -it.bias * scaler;
                index++;
            }
            matA.setFromTriplets(tripletList.begin(), tripletList.end());
        }
        {
            tripletList.clear();
            tripletList.reserve(hard_triplet_count);
            int index = 0;
            for (auto &it : _constraints) {
                if (it.eliminate_variable != -1 || it.items.size() == 0) {
                    continue;
                }
                if (it.strength != kSolverStrength_HARD) {
                    continue;
                }
                number_t scaler = strength_to_scale(it.strength);
                for (auto item : it.items) {
                    int variable_index = item.first;
                    number_t weight = item.second;
                    tripletList.push_back(
                        Triplet<double>(index, _variables[item.first].location, weight * scaler));
                }
                rhsD[index] = -it.bias * scaler;
                index++;
            }
            matC.setFromTriplets(tripletList.begin(), tripletList.end());
        }

        for (auto &it : _variables) {
            if (!it.eliminated) {
                x0[it.location] = it.variable_ptr->value;
            }
        }

        Index num_iterations = (variable_count + hard_constraint_count) * 2;
        double tol_error = NumTraits<number_t>::epsilon();

        if (_attributes.find(kSolverAttribute_TOLERANCE) != _attributes.end()) {
            tol_error = _attributes[kSolverAttribute_TOLERANCE].float_value;
        }
        if (_attributes.find(kSolverAttribute_MAX_ITERATIONS) != _attributes.end()) {
            num_iterations = _attributes[kSolverAttribute_MAX_ITERATIONS].int_value;
        }

        constrained_least_square_conjugate_gradient(
            matA, matC, rhsB, rhsD, x0, CLSCGIdentityPreconditioner(), num_iterations, tol_error);

        x = x0;

        _attributes[kSolverAttribute_NUM_VARIABLES].int_value = variable_count;
        _attributes[kSolverAttribute_NUM_CONSTRAINTS].int_value = constraint_count;
        _attributes[kSolverAttribute_NUM_ITERATIONS].int_value = num_iterations;
        _attributes[kSolverAttribute_ERROR].float_value = tol_error;

        for (auto &it : _variables) {
            if (!it.eliminated) {
                it.variable_ptr->value = x[it.location];
            }
        }
        for (int i = (int)_eliminate_order.size() - 1; i >= 0; --i) {
            auto &c = _constraints[_eliminate_order[i]];
            number_t other = c.bias;
            number_t scale = c.items[c.eliminate_variable];
            for (auto &it : c.items) {
                if (it.first != c.eliminate_variable) {
                    other += it.second * _variables[it.first].variable_ptr->value;
                }
            }
            _variables[c.eliminate_variable].variable_ptr->value = -other / scale;
        }
    }
};

void SolverImpl::solveWithReduce() {
    Snapshot snapshot(_attributes);
    snapshot.initialize(_variables, _constraints);
    snapshot.eliminate();
    if (_attributes[kSolverAttribute_FLAGS].int_value & kSolverFlag_LAGRANGE) {
        snapshot.solveLagrange();
    } else {
        snapshot.solve();
    }
}
