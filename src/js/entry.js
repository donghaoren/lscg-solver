var internals;
var isInitialized = false;
var initializedHooks = [];

var useWASM = true;
if (typeof (WebAssembly) == "undefined") {
    useWASM = false;
}

if (useWASM) {
    internals = require("../../build/internals.js")();
    internals.onRuntimeInitialized = function () {
        isInitialized = true;
        for (let i = 0; i < initializedHooks.length; i++) {
            initializedHooks[i]();
        }
    };
} else {
    console.warn("WASMSolver: not using web assembly");
    internals = require("../../build/internals_nowasm.js")();
    isInitialized = true;
}

function initialize() {
    return new Promise(function (resolve, reject) {
        if (isInitialized) {
            resolve();
        } else {
            initializedHooks.push(function () {
                resolve();
            });
        }
    });
}

// Wrap the functions and constants from api.h
var memory_alloc = internals.cwrap("memory_alloc", "number", ["number"]);
var memory_free = internals.cwrap("memory_free", null, ["number"]);

var getUint8Array = function (ptr, length) {
    return new Uint8Array(internals.HEAPU8.buffer, ptr, length);
};
var getUint16Array = function (ptr, length) {
    return new Uint16Array(internals.HEAPU16.buffer, ptr, length);
};
var getUint32Array = function (ptr, length) {
    return new Uint32Array(internals.HEAPU32.buffer, ptr, length);
};
var getFloat64Array = function (ptr, length) {
    return new Float64Array(internals.HEAPF64.buffer, ptr, length);
};
var getFloat32Array = function (ptr, length) {
    return new Float32Array(internals.HEAPF32.buffer, ptr, length);
};

var kSolverStrength_HARD = 0;
var kSolverStrength_STRONG = 1;
var kSolverStrength_MEDIUM = 2;
var kSolverStrength_WEAK = 3;
var kSolverStrength_WEAKER = 4;
var kSolverStrength_DISABLED = 10;

var kSolverAttribute_NUM_ITERATIONS = 1;
var kSolverAttribute_TOLERANCE = 2;
var kSolverAttribute_FLAGS = 3;

var kSolverAttribute_NUM_VARIABLES = 10;
var kSolverAttribute_NUM_CONSTRAINTS = 11;
var kSolverAttribute_MAX_ITERATIONS = 12;
var kSolverAttribute_ERROR = 13;

var kSolverAttribute_HARD_LOSS = 20;
var kSolverAttribute_SOFT_LOSS = 21;

var kSolverFlag_DEFAULT = 0;
var kSolverFlag_REDUCE = 1 << 1;
var kSolverFlag_LAGRANGE = 1 << 2;

var solver_create = internals.cwrap("solver_create", "number", []);
var solver_destroy = internals.cwrap("solver_destroy", null, ["number"]);
var solver_add_variable = internals.cwrap("solver_add_variable", null, ["number", "number", "number", "boolean"]);
var solver_make_constant = internals.cwrap("solver_make_constant", null, ["number", "number"]);
var solver_add_constraint = internals.cwrap("solver_add_constraint", "number", ["number", "number", "number", "number", "number", "number"]);
var solver_add_constraint_coefficient = internals.cwrap("solver_add_constraint_coefficient", null, ["number", "number", "number", "number"]);
var solver_set_constraint_strength = internals.cwrap("solver_set_constraint_strength", null, ["number", "number", "number"]);
var solver_set_constraint_bias = internals.cwrap("solver_set_constraint_bias", null, ["number", "number", "number"]);
var solver_clear_constraint_coefficients = internals.cwrap("solver_clear_constraint_coefficients", null, ["number", "number"]);
var solver_set_values = internals.cwrap("solver_set_values", null, ["number", "number", "number", "number"]);
var solver_get_values = internals.cwrap("solver_get_values", null, ["number", "number", "number", "number"]);
var solver_set_value = internals.cwrap("solver_set_value", null, ["number", "number", "number"]);
var solver_get_value = internals.cwrap("solver_get_value", "number", ["number", "number"]);
var solver_solve = internals.cwrap("solver_solve", null, ["number"]);
var solver_compute_loss = internals.cwrap("solver_compute_loss", null, ["number"]);
var solver_set_attribute_i = internals.cwrap("solver_set_attribute_i", null, ["number", "number", "number"]);
var solver_set_attribute_f = internals.cwrap("solver_set_attribute_f", null, ["number", "number", "number"]);
var solver_get_attribute_i = internals.cwrap("solver_get_attribute_i", "number", ["number", "number"]);
var solver_get_attribute_f = internals.cwrap("solver_get_attribute_f", "number", ["number", "number"]);

var linalg_matrix_create = internals.cwrap("linalg_matrix_create", "number", []);
var linalg_matrix_init = internals.cwrap("linalg_matrix_init", null, ["number", "number", "number", "number"]);
var linalg_matrix_fill = internals.cwrap("linalg_matrix_fill", null, ["number", "number"]);
var linalg_matrix_destroy = internals.cwrap("linalg_matrix_destroy", null, ["number"]);

var linalg_matrix_data = internals.cwrap("linalg_matrix_data", null, ["number"]);
var linalg_matrix_size = internals.cwrap("linalg_matrix_size", "number", ["number"]);
var linalg_matrix_rows = internals.cwrap("linalg_matrix_rows", "number", ["number"]);
var linalg_matrix_cols = internals.cwrap("linalg_matrix_cols", "number", ["number"]);
var linalg_matrix_row_stride = internals.cwrap("linalg_matrix_row_stride", "number", ["number"]);
var linalg_matrix_col_stride = internals.cwrap("linalg_matrix_col_stride", "number", ["number"]);
var linalg_matrix_norm = internals.cwrap("linalg_matrix_norm", null, ["number"]);
var linalg_matrix_l1_norm = internals.cwrap("linalg_matrix_l1_norm", null, ["number"]);

var linalg_matrix_add = internals.cwrap("linalg_matrix_add", null, ["number", "number", "number"]);
var linalg_matrix_sub = internals.cwrap("linalg_matrix_sub", null, ["number", "number", "number"]);
var linalg_matrix_scale = internals.cwrap("linalg_matrix_scale", null, ["number", "number", "number"]);
var linalg_matrix_add_scale = internals.cwrap("linalg_matrix_add_scale", null, ["number", "number", "number", "number", "number"]);
var linalg_matrix_emul = internals.cwrap("linalg_matrix_emul", null, ["number", "number", "number"]);
var linalg_matrix_ediv = internals.cwrap("linalg_matrix_ediv", null, ["number", "number", "number"]);
var linalg_matrix_mmul = internals.cwrap("linalg_matrix_mmul", null, ["number", "number", "number"]);

var linalg_solve_linear_system = internals.cwrap("linalg_solve_linear_system", null, ["number", "number", "number", "number"]);

// Helper classes
var WASMBuffer = function () {
    this.byte_size = 64;
    this.pointer = memory_alloc(this.byte_size);
    this._setupArrays();
}
function powerOf8(sz) {
    if (sz % 8 == 0) return sz;
    return sz + (8 - sz % 8);
}
WASMBuffer.prototype._setupArrays = function () {
    this.U32 = getUint32Array(this.pointer, this.byte_size / 4);
    this.F64 = getFloat64Array(this.pointer, this.byte_size / 8);
}
WASMBuffer.prototype.destroy = function () {
    this.byte_size = 0;
    this.U32 = null;
    this.F64 = null;
    memory_free(this.pointer);
};
WASMBuffer.prototype.reserve = function (bytes) {
    if (this.byte_size >= bytes) return;
    else {
        memory_free(this.pointer);
        this.byte_size = powerOf8(bytes * 2);
        this.pointer = memory_alloc(this.byte_size);
        this._setupArrays();
    }
};

// ConstraintSolver class
var ConstraintSolver = function () {
    this.solver = solver_create();
    this.buffer = new WASMBuffer();
};

ConstraintSolver.STRENGTH_HARD = kSolverStrength_HARD;
ConstraintSolver.STRENGTH_STRONG = kSolverStrength_STRONG;
ConstraintSolver.STRENGTH_MEDIUM = kSolverStrength_MEDIUM;
ConstraintSolver.STRENGTH_WEAK = kSolverStrength_WEAK;
ConstraintSolver.STRENGTH_WEAKER = kSolverStrength_WEAKER;
ConstraintSolver.STRENGTH_DISABLED = kSolverStrength_DISABLED;
ConstraintSolver.FLAG_DEFAULT = kSolverFlag_DEFAULT;
ConstraintSolver.FLAG_REDUCE = kSolverFlag_REDUCE;
ConstraintSolver.FLAG_LAGRANGE = kSolverFlag_LAGRANGE;

ConstraintSolver.prototype.destroy = function () {
    solver_destroy(this.solver);
    this.buffer.destroy();
};
ConstraintSolver.prototype.addVariable = function (variable_name, value, edit) {
    solver_add_variable(this.solver, variable_name, value, edit);
};
ConstraintSolver.prototype.makeConstant = function (variable_name) {
    solver_make_constant(this.solver, variable_name);
};
ConstraintSolver.prototype.addConstraint = function (strength, bias, variable_names, weights) {
    if (variable_names != null && weights != null) {
        this.buffer.reserve(16 * weights.length);
        let a1 = this.buffer.U32;
        let a2 = this.buffer.F64;
        for (var i = 0; i < weights.length; i++) {
            a1[i] = variable_names[i];
            a2[i + weights.length] = weights[i];
        }
        return solver_add_constraint(this.solver, strength, bias, weights.length, this.buffer.pointer, this.buffer.pointer + (8 * weights.length));
    } else {
        return solver_add_constraint(this.solver, strength, bias, 0, 0, 0);
    }
};
ConstraintSolver.prototype.addConstraintCoefficient = function (constraint, variable_name, weight) {
    solver_add_constraint_coefficient(this.solver, constraint, variable_name, weight);
};
ConstraintSolver.prototype.setConstraintStrength = function (constraint, strength) {
    solver_set_constraint_strength(this.solver, constraint, strength);
};
ConstraintSolver.prototype.setConstraintBias = function (constraint, bias) {
    solver_set_constraint_bias(this.solver, constraint, bias);
};
ConstraintSolver.prototype.clearConstraintCoefficients = function (constraint) {
    solver_clear_constraint_coefficients(this.solver, constraint);
};
ConstraintSolver.prototype.solve = function () {
    solver_solve(this.solver);
};
ConstraintSolver.prototype.computeLoss = function () {
    solver_compute_loss(this.solver);
};
ConstraintSolver.prototype.setValue = function (variable_name, value) {
    solver_set_value(this.solver, variable_name, value);
};
ConstraintSolver.prototype.getValue = function (variable_name) {
    return solver_get_value(this.solver, variable_name);
};
function defineSolverProperty(name, string_name, type, ro) {
    desc = {}
    if (type == "i") {
        desc.get = function () { return solver_get_attribute_i(this.solver, name); };
        desc.set = function (value) { solver_set_attribute_i(this.solver, name, value); };
    }
    if (type == "f") {
        desc.get = function () { return solver_get_attribute_f(this.solver, name); };
        desc.set = function (value) { solver_set_attribute_f(this.solver, name, value); };
    }
    if (ro) {
        desc.set = function () { throw new Error("property " + string_name + " is readonly"); };
    }
    Object.defineProperty(ConstraintSolver.prototype, string_name, desc);
}
defineSolverProperty(kSolverAttribute_MAX_ITERATIONS, "maxIterations", "i");
defineSolverProperty(kSolverAttribute_TOLERANCE, "tolerance", "f");
defineSolverProperty(kSolverAttribute_FLAGS, "flags", "i");
defineSolverProperty(kSolverAttribute_NUM_CONSTRAINTS, "numConstraints", "i", true);
defineSolverProperty(kSolverAttribute_NUM_VARIABLES, "numVariables", "i", true);
defineSolverProperty(kSolverAttribute_NUM_ITERATIONS, "numIterations", "i", true);
defineSolverProperty(kSolverAttribute_ERROR, "error", "f", true);
defineSolverProperty(kSolverAttribute_HARD_LOSS, "hardLoss", "f", true);
defineSolverProperty(kSolverAttribute_SOFT_LOSS, "softLoss", "f", true);

function Matrix() {
    this.matrix = linalg_matrix_create();
}
Matrix.prototype.destroy = function () {
    linalg_matrix_destroy(this.matrix);
};
Matrix.prototype.init = function (rows, cols) {
    linalg_matrix_init(this.matrix, rows, cols, 0);
};
Matrix.prototype.fill = function (value) {
    linalg_matrix_fill(this.matrix, value);
};
Matrix.prototype.data = function () {
    var length = linalg_matrix_size(this.matrix);
    return getFloat64Array(linalg_matrix_data(this.matrix), length);
};
Matrix.prototype.norm = function () {
    return linalg_matrix_norm(this.matrix);
};
Matrix.prototype.l1Norm = function () {
    return linalg_matrix_l1_norm(this.matrix);
};
Matrix.Add = function (dest, a, b) {
    linalg_matrix_add(dest.matrix, a.matrix, b.matrix);
};
Matrix.Sub = function (dest, a, b) {
    linalg_matrix_sub(dest.matrix, a.matrix, b.matrix);
};
Matrix.Scale = function (dest, a, s) {
    linalg_matrix_scale(dest.matrix, a.matrix, s);
};
Matrix.AddScale = function (dest, a, sa, b, sb) {
    linalg_matrix_add_scale(dest.matrix, a.matrix, sa, b.matrix, sb);
};
Matrix.EMul = function (dest, a, b) {
    linalg_matrix_emul(dest.matrix, a.matrix, b.matrix);
};
Matrix.EDiv = function (dest, a, b) {
    linalg_matrix_ediv(dest.matrix, a.matrix, b.matrix);
};
Matrix.MMul = function (dest, a, b) {
    linalg_matrix_mmul(dest.matrix, a.matrix, b.matrix);
};
Object.defineProperty(Matrix.prototype, "rows", {
    get: function () { return linalg_matrix_rows(this.matrix); },
    set: function () { throw new Error("property rows is readonly"); }
});
Object.defineProperty(Matrix.prototype, "rowStride", {
    get: function () { return linalg_matrix_row_stride(this.matrix); },
    set: function () { throw new Error("property row_stride is readonly"); }
});
Object.defineProperty(Matrix.prototype, "cols", {
    get: function () { return linalg_matrix_cols(this.matrix); },
    set: function () { throw new Error("property cols is readonly"); }
});
Object.defineProperty(Matrix.prototype, "colStride", {
    get: function () { return linalg_matrix_col_stride(this.matrix); },
    set: function () { throw new Error("property col_stride is readonly"); }
});

Matrix.SolveLinearSystem = function (X, ker, A, B) {
    linalg_solve_linear_system(X.matrix, ker.matrix, A.matrix, B.matrix);
};

// Exports
exports.memoryAlloc = memory_alloc;
exports.memoryFree = memory_free;

exports.getUint8Array = getUint8Array;
exports.getUint16Array = getUint16Array;
exports.getUint32Array = getUint32Array;
exports.getFloat64Array = getFloat64Array;
exports.getFloat32Array = getFloat32Array;

exports.ConstraintSolver = ConstraintSolver;

exports.Matrix = Matrix;

exports.initialize = initialize;