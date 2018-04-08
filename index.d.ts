export type Pointer = number;
export type Constraint = Pointer;

export class ConstraintSolver {
    /** Hard constraint */
    static STRENGTH_HARD: number;
    /** Strong constraint */
    static STRENGTH_STRONG: number;
    /** Medium constraint */
    static STRENGTH_MEDIUM: number;
    /** Weak constraint */
    static STRENGTH_WEAK: number;
    /** Weaker constraint */
    static STRENGTH_WEAKER: number;

    /** Default flag (0) */
    static FLAG_DEFAULT: number;
    /** Set this to enable hard constraint reducing */
    static FLAG_REDUCE: number;
    /** Set this to enable lagrange method to enforce hard constraints */
    static FLAG_LAGRANGE: number;

    constructor();

    /** Destroy the solver. Must be called after use to release memory. */
    destroy(): void;

    /** Add a variable to the solver */
    addVariable(variable_name: number, value: number, edit: boolean): void;
    /** Mark a variable as constant */
    makeConstant(variable_name: number): void;
    /** Add a new constraint to the solver */
    addConstraint(strength: number, bias: number): Constraint;
    /** Add a new constraint to the solver */
    addConstraint(strength: number, bias: number, variable_names: number[], weights: number[]): Constraint;
    /** Add a coefficient to a given constraint */
    addConstraintCoefficient(constraint: Constraint, variable_name: number, weight: number): void;
    /** Set the value of a variable */
    setValue(variable_name: number, value: number): void;
    /** Get the value of a variable */
    getValue(variable_name: number): number;
    /** Solve the ocosntraints */
    solve(): void;
    /** Compute the loss of the current state */
    computeLoss(): void;


    /** Parameter: max iterations, default to 2 * number of variables */
    maxIterations: number;
    /** Parameter: max iterations, default to the platform's epsilon */
    tolerance: number;

    /** Parameter: flags, default to FLAG_DEFAULT (0) */
    flags: number;

    /** The number of variables (after elimination in MODE_REDUCE) */
    readonly numVariables: number;
    /** The number of constraints (after elimination in MODE_REDUCE) */
    readonly numConstraints: number;
    /** The number of iterations taken */
    readonly numIterations: number;
    /** The residue error */
    readonly error: number;
    /** The loss of hard constraints */
    readonly hardLoss: number;
    /** The loss of soft constraints */
    readonly softLoss: number;
}

/** Allocate memory */
export function memoryAlloc(byte_size: number): Pointer;
/** Free memory */
export function memoryFree(pointer: Pointer): void;

/** Access value from pointer */
export function getUint8Array(ptr: Pointer, length: number): Uint8Array;
export function getUint16Array(ptr: Pointer, length: number): Uint16Array;
export function getUint32Array(ptr: Pointer, length: number): Uint32Array;
export function getFloat32Array(ptr: Pointer, length: number): Float32Array;
export function getFloat64Array(ptr: Pointer, length: number): Float64Array;

/** Matrix class */
export class Matrix {
    constructor();
    /** Free the memory. Must be called to prevent memory leak */
    destroy(): void;

    /** Init with rows and cols, all zero */
    init(rows: number, cols: number): void;
    /** Fill with value */
    fill(value: number): void;
    /** Get the data array */
    data(): Float64Array;

    /** Number of rows */
    readonly rows: number;
    /** Row stride */
    readonly rowStride: number;
    /** Number of columns */
    readonly cols: number;
    /** Column stride */
    readonly colStride: number;

    /** Compute the norm */
    norm(): number;
    /** Compute the l1 norm */
    l1Norm(): number;

    /** dest = a + b */
    static Add(dest: Matrix, a: Matrix, b: Matrix): void;
    /** dest = a - b */
    static Sub(dest: Matrix, a: Matrix, b: Matrix): void;
    /** dest = a * b, element-wise */
    static EMul(dest: Matrix, a: Matrix, b: Matrix): void;
    /** dest = a / b, element-wise */
    static EDiv(dest: Matrix, a: Matrix, b: Matrix): void;
    /** dest = a * scale */
    static Scale(dest: Matrix, a: Matrix, scale: number): void;
    /** dest = a * scale_a + b * scale_b */
    static AddScale(dest: Matrix, a: Matrix, scale_a: number, b: Matrix, scale_b: number): void;

    /** dest = a * b, matrix multiplication */
    static MMul(dest: Matrix, a: Matrix, b: Matrix): void;

    /** Solve Ax = B, output X and the null kernel */
    static SolveLinearSystem(X: Matrix, ker: Matrix, A: Matrix, B: Matrix): void;
}

/** Initialize the library. Must be called before using. */
export function initialize(): Promise<void>;