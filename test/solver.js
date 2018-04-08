let m = require("..");
let mlog = require("mocha-logger");

// Initialize the module before testing
before(() => m.initialize());

// Test if a == b with 1e-5 tolerance
function assert_equal(a, b, message) {
    if (Math.abs(a - b) > 1e-5) {
        throw new Error("Assertion failed (" + message + "): " + a + " != " + b);
    }
}

function problem_chain(N) {
    return {
        name: "chain" + N,
        setup: (solver) => {
            let VAR_S = N + 1;
            solver.addVariable(VAR_S, 1, true);
            for (let i = 0; i < N; i++) {
                solver.addVariable(i, 0, true);
                if (i >= 1) {
                    let c = solver.addConstraint(m.ConstraintSolver.STRENGTH_HARD, 0);
                    solver.addConstraintCoefficient(c, i - 1, 1);
                    solver.addConstraintCoefficient(c, i, -1);
                    solver.addConstraintCoefficient(c, VAR_S, 1);
                }
            }
            solver.addConstraint(m.ConstraintSolver.STRENGTH_HARD, -1, [0], [1]);
            solver.addConstraint(m.ConstraintSolver.STRENGTH_MEDIUM, -1000, [N - 1], [1]);
            solver.addConstraint(m.ConstraintSolver.STRENGTH_MEDIUM, -500, [N - 1], [1]);
        },
        verify: (solver) => {
            let t = (750 - 1) / (N - 1);
            for (let i = 0; i < N; i++) {
                assert_equal(solver.getValue(i), 1 + t * i, "at " + i);
            }
        }
    };
}

function problem_chain_soft(N) {
    return {
        name: "chainsoft" + N,
        setup: (solver) => {
            let VAR_S = N + 1;
            solver.addVariable(VAR_S, 1, true);
            for (let i = 0; i < N; i++) {
                solver.addVariable(i, 0, true);
                if (i >= 1) {
                    solver.addConstraint(m.ConstraintSolver.STRENGTH_HARD, 0, [i - 1, i, VAR_S], [1, -1, 1]);
                }
            }
            solver.addConstraint(m.ConstraintSolver.STRENGTH_STRONG, -1, [VAR_S], [1]);
            solver.addConstraint(m.ConstraintSolver.STRENGTH_STRONG, -100, [0], [1]);
            solver.addConstraint(m.ConstraintSolver.STRENGTH_STRONG, -200, [N - 1], [1]);
        },
        verify: (solver) => {
            let x = 100, y = 200, n = N;
            let a1 = (2 * x + y + 1 + x * n * n - 2 * x * n - n) / (n * n - 2 * n + 3);
            let s = (-n * x + x + (n - 1) * y + 2) / (n * n - 2 * n + 3);
            for (let i = 0; i < N; i++) {
                assert_equal(solver.getValue(i), a1 + i * s, "at " + i);
            }
        }
    };
}

function test_solver(flags, problem) {
    let flag_names = [];
    if (flags == 0) flag_names.push("DEFAULT");
    if (flags & m.ConstraintSolver.FLAG_REDUCE) flag_names.push("REDUCE");
    if (flags & m.ConstraintSolver.FLAG_LAGRANGE) flag_names.push("LAGRANGE");
    it(`Solver (${flag_names.join(",")}) ` + problem.name, () => {
        let solver = new m.ConstraintSolver();
        solver.flags = flags;
        problem.setup(solver);

        t0 = new Date().getTime();
        solver.solve(solver);
        t1 = new Date().getTime();

        solver.computeLoss();

        mlog.log(`time(ms) = ${t1 - t0}\titerations = ${solver.num_iterations}\thardLoss = ${solver.hardLoss.toFixed(2)}\tsoftLoss = ${solver.softLoss.toFixed(2)}`);

        problem.verify(solver);
        solver.destroy(solver);
    });
}

describe("Solver", () => {
    [
        0,
        m.ConstraintSolver.FLAG_REDUCE,
        m.ConstraintSolver.FLAG_LAGRANGE,
        m.ConstraintSolver.FLAG_REDUCE | m.ConstraintSolver.FLAG_LAGRANGE
    ].forEach((flag) => {
        test_solver(flag, problem_chain(10));
        test_solver(flag, problem_chain(100));
        test_solver(flag, problem_chain(1000));
    });
    [
        m.ConstraintSolver.FLAG_LAGRANGE,
        m.ConstraintSolver.FLAG_REDUCE | m.ConstraintSolver.FLAG_LAGRANGE
    ].forEach((flag) => {
        test_solver(flag, problem_chain_soft(10));
    });
});

