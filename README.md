lscg-solver
====

A linear least square solver using the [Conjugate Gradient](https://en.wikipedia.org/wiki/Conjugate_gradient_method) method. Uses the Eigen library, and compiled as WebAssembly for performance.

# Install

    npm install lscg-solver

# Usage

You may `require` the solver as a node(commonjs) module:

    let LSCGSolver = require("lscg-solver");

You may also use the solver directly as a UMD module:

    <script src="node_modules/lscg-solver/dist/lscg-solver-umd.js"></script>

Once you have the solver imported, call the `initialize` method before using any of it:

    LSCGSolver.initialize().then(() => {
        // Code that use the solver
    });

# Documentation

- See `index.d.ts` for a full declaration of the methods available.
- See `test` for examples on how to use the solver.

License
----

MIT License