lscg-solver
====

A linear least square solver using the [Conjugate Gradient](https://en.wikipedia.org/wiki/Conjugate_gradient_method) method.

- Solve sparse linear least square problems.
- Solve sparse constrained linear least square problems using a lagrange method.
- Decompose dense matrices with full pivot LU to obtain solution and kernel (null space).
- Uses the Eigen library for sparse matrix computation and the conjugate gradient method.
- Compiled into WebAssembly for performance.

## Install

```bash
npm install lscg-solver
```

## Usage

You may `require` the solver as a node(commonjs) module:

```javascript
let LSCGSolver = require("lscg-solver");
```

You may also use the solver directly as a UMD module:

```html
<script src="node_modules/lscg-solver/dist/lscg-solver-umd.js"></script>
```

Once you have the solver imported, call the `initialize` method before using any of it:

```javascript
LSCGSolver.initialize().then(() => {
    // Code that uses the solver
});
```

## Documentation

- See `index.d.ts` for a full declaration of the methods available.
- See `test` for examples on how to use the solver.

## Development

To build the solver, you'll need to first install [Emscripten](http://kripken.github.io/emscripten-site/). Make sure the `emcc` command is available in your `PATH`.

The build process was tested on macOS 10.13.4 with Emscripten 1.37.35.

Once the dependencies has been installed, simply clone the repository and run:

```bash
# Get the submodules (Eigen)
git submodule init
git submodule update

# Install required node modules
npm install

# Build
npm run build

# Test
npm run test
```

## License

MIT License

lscg-solver uses the Eigen library, which can be obtained from <http://eigen.tuxfamily.org/index.php?title=Main_Page>.
