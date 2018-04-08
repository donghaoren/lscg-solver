module.exports = {
    entry: "./src/js/entry.js",
    output: {
        filename: "lscg-solver-umd.js",
        path: __dirname + "/dist/",
        // Export the app as a global variable "Charticulator"
        libraryTarget: "umd",
        library: "LSCGSolver"
    }
};
