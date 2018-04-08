let multirun = require("multirun");
let fs = require("fs");

let files = [
    "src/cpp/api.cpp",
    "src/cpp/linalg.cpp",
    "src/cpp/solver.cpp",
    "src/cpp/impl_solver.cpp"
];

let flags = [
    ...files,
    "-std=c++11",
    "-s", "MODULARIZE=1",
    "-s", "NO_FILESYSTEM=1",
    "-s", `EXTRA_EXPORTED_RUNTIME_METHODS=["cwrap"]`,
    "-s", "ALLOW_MEMORY_GROWTH=1",
    "-s", "SINGLE_FILE=1",
    "-O3"
];


let COMMANDS = {};

function makeSingleCommand(command, args) {
    return command + " " + args.map(x => x.replace(/\"/g, "\\\"").replace(/\'/g, "\\\'").replace(/ /g, "\\ ")).join(" ");
}

COMMANDS["internals"] = makeSingleCommand("emcc", [...flags, "-o", "build/internals.js", "-s", "WASM=1"]);

COMMANDS["internals_nowasm"] = makeSingleCommand("emcc", [...flags, "-o", "build/internals_nowasm.js"])

COMMANDS["postprocess"] = () => {
    return new Promise((resolve, reject) => {
        let i1 = fs.readFileSync("build/internals.js", "utf-8");
        let i2 = fs.readFileSync("build/internals_nowasm.js", "utf-8");
        i1 = i1.replace(/require\(\".*?\"\)/g, "null");
        i2 = i2.replace(/require\(\".*?\"\)/g, "null");
        fs.writeFileSync("build/internals.js", i1, "utf-8");
        fs.writeFileSync("build/internals_nowasm.js", i2, "utf-8");
        resolve();
    });
};

COMMANDS["webpack"] = "webpack";

let sequence = process.argv.slice(2);
if (sequence.length == 0) {
    sequence = ["internals", "internals_nowasm", "postprocess", "webpack"];
}
multirun.runCommands(COMMANDS, sequence, "Build");