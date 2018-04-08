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

describe("Matrix", () => {
    it("fill", () => {
        let v = new m.Matrix();
        v.init(10, 10);
        v.fill(2);
        assert_equal(v.norm(), Math.sqrt(4 * 100));
        assert_equal(v.l1Norm(), 2 * 100);
        let d = v.data();
        let sum_sqr = 0;
        for (let i = 0; i < 100; i++) {
            d[i] = i + 1;
            sum_sqr += (i + 1) * (i + 1);
        }
        assert_equal(v.norm(), Math.sqrt(sum_sqr));
        assert_equal(v.l1Norm(), 5050);
        v.destroy();
    });

    it("add,sub,scale,addScale", () => {
        let v1 = new m.Matrix();
        let v2 = new m.Matrix();
        let v3 = new m.Matrix();
        v1.init(10, 10); let d1 = v1.data();
        v2.init(10, 10); let d2 = v2.data();
        v3.init(10, 10); let d3 = v3.data();
        for (let i = 0; i < 100; i++) {
            d1[i] = Math.random();
            d2[i] = Math.random();
        }
        m.Matrix.Add(v3, v1, v2);
        for (let i = 0; i < 100; i++) {
            assert_equal(d3[i], d1[i] + d2[i]);
        }
        m.Matrix.Sub(v3, v1, v2);
        for (let i = 0; i < 100; i++) {
            assert_equal(d3[i], d1[i] - d2[i]);
        }
        m.Matrix.Scale(v3, v1, 2001);
        for (let i = 0; i < 100; i++) {
            assert_equal(d3[i], d1[i] * 2001);
        }
        m.Matrix.AddScale(v3, v1, 10, v2, 13);
        for (let i = 0; i < 100; i++) {
            assert_equal(d3[i], d1[i] * 10 + d2[i] * 13);
        }
        v1.destroy();
        v2.destroy();
        v3.destroy();
    });

    it("solve_linear_system", () => {
        let nA = 400;
        let nX = 1000;
        let nB = 300;
        let A = new m.Matrix();
        A.init(nA, nX);
        var d = A.data();
        for (let i = 0; i < d.length; i++) {
            d[i] = Math.random();
        }
        let B = new m.Matrix();
        B.init(nA, nB);
        var d = B.data();
        for (let i = 0; i < d.length; i++) {
            d[i] = Math.random();
        }
        let X = new m.Matrix();
        let ker = new m.Matrix();

        let t0 = new Date().getTime();
        m.Matrix.SolveLinearSystem(X, ker, A, B);
        let t1 = new Date().getTime();
        mlog.log(`time(ms) = ${t1 - t0}`);

        let AX = new m.Matrix();
        m.Matrix.MMul(AX, A, X);
        m.Matrix.Sub(B, AX, B);

        assert_equal(B.norm(), 0);

        A.destroy();
        B.destroy();
        ker.destroy();
        X.destroy();
        AX.destroy();
    });
});