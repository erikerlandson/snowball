/*
Copyright 2018 Erik Erlandson
Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at
http://www.apache.org/licenses/LICENSE-2.0
Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/

package com.manyangled.snowball.analysis.interpolation;

import java.util.ArrayList;
import static java.util.Arrays.binarySearch;

import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.DiagonalMatrix;

import org.apache.commons.math3.optim.OptimizationData;
import org.apache.commons.math3.optim.PointValuePair;
import org.apache.commons.math3.optim.InitialGuess;
import org.apache.commons.math3.optim.nonlinear.scalar.ObjectiveFunction;

import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction;
import org.apache.commons.math3.analysis.polynomials.PolynomialFunction;

import static com.manyangled.gibbous.optim.convex.ConvexOptimizer.feasiblePoint;
import com.manyangled.gibbous.optim.convex.BarrierOptimizer;
import com.manyangled.gibbous.optim.convex.QuadraticFunction;
import com.manyangled.gibbous.optim.convex.LinearInequalityConstraint;
import com.manyangled.gibbous.optim.convex.LinearEqualityConstraint;
import com.manyangled.gibbous.optim.convex.SVDSchurKKTSolver;

import com.manyangled.gibbous.optim.convex.QuadraticFunction;
import com.manyangled.gibbous.optim.convex.LinearInequalityConstraint;

class MSISupport {
    /** basis B<sub>3</sub>(t) from Eq(3) */
    public static double B3(double t) {
        if (t < 0.0) return 0.0;
        if (t < 1.0) {
            // N33(t)
            return t * t * t / 6.0;
        }
        if (t < 2.0) {
            // N23(t-1)
            t -= 1.0;
            double t2 = t * t;
            double t3 = t * t2;
            return (1.0 + (3.0 * t) + (3.0 * t2) - (3.0 * t3)) / 6.0;
        }
        if (t < 3.0) {
            // N13(t-2)
            t -= 2.0;
            double t2 = t * t;
            double t3 = t * t2;
            return (4.0 - (6.0 * t2) + (3.0 * t3)) / 6.0;
        }
        if (t < 4.0) {
            // N03(t-3)
            t = 4.0 - t;  // 1 - (t-3)
            return t * t * t / 6.0;
        }
        return 0.0;
    }

    public static RealMatrix matrixBT(double alpha, double[] tk, double[] u) {
        // tk is array of knot points, dimension M
        // u is vector of data "x" points; in same axis as knot points tk
        // note this returns B-transposed relative to Eq(12).
        int n = u.length;
        int M = tk.length;
        double[][] bt = new double[n][M];
        for (int j = 0; j < n; ++j) {
            for (int k = 0; k < M; ++k) {
                bt[j][k] = B3(alpha * (u[j] - tk[k]));
            }
        }
        return new Array2DRowRealMatrix(bt, false);
    }

    public static double[][] matrixRinf(int M) {
        // returns matrix Rinf Eq(26) for size MxM
        // unsure if this is hard minimum but makes logic safe for now by
        // guaranteeing no overlap of special sections
        if (M < 7) throw new IllegalArgumentException("R must have >= 7 rows");
        double[] band = { 1.0 / 6.0, 0.0, -9.0 / 6.0, 16.0 / 6.0, -9.0 / 6.0, 0.0, 1.0 / 6.0 };
        int b = band.length;
        double[][] rinf = new double[M][M];
        // first 3 rows: clip band on the left
        for (int j = 0; j < 3; ++j) {
            int d = 3 - j;
            for (int k = d; k < b; ++k) rinf[j][k - d] = band[k];
            for (int k = b - d; k < M; ++k) rinf[j][k] = 0.0;
        }
        // the middle part: all of band will fit
        for (int j = 3; j < M - 3; ++j) {
            int z = j - 3;
            for (int k = 0; k < z; ++k) rinf[j][k] = 0.0;
            for (int k = 0; k < b; ++k) rinf[j][z + k] = band[k];
            for (int k = z + b; k < M; ++k) rinf[j][k] = 0.0;
        }
        // last 3 rows: clip band on the right
        for (int d = 1; d <= 3; ++d) {
            int j = M - 3 + (d - 1);
            int z = M - (b - d);
            for (int k = 0; k < z; ++k) rinf[j][k] = 0.0;
            for (int k = 0; k < (b - d); ++k) rinf[j][z + k] = band[k];
        }
        return rinf;
    }

    public static RealMatrix lambdaQ(int M, double lambda, double alpha) {
        // From Eq(27), upper left corner only, delaying 1/6 factor
        double[][] rm = {
            { 14.0,  -6.0,   0.0 },
            { -6.0,   8.0,  -3.0 },
            { 0.0,  -3.0,   2.0 }
        };
        // From Eq(28), lower right corner only
        double[][] rp = {
            {  2.0,  -3.0,   0.0 },
            { -3.0,   8.0,  -6.0 },
            {  0.0,  -6.0,   14.0 }
        };
        double[][] lq = matrixRinf(M); // Rinf, from Eq(26)
        // From Eq(25): subtract rm and rp to get R
        for (int j = 0; j < 3; ++j) {
            for (int k = 0; k < 3; ++k) {
                lq[j][k] -= rm[j][k] / 6.0;
                lq[M - 3 + j][M - 3 + k] -= rp[j][k] / 6.0;
            }
        }
        // Eq(15) & Eq(22) - matrix lambda Q from R
        double f = lambda * alpha * alpha * alpha;
        for (int j = 0; j < M; ++j)
            for (int k = 0; k < M; ++k)
                lq[j][k] *= f;
        return new Array2DRowRealMatrix(lq, false);
    }

    public static QuadraticFunction quadraticObjective(
        double[] u,
        double[] d,
        double[] tk,
        double[] w,
        double lambda,
        double alpha) {
        // constructs G, g and r from Eq(14)
        // Warning: Eq(68) is presented as a generalization but its (tau).(g) term is
        // positive not negative, and it is missing the (r) term.
        RealMatrix BT = matrixBT(alpha, tk, u);
        RealMatrix B = BT.transpose();
        RealMatrix lambdaQ = lambdaQ(tk.length, lambda, alpha);
        RealMatrix W = new DiagonalMatrix(w, false);
        // Eq(15)
        RealMatrix G = lambdaQ.add(B.multiply(W).multiply(BT));
        // Eq(16); 'operate' method is "multiply by vector"
        // This adds a (-1) factor for the subtraction from Eq(14)
        double[] rt = W.operate(d);
        double[] g = B.operate(rt);
        for (int j = 0; j < g.length; ++j) g[j] *= -1.0;
        // Eq(17)
        double r = 0.0;
        for (int j = 0; j < d.length; ++j) r += rt[j] * d[j];
        return new QuadraticFunction(G.getData(), g, r);
    }

    public static LinearInequalityConstraint monotoneConstraints(int m, int M) {
        // negative of Eq(64), because barrier method expects constraints of
        // form Hx < 0 instead of Hx >= 0
        double[] Fm3 = { 1.0, 0.0, -1.0,  0.0 };
        double[] Fm2 = { 1.0, 2.0, -3.0,  0.0 };
        double[] Fm1 = { 0.0, 1.0,  0.0, -1.0 };
        double[] Fm0 = { 0.0, 3.0, -2.0, -1.0 };
        int f = Fm0.length;
        
        // apply constraints over the interpolation interval (starting at knot interval [t0, t1])
        // Eq(66), replicated for each interval starting at t0, t1, ... t(m-1)
        double[][] H = new double[4 * m][M];
        int j = 0;
        for (int z = 0; z < m; ++z) {
            for (int k = 0; k < z; ++k) H[j][k] = 0.0;
            for (int k = 0; k < f; ++k) H[j][z + k] = Fm3[k];
            for (int k = z + f; k < M; ++k) H[j][k] = 0.0;
            j += 1;
            for (int k = 0; k < z; ++k) H[j][k] = 0.0;
            for (int k = 0; k < f; ++k) H[j][z + k] = Fm2[k];
            for (int k = z + f; k < M; ++k) H[j][k] = 0.0;
            j += 1;
            for (int k = 0; k < z; ++k) H[j][k] = 0.0;
            for (int k = 0; k < f; ++k) H[j][z + k] = Fm1[k];
            for (int k = z + f; k < M; ++k) H[j][k] = 0.0;
            j += 1;
            for (int k = 0; k < z; ++k) H[j][k] = 0.0;
            for (int k = 0; k < f; ++k) H[j][z + k] = Fm0[k];
            for (int k = z + f; k < M; ++k) H[j][k] = 0.0;
            j += 1;            
        }

        double[] h = new double[4 * m];
        for (int k = 0; k < h.length; ++k) h[k] = 0.0;

        return new LinearInequalityConstraint(H, h);
    }

    public static PolynomialSplineFunction polynomialSplineFunction(double[] tau, double alpha, double xmin) {
        assert tau.length >= 4;
        int m = tau.length - 3;

        PolynomialFunction[] poly = new PolynomialFunction[m];
        double[] knot = new double[m + 1];
        for (int j = 0; j < m; ++j) {
            knot[j] = xmin + ((double)j / alpha);
            poly[j] = new PolynomialFunction(standardCoefficientsB3(3 + j, tau, alpha));
        }
        knot[m] = xmin + ((double)m / alpha);

        return new PolynomialSplineFunction(knot, poly);
    }

    // coefficients for knot interval [K[j], K[j+1]]
    // http://erikerlandson.github.io/blog/2018/09/02/putting-cubic-b-splines-into-standard-polynomial-form/
    public static double[] standardCoefficientsB3(int j, double[] tau, double alpha) {
        double[] c = new double[4];
        double a = 1.0 / 6.0;
        c[0] = a * (tau[j-1] + (4.0 * tau[j-2]) + tau[j-3]);
        a *= alpha;
        c[1] = a * ((3.0 * tau[j-1]) - (3.0 * tau[j-3]));
        a *= alpha;
        c[2] = a * ((3.0 * tau[j-1]) - (6.0 * tau[j-2]) + (3.0 * tau[j-3]));
        a *= alpha;
        c[3] = a * (tau[j] - (3.0 * tau[j-1]) + (3.0 * tau[j-2]) - tau[j-3]);
        return c;
    }

    public static int queryKj(double x, double[] K) {
        int j = binarySearch(K, x);
        if (j < 0) j = -(j + 2);
        if ((j < 0) || (j >= K.length))
            throw new IndexOutOfBoundsException("binary search landed outside of array");
        assert x >= K[j];
        if ((j + 1) < K.length) assert x < K[j + 1];
        return j;
    }

    public static LinearEqualityConstraint linearEqualityConstraint(
        double[] K,
        double alpha,
        double xmin,
        double xmax,
        double[] xC,
        double[] yC,
        double[] xgC,
        double[] ygC) {
        int n = xC.length + xgC.length;
        int M = K.length;
        double[][] A = new double[n][M];
        double[] b = new double[n];
        int jA = 0;
        for (int j = 0; j < xC.length; ++j) {
            double x = xC[j];
            if ((x < xmin) || (x > xmax))
                throw new IllegalArgumentException("equality constraint declared outside the interpolation domain");
            int q = queryKj(x, K);
            assert q >= 3;
            double t = alpha * (x - K[q]);
            b[jA] = yC[j];
            int z = q - 3;
            for (int k = 0; k < z; ++k) A[jA][k] = 0.0;
            // http://erikerlandson.github.io/blog/2018/09/08/equality-constraints-for-cubic-b-splines/
            A[jA][q - 3] = (1.0 - (3.0 * t) + (3.0 * t * t) - (t * t * t)) / 6.0;
            A[jA][q - 2] = (4.0 - (6.0 * t * t) + (3.0 * t * t * t)) / 6.0;
            A[jA][q - 1] = (1.0 + (3.0 * t) + (3.0 * t * t) - (3.0 * t * t * t)) / 6.0;
            A[jA][q]     = (t * t * t) / 6.0;
            for (int k = q + 1; k < M; ++k) A[jA][k] = 0.0;
            ++jA;
        }
        for (int j = 0; j < xgC.length; ++j) {
            double x = xgC[j];
            if ((x < xmin) || (x > xmax))
                throw new IllegalArgumentException("equality constraint declared outside the interpolation domain");
            int q = queryKj(x, K);
            assert q >= 3;
            double t = alpha * (x - K[q]);
            b[jA] = ygC[j];
            int z = q - 3;
            for (int k = 0; k < z; ++k) A[jA][k] = 0.0;
            // http://erikerlandson.github.io/blog/2018/09/08/equality-constraints-for-cubic-b-splines/
            A[jA][q - 3] = alpha * (-3.0 + (6.0 * t) - (3.0 * t * t)) / 6.0;
            A[jA][q - 2] = alpha * ((-12.0 * t) + (9.0 * t * t)) / 6.0;
            A[jA][q - 1] = alpha * (3.0 + (6.0 * t) - (9.0 * t * t)) / 6.0;
            A[jA][q]     = alpha * (3.0 * t * t) / 6.0;
            for (int k = q + 1; k < M; ++k) A[jA][k] = 0.0;
            ++jA;
        }
        return new LinearEqualityConstraint(new Array2DRowRealMatrix(A, false), new ArrayRealVector(b, false));
    }

    public static LinearInequalityConstraint linearInequalityConstraint(
        double[] K,
        double alpha,
        double xmin,
        double xmax,
        double[] xltC,
        double[] yltC,
        double[] ltCF) {
        int n = xltC.length;
        int M = K.length;
        double[][] A = new double[n][M];
        double[] b = new double[n];
        for (int j = 0; j < n; ++j) {
            double x = xltC[j];
            if ((x < xmin) || (x > xmax))
                throw new IllegalArgumentException("equality constraint declared outside the interpolation domain");
            int q = queryKj(x, K);
            assert q >= 3;
            double t = alpha * (x - K[q]);
            b[j] = ltCF[j] * yltC[j];
            int z = q - 3;
            for (int k = 0; k < z; ++k) A[j][k] = 0.0;
            A[j][q - 3] = ltCF[j] * (1.0 - (3.0 * t) + (3.0 * t * t) - (t * t * t)) / 6.0;
            A[j][q - 2] = ltCF[j] * (4.0 - (6.0 * t * t) + (3.0 * t * t * t)) / 6.0;
            A[j][q - 1] = ltCF[j] * (1.0 + (3.0 * t) + (3.0 * t * t) - (3.0 * t * t * t)) / 6.0;
            A[j][q]     = ltCF[j] * (t * t * t) / 6.0;
            for (int k = q + 1; k < M; ++k) A[j][k] = 0.0;
        }
        return new LinearInequalityConstraint(new Array2DRowRealMatrix(A, false), new ArrayRealVector(b, false));
    }

    public static PolynomialSplineFunction fitMonotoneSpline(
        double[] x,
        double[] y,
        int m,
        double xmin,
        double xmax,
        double lambda,
        double[] w,
        double[] xC,
        double[] yC,
        double[] xgC,
        double[] ygC,
        double[] xltC,
        double[] yltC,
        double[] ltCF,
        ArrayList<OptimizationData> fitOpts)
    {
        final double alpha = (double)m / (xmax - xmin);
        final int M = m + 3;
        final double[] K = new double[M];
        for (int j = -3; j < m; ++j) K[3+j] = xmin + ((double)j / alpha);

        ArrayList<OptimizationData> optArgs = new ArrayList<OptimizationData>();

        // Equality constraints seem to be producing positive-semi-definite matrix
        // and Cholesky solver isn't kidding about wanting strict positive definite.
        // I'm guessing that the hyperplanar constraints are reducing the rank.
        // SVD was born to solvev matrices of less than full rank, and seems to be working.
        optArgs.add(new SVDSchurKKTSolver());

        // include user supplied options
        // Any options I add after this will override any user settings.
        optArgs.addAll(fitOpts);

        if ((xC.length + xgC.length) > 0) {
            LinearEqualityConstraint eqc = linearEqualityConstraint(K, alpha, xmin, xmax, xC, yC, xgC, ygC);
            optArgs.add(eqc);
        }

        if (xltC.length > 0) {
            LinearInequalityConstraint iqc = linearInequalityConstraint(K, alpha, xmin, xmax, xltC, yltC, ltCF);
            optArgs.add(iqc);
        }

        LinearInequalityConstraint monotone = monotoneConstraints(m, M);
        optArgs.add(monotone);

        PointValuePair fpvp = feasiblePoint(optArgs.toArray(new OptimizationData[0]));
        if (fpvp.getSecond() >= 0.0)
            throw new RuntimeException("Unable to find an initial point in the feasible region");
        double[] ig = fpvp.getFirst();

        QuadraticFunction qf = quadraticObjective(x, y, K, w, lambda, alpha);
        optArgs.add(new ObjectiveFunction(qf));
        optArgs.add(new InitialGuess(ig));

        PointValuePair pvp = (new BarrierOptimizer()).optimize(optArgs.toArray(new OptimizationData[0]));

        double[] tau = pvp.getFirst();
        return polynomialSplineFunction(tau, alpha, xmin);
    }
}
