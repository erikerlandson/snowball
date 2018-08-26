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

import org.apache.commons.math3.exception.DimensionMismatchException;

import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;

import org.apache.commons.math3.analysis.interpolation.UnivariateInterpolator;
import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction;

public class MonotonicSplineInterpolator implements UnivariateInterpolator {
    private int m = M_DEFAULT;
    private double lambda = LAMBDA_DEFAULT;
    private double[] w = null;
    private double umin = UNSET_DOUBLE;
    private double umax = UNSET_DOUBLE;

    public PolynomialSplineFunction interpolate(double x[], double y[]) {
        final int n = x.length;
        if (n < m) throw new IllegalArgumentException("data length (n) must be >= m");
        if (y.length != n) throw new DimensionMismatchException(y.length, n);
        if (w == null) {
            double[] w = new double[n];
            for (int j = 0; j < n; ++j) w[j] = 1.0;
        }
        if (w.length != n) throw new DimensionMismatchException(w.length, n);
        for (int j = 0; j < n; ++j)
            if (w[j] <= 0.0) throw new IllegalArgumentException("weights (w) must be > 0");
        if (umin == UNSET_DOUBLE) {
            double z = x[0];
            for (int j = 1; j < n; ++j) if (x[j] < z) z = x[j];
            umin = z;
        }
        if (umax == UNSET_DOUBLE) {
            double z = x[0];
            for (int j = 1; j < n; ++j) if (x[j] > z) z = x[j];
            umax = z;
        }
        if (umax <= umin) throw new IllegalArgumentException("xMin must be < xMax");
        final double alpha = (double)m / (umax - umin);
        final int M = m + 3;
        final double[] tk = new double[M];
        for (int j = -3; j < m; ++j) tk[3+j] = umin + ((double)j / alpha);
        return null;
    }

    public void setM(int m) {
        if (m < M_MINIMUM)
            throw new IllegalArgumentException("m was too small");
        this.m = m;
    }

    public void setLambda(double lambda) {
        if (lambda <= 0.0)
            throw new IllegalArgumentException("lambda must be > 0");
        this.lambda = lambda;
    }

    public void setBounds(double xMin, double xMax) {
        if (xMax <= xMin)
            throw new IllegalArgumentException("xMin must be < xMax");
        umin = xMin;
        umax = xMax;
    }

    public void setW(double... w) {
        for (int j = 0; j < w.length; ++j)
            if (w[j] <= 0.0)
                throw new IllegalArgumentException("elements of w must be > 0");
        this.w = w;
    }

    public static final double UNSET_DOUBLE = Double.NaN;

    public static final double LAMBDA_DEFAULT = 1.0;

    public static final int M_DEFAULT = 5;

    private static final int M_MINIMUM = 4;

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

    private static RealMatrix matrixBT(double alpha, double[] tk, double[] u) {
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

    private static double[][] matrixRinf(int M) {
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
}
