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

    /** N<sub>0,3</sub>(t) from Table 1 */
    public static double N03(double t) {
        double omt = 1.0 - t;
        return omt * omt * omt / 6.0;
    }

    /** N<sub>1,3</sub>(t) from Table 1 */
    public static double N13(double t) {
        double t2 = t * t;
        double t3 = t * t2;
        return (4.0 - (6.0 * t2) + (3.0 * t3)) / 6.0;
    }

    /** N<sub>2,3</sub>(t) from Table 1 */
    public static double N23(double t) {
        double t2 = t * t;
        double t3 = t * t2;
        return (1.0 + (3.0 * t) + (3.0 * t2) - (3.0 * t3)) / 6.0;
    }

    /** N<sub>3,3</sub>(t) from Table 1 */
    public static double N33(double t) {
        return t * t * t / 6.0;
    }

    /** basis B<sub>3</sub>(t) from Eq(3) */
    public static double B3(double t) {
        if (t < 0.0) return 0.0;
        if (t < 1.0) return N33(t);
        if (t < 2.0) return N23(t - 1.0);
        if (t < 3.0) return N13(t - 2.0);
        if (t < 4.0) return N03(t - 3.0);
        return 0.0;
    }
}
