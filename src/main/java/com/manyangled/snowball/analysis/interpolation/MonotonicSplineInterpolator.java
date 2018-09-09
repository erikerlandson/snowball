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

import org.apache.commons.math3.exception.DimensionMismatchException;

import org.apache.commons.math3.analysis.interpolation.UnivariateInterpolator;
import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction;

import static com.manyangled.snowball.analysis.interpolation.MSISupport.fitMonotoneSpline;

/**
 * Interpolates data using a spline that is constrained to be monotonic non-decreasing.
 */
public class MonotonicSplineInterpolator implements UnivariateInterpolator {
    private int m = M_DEFAULT;
    private double lambda = LAMBDA_DEFAULT;
    private double[] w = null;
    private double xmin = Double.NaN;
    private double xmax = Double.NaN;
    private ArrayList<Double> constraintX = new ArrayList<Double>();
    private ArrayList<Double> constraintY = new ArrayList<Double>();
    private ArrayList<Double> gConstraintX = new ArrayList<Double>();
    private ArrayList<Double> gConstraintY = new ArrayList<Double>();

    /**
     * Given data (x1, y1), (x2, y2)..., fit an interpolating spline that is constrained to be monotonic.
     * 
     * @param x the x data x1, x2, ...
     * @param y the y data y1, y2, ...
     * @return a polynomial spline that interpolates the data, and is monotonic non-decreasing over its 
     * interpolation domain.
     * <p>
     * NOTE: the number of data provided must be &ge; (m + 3), where (m) is the number of 
     * spline intervals configured. See the setM method below.
     */
    public PolynomialSplineFunction interpolate(double x[], double y[]) {
        final int n = x.length;
        final int M = m + 3;
        if (n < M) throw new IllegalArgumentException(String.format("data length (%d) must be >= %d", n, M));
        if (y.length != n) throw new DimensionMismatchException(y.length, n);
        if (w == null) {
            w = new double[n];
            for (int j = 0; j < n; ++j) w[j] = 1.0;
        }
        if (w.length != n) throw new DimensionMismatchException(w.length, n);
        for (int j = 0; j < n; ++j)
            if (w[j] <= 0.0) throw new IllegalArgumentException("weights (w) must be > 0");
        if (Double.isNaN(xmin)) {
            double z = x[0];
            for (int j = 1; j < n; ++j) if (x[j] < z) z = x[j];
            xmin = z;
        }
        if (Double.isNaN(xmax)) {
            double z = x[0];
            for (int j = 1; j < n; ++j) if (x[j] > z) z = x[j];
            xmax = z;
        }
        if (xmax <= xmin) throw new IllegalArgumentException("xMin must be < xMax");
        int nC = constraintX.size();
        double[] xC = new double[nC];
        double[] yC = new double[nC];
        for (int j = 0; j < nC; ++j) {
            xC[j] = constraintX.get(j);
            yC[j] = constraintY.get(j);
        }
        nC = gConstraintX.size();
        double[] xgC = new double[nC];
        double[] ygC = new double[nC];
        for (int j = 0; j < nC; ++j) {
            xgC[j] = gConstraintX.get(j);
            ygC[j] = gConstraintY.get(j);
        }

        return fitMonotoneSpline(x, y, m, xmin, xmax, lambda, w, xC, yC, xgC, ygC);
    }

    /**
     * Set the number of piecewise polynomial intervals over the interpolation domain.
     * @param m the number of piecewise intervals.
     * <p>
     * NOTE: m is currently required to be &ge; 4, due to internal numeric considerations.
     * (m + 3) is also expected to be &le; the number of data points provided for interpolation.
     * for example, if m is set to 5, then at least 8 data points must be provided for interpolation.
     */
    public void setM(int m) {
        if (m < M_MINIMUM)
            throw new IllegalArgumentException(String.format("m must be >= %d", M_MINIMUM));
        this.m = m;
    }

    /**
     * Set the smoothing parameter for the spline fitting
     * @param lambda the smoothing parameter. lambda is &gt; 0. Larger values increase bias toward smoothing.
     * Defaults to 1.
     */
    public void setLambda(double lambda) {
        if (lambda <= 0.0)
            throw new IllegalArgumentException("lambda must be > 0");
        this.lambda = lambda;
    }

    /**
     * Set the interpolation domain for the fitting. Values outside this domain will be considered illegal.
     * @param xMin the lower bound of the domain. Defaults to minimum x data value.
     * @param xMax the upper bound of the domain. Defaults to maximum x data value.
     */
    public void setBounds(double xMin, double xMax) {
        if (xMax <= xMin)
            throw new IllegalArgumentException("xMin must be < xMax");
        xmin = xMin;
        xmax = xMax;
    }

    /**
     * Set the weights for data points. Higher weights at (x, y) increase bias toward fitting an interpolation
     * that passes close to (x, y).
     * @param w the weight vector. Must be same length as data. Defaults to [1, 1, 1... ].
     */
    public void setW(double... w) {
        for (int j = 0; j < w.length; ++j)
            if (w[j] <= 0.0)
                throw new IllegalArgumentException("elements of w must be > 0");
        this.w = w;
    }

    /**
     * Add a hard equality constraint that the interpolation s(x) = y.
     * @param x the x value of the constraint
     * @param y the value that the interpolation s(x) is constrained to equal.
     */
    public void addEqualityConstraint(double x, double y) {
        constraintX.add(x);
        constraintY.add(y);
    }

    /**
     * Add a hard equality constraint that the derivative of interpolation s'(x) = dydx.
     * @param x the x value of the constraint
     * @param dydx the value that the interpolation derivative s'(x) is constrained to equal.
     * <p>
     * NOTE: dydx must be &ge; 0, to be consistent with monotonic interpolation.
     */
    public void addGradientEqualityConstraint(double x, double dydx) {
        // If and when I support monotone decreasing I'll need to defer this check.
        if (dydx < 0.0)
            throw new IllegalArgumentException("dydx cannot be negative for monotone spline fitting");
        gConstraintX.add(x);
        gConstraintY.add(dydx);
    }

    /** The default value for smoothing parameter lambda */
    public static final double LAMBDA_DEFAULT = 1.0;

    /** The default value for number of piecewise polynomial intervals (m) */
    public static final int M_DEFAULT = 5;

    /** The minimum number of piecewise intervals */
    private static final int M_MINIMUM = 4;
}
