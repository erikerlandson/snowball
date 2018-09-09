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

    public PolynomialSplineFunction interpolate(double x[], double y[]) {
        final int n = x.length;
        if (n < m) throw new IllegalArgumentException("data length (n) must be >= m");
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
        xmin = xMin;
        xmax = xMax;
    }

    public void setW(double... w) {
        for (int j = 0; j < w.length; ++j)
            if (w[j] <= 0.0)
                throw new IllegalArgumentException("elements of w must be > 0");
        this.w = w;
    }

    public void addEqualityConstraint(double x, double y) {
        constraintX.add(x);
        constraintY.add(y);
    }

    public void addGradientEqualityConstraint(double x, double dydx) {
        // If and when I support monotone decreasing I'll need to defer this check.
        if (dydx < 0.0)
            throw new IllegalArgumentException("dydx cannot be negative for monotone spline fitting");
        gConstraintX.add(x);
        gConstraintY.add(dydx);
    }

    public static final double LAMBDA_DEFAULT = 1.0;

    public static final int M_DEFAULT = 5;

    private static final int M_MINIMUM = 4;
}
