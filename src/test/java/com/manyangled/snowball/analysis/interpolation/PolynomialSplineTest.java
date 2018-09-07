package com.manyangled.snowball.analysis.interpolation;

import org.junit.Test;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertArrayEquals;

import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction;
import org.apache.commons.math3.analysis.polynomials.PolynomialFunction;

import static com.manyangled.snowball.analysis.interpolation.MSISupport.*;

public class PolynomialSplineTest {
    public static double refSpline(double x, double alpha, double[] tau, double[] K) {
        assert tau.length == K.length;
        int m = tau.length;
        double s = 0.0;
        for (int j = 0; j < m; ++j) {
            s += tau[j] * B3(alpha * (x - K[j]));
        }
        return s;
    }

    @Test
    public void test1() {
        double[] tau = { 1.0, 2.0, 3.0, 4.0 };
        int m = tau.length - 3;
        double xmin = 10;
        double xmax = 20;
        double alpha = (double)m / (xmax - xmin);
        double[] K = new double[tau.length];
        for (int j = -3; j < m; ++j) {
            K[j + 3] = xmin + ((double)j / alpha);
        }
        PolynomialSplineFunction spline = polynomialSplineFunction(tau, alpha, xmin);
        for (double x = xmin; x <= xmax; x += (xmax - xmin)/100.0) {
            double y = spline.value(x);
            double yRef = refSpline(x, alpha, tau, K);
            assertEquals(y, yRef, 1e-9);
        }
    }

    @Test
    public void test2() {
        double[] tau = { 1.0, 2.0, 4.0, 8.0, 16.0, 32.0 };
        int m = tau.length - 3;
        double xmin = -10;
        double xmax = 10;
        double alpha = (double)m / (xmax - xmin);
        double[] K = new double[tau.length];
        for (int j = -3; j < m; ++j) {
            K[j + 3] = xmin + ((double)j / alpha);
        }
        PolynomialSplineFunction spline = polynomialSplineFunction(tau, alpha, xmin);
        for (double x = xmin; x <= xmax; x += (xmax - xmin)/1000.0) {
            double y = spline.value(x);
            double yRef = refSpline(x, alpha, tau, K);
            assertEquals(y, yRef, 1e-9);
        }
    }
}
