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

import org.junit.Test;
import org.junit.Rule;
import org.junit.rules.ExpectedException;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertArrayEquals;

import static org.hamcrest.MatcherAssert.assertThat;
import static org.hamcrest.number.OrderingComparison.*;
import static org.hamcrest.number.IsCloseTo.closeTo;

import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction;
import org.apache.commons.math3.optim.MaxIter;
import org.apache.commons.math3.optim.OptimizationData;

public class MonotonicSplineTest {
    public static void testMonotone(PolynomialSplineFunction s) {
        double[] K = s.getKnots();
        double xmin = K[0];
        double xmax = K[K.length - 1];
        assertThat(xmin, lessThan(xmax));
        PolynomialSplineFunction ds = s.polynomialSplineDerivative();
        double dx = (xmax - xmin) * 1e-5;
        double xprv = xmin;
        for (double x = xmin; x <= xmax; x += dx) {
            // directly test that s(x) is non-decreasing
            assertThat(s.value(x), greaterThanOrEqualTo(s.value(xprv)));
            // test that the derivative of s(x) is always >= 0
            assertThat(ds.value(x), greaterThanOrEqualTo(0.0));
            xprv = x;
        }
    }

    @Rule
    public ExpectedException thrown = ExpectedException.none();

    @Test
    public void test1() {
        double[] x = { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0 };
        double[] y = { 0.0, 0.05, 0.02, 0.3, 0.5, 0.7, 0.99, 0.95, 1.0 };
        MonotonicSplineInterpolator interpolator = new MonotonicSplineInterpolator();
        PolynomialSplineFunction s = interpolator.interpolate(x, y);
        testMonotone(s);
    }

    @Test
    public void test2() {
        // These data turned up a subtle bug in gibbous-0.1.0 that is fixed in gibbous-0.1.1
        double[] x = { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0 };
        double[] y = { 0.0, 0.15, 0.05, 0.3, 0.5, 0.7, 0.95, 0.98, 1.0 };
        MonotonicSplineInterpolator interpolator = new MonotonicSplineInterpolator();
        PolynomialSplineFunction s = interpolator.interpolate(x, y);
        testMonotone(s);
    }

    @Test
    public void test3() {
        // These data turned up a subtle bug in gibbous-0.1.0 that is fixed in gibbous-0.1.1
        double[] x = { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0 };
        double[] y = { 9.0, 8.0, 7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0 };
        MonotonicSplineInterpolator interpolator = new MonotonicSplineInterpolator();
        PolynomialSplineFunction s = interpolator.interpolate(x, y);
        testMonotone(s);
    }

    @Test
    public void testEqualityConstraint1() {
        double[] x = { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0 };
        double[] y = { 0.0, 0.2, 0.1, 0.4, 0.5, 0.6, 0.9, 0.8, 1.0 };
        MonotonicSplineInterpolator interpolator = new MonotonicSplineInterpolator();
        interpolator.addEqualityConstraint(5.0, 0.7);
        PolynomialSplineFunction s = interpolator.interpolate(x, y);
        testMonotone(s);
        assertThat(s.value(5.0), closeTo(0.7, 1e-9));
    }

    @Test
    public void testEqualityConstraint2() {
        double[] x = { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0 };
        double[] y = { 0.0, 0.2, 0.1, 0.4, 0.5, 0.6, 0.9, 0.8, 1.0 };
        MonotonicSplineInterpolator interpolator = new MonotonicSplineInterpolator();
        interpolator.addEqualityConstraint(1.0, 0.0);
        interpolator.addEqualityConstraint(9.0, 1.0);
        PolynomialSplineFunction s = interpolator.interpolate(x, y);
        testMonotone(s);
        assertThat(s.value(1.0), closeTo(0.0, 1e-9));
        assertThat(s.value(9.0), closeTo(1.0, 1e-9));
    }

    @Test
    public void testGradientConstraint1() {
        double[] x = { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0 };
        double[] y = { 0.0, 0.2, 0.1, 0.4, 0.5, 0.6, 0.9, 0.8, 1.0 };
        MonotonicSplineInterpolator interpolator = new MonotonicSplineInterpolator();
        interpolator.addGradientEqualityConstraint(5.0, 2.0);
        PolynomialSplineFunction s = interpolator.interpolate(x, y);
        testMonotone(s);
        PolynomialSplineFunction ds = s.polynomialSplineDerivative();
        assertThat(ds.value(5.0), closeTo(2.0, 1e-9));
    }

    @Test
    public void testGradientConstraint2() {
        double[] x = { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0 };
        double[] y = { 0.0, 0.2, 0.1, 0.4, 0.5, 0.6, 0.9, 0.8, 1.0 };
        MonotonicSplineInterpolator interpolator = new MonotonicSplineInterpolator();
        interpolator.addEqualityConstraint(1.0, 0.0);
        interpolator.addEqualityConstraint(9.0, 1.0);
        interpolator.addGradientEqualityConstraint(5.0, 0.321);
        PolynomialSplineFunction s = interpolator.interpolate(x, y);
        testMonotone(s);
        assertThat(s.value(1.0), closeTo(0.0, 1e-9));
        assertThat(s.value(9.0), closeTo(1.0, 1e-9));
        PolynomialSplineFunction ds = s.polynomialSplineDerivative();
        assertThat(ds.value(5.0), closeTo(0.321, 1e-9));
    }

    @Test
    public void testInequalityConstraint1() {
        double[] x = { 1.0, 2.0, 3.0,  4.0, 5.0, 6.0, 7.0,  8.0, 9.0 };
        double[] y = { 0.0, 0.2, 0.05, 0.3, 0.5, 0.7, 0.95, 0.8, 1.0 };
        MonotonicSplineInterpolator interpolator = new MonotonicSplineInterpolator();
        interpolator.addGreaterThanConstraint(1.0, 0.1);
        interpolator.addLessThanConstraint(9.0, 0.9);
        PolynomialSplineFunction s = interpolator.interpolate(x, y);
        testMonotone(s);
        assertThat(s.value(1.0), closeTo(0.1, 1e-9));
        assertThat(s.value(9.0), closeTo(0.9, 1e-9));
    }

    @Test
    public void testInequalityConstraint2() {
        double[] x = { 1.0, 2.0, 3.0,  4.0, 5.0, 6.0, 7.0,  8.0, 9.0 };
        double[] y = { 0.0, 0.2, 0.05, 0.3, 0.5, 0.7, 0.95, 0.8, 1.0 };
        MonotonicSplineInterpolator interpolator = new MonotonicSplineInterpolator();
        interpolator.addGreaterThanConstraint(1.0, 0.1);
        interpolator.addLessThanConstraint(9.0, 0.9);
        interpolator.addEqualityConstraint(5.0, 0.555);
        PolynomialSplineFunction s = interpolator.interpolate(x, y);
        testMonotone(s);
        assertThat(s.value(1.0), closeTo(0.1, 1e-9));
        assertThat(s.value(9.0), closeTo(0.9, 1e-9));
        assertThat(s.value(5.0), closeTo(0.555, 1e-9));
    }

    @Test
    public void testInterpOptions() {
        double[] x = { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0 };
        double[] y = { 0.0, 0.05, 0.02, 0.3, 0.5, 0.7, 0.99, 0.95, 1.0 };
        MonotonicSplineInterpolator interpolator = new MonotonicSplineInterpolator();
        interpolator.addInterpolationOptions(new org.apache.commons.math3.optim.MaxIter(1));
        // Very small max iterations setting should be passed into fitting and cause an exception
        // when the optimization exceeds that number of iterations:
        thrown.expect(org.apache.commons.math3.exception.TooManyIterationsException.class);
        PolynomialSplineFunction s = interpolator.interpolate(x, y);
    }

    @Test
    public void testReproData1() {
        // this setup caused an infinite loop on v0.2.0 of snowball and gibbous
        double[] x = { 2.417586222325241E-6, 0.01833989628165254, 0.03881204328343747, 0.06027210493545283, 0.08641191222046411, 0.08920708170338688, 0.11856432752175541, 0.1271093283793493, 0.15611716035693607, 0.17812340282854217, 0.2003336549968174, 0.21010068304660953, 0.24461903435787338, 0.25330687443961475, 0.2938177183897125, 0.3087980028802898, 0.3234689455608389, 0.3234689455608389, 0.3234689455608389, 0.3874888123642846, 0.4037770547454288, 0.4037770547454288, 0.4451537836690316, 0.45466841689883003, 0.4922791775622433, 0.5049112468735633, 0.5097110552769382, 0.542354984673194, 0.5625051862341103, 0.5812081076929951, 0.5812081076929951, 0.624748141181579, 0.6457634498888773, 0.649465590009825, 0.6792874416836917, 0.7015743319002873, 0.7250186518350656, 0.7257475192866188, 0.7638636260578827, 0.7751265901018766, 0.7751265901018766, 0.8204162734846717, 0.844432325829453, 0.8621117530632927, 0.8621117530632927, 0.8984423523978251, 0.9137609097873816, 0.9453317721241555, 0.9532748155064319, 0.9794713197145162, 0.9999049187132192 };
        double[] y = { 0.0, 0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2, 0.22, 0.24, 0.26, 0.28, 0.3, 0.32, 0.34, 0.36, 0.38, 0.4, 0.42, 0.44, 0.46, 0.48, 0.5, 0.52, 0.54, 0.56, 0.58, 0.6, 0.62, 0.64, 0.66, 0.68, 0.7000000000000001, 0.72, 0.74, 0.76, 0.78, 0.8, 0.8200000000000001, 0.84, 0.86, 0.88, 0.9, 0.92, 0.9400000000000001, 0.96, 0.98, 1.0 };

        double xmin = x[0];
        double xmax = x[x.length - 1];
        double eps = 1e-9;

        MonotonicSplineInterpolator interpolator = new MonotonicSplineInterpolator();
        interpolator.setBounds(xmin, xmax);
        interpolator.setM(20);
        interpolator.addEqualityConstraint(xmin, eps);
        interpolator.addGreaterThanConstraint(xmin, 0.0);
        interpolator.addEqualityConstraint(xmax, 1.0 - eps);
        interpolator.addLessThanConstraint(xmax, 1.0);
        interpolator.addInterpolationOptions(new MaxIter(25));

        PolynomialSplineFunction s = interpolator.interpolate(x, y);
        testMonotone(s);
        assertThat(s.value(xmin), closeTo(0.0, 1e-8));
        assertThat(s.value(xmax), closeTo(1.0, 1e-8));
    }
}
