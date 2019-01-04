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
}
