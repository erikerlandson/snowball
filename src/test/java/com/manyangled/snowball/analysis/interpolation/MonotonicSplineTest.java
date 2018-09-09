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

    @Test
    public void test1() {
        double[] x = { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0 };
        double[] y = { 0.0, 0.3, 0.6, 0.5, 0.7, 1.0 };
        MonotonicSplineInterpolator interpolator = new MonotonicSplineInterpolator();
        PolynomialSplineFunction s = interpolator.interpolate(x, y);
        testMonotone(s);
    }

    @Test
    public void testEqualityConstraint1() {
        double[] x = { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0 };
        double[] y = { 0.0, 0.3, 0.6, 0.5, 0.7, 1.0 };
        MonotonicSplineInterpolator interpolator = new MonotonicSplineInterpolator();
        interpolator.addEqualityConstraint(3.5, 0.7);
        PolynomialSplineFunction s = interpolator.interpolate(x, y);
        testMonotone(s);
        assertThat(s.value(3.5), closeTo(0.7, 1e-9));
    }

    @Test
    public void testEqualityConstraint2() {
        double[] x = { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0 };
        double[] y = { 0.0, 0.3, 0.6, 0.5, 0.7, 1.0 };
        MonotonicSplineInterpolator interpolator = new MonotonicSplineInterpolator();
        interpolator.addEqualityConstraint(1.0, 0.0);
        interpolator.addEqualityConstraint(6.0, 1.0);
        PolynomialSplineFunction s = interpolator.interpolate(x, y);
        testMonotone(s);
        assertThat(s.value(1.0), closeTo(0.0, 1e-9));
        assertThat(s.value(6.0), closeTo(1.0, 1e-9));
    }

    @Test
    public void testGradientConstraint1() {
        double[] x = { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0 };
        double[] y = { 0.0, 0.3, 0.6, 0.5, 0.7, 1.0 };
        MonotonicSplineInterpolator interpolator = new MonotonicSplineInterpolator();
        interpolator.addGradientEqualityConstraint(3.5, 2.0);
        PolynomialSplineFunction s = interpolator.interpolate(x, y);
        testMonotone(s);
        PolynomialSplineFunction ds = s.polynomialSplineDerivative();
        assertThat(ds.value(3.5), closeTo(2.0, 1e-9));
    }

    @Test
    public void testGradientConstraint2() {
        double[] x = { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0 };
        double[] y = { 0.0, 0.3, 0.6, 0.5, 0.7, 1.0 };
        MonotonicSplineInterpolator interpolator = new MonotonicSplineInterpolator();
        interpolator.addEqualityConstraint(1.0, 0.0);
        interpolator.addEqualityConstraint(6.0, 1.0);
        interpolator.addGradientEqualityConstraint(3.5, 0.4567);
        PolynomialSplineFunction s = interpolator.interpolate(x, y);
        testMonotone(s);
        assertThat(s.value(1.0), closeTo(0.0, 1e-9));
        assertThat(s.value(6.0), closeTo(1.0, 1e-9));
        PolynomialSplineFunction ds = s.polynomialSplineDerivative();
        assertThat(ds.value(3.5), closeTo(0.4567, 1e-9));
    }
}
