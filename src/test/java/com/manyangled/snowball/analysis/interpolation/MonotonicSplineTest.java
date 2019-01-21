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

import java.util.Arrays;

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

import org.apache.commons.math3.distribution.RealDistribution;
import org.apache.commons.math3.distribution.UniformRealDistribution;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.distribution.GammaDistribution;

public class MonotonicSplineTest {
    static final double epsEq = 1e-8;

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
        assertThat(s.value(5.0), closeTo(0.7, epsEq));
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
        assertThat(s.value(1.0), closeTo(0.0, epsEq));
        assertThat(s.value(9.0), closeTo(1.0, epsEq));
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
        assertThat(ds.value(5.0), closeTo(2.0, epsEq));
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
        assertThat(s.value(1.0), closeTo(0.0, epsEq));
        assertThat(s.value(9.0), closeTo(1.0, epsEq));
        PolynomialSplineFunction ds = s.polynomialSplineDerivative();
        assertThat(ds.value(5.0), closeTo(0.321, epsEq));
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
        assertThat(s.value(1.0), closeTo(0.1, epsEq));
        assertThat(s.value(9.0), closeTo(0.9, epsEq));
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
        assertThat(s.value(1.0), closeTo(0.1, epsEq));
        assertThat(s.value(9.0), closeTo(0.9, epsEq));
        assertThat(s.value(5.0), closeTo(0.555, epsEq));
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
        assertThat(s.value(xmin), closeTo(0.0, epsEq));
        assertThat(s.value(xmax), closeTo(1.0, epsEq));
    }

    @Test
    public void testReproData2() {
        // this setup caused an infinite loop on v0.2.1 of snowball and gibbous
        double[] x = { -4.184507076820324, -2.1146902974778965, -1.7709154201703174, -1.5333773805478366, -1.3456728173480172, -1.2437838008041608, -1.233006036263326, -1.073638975747641, -0.9803883478450147, -0.9033586284461634, -0.8916327289429522, -0.7504705542249223, -0.695307246427601, -0.6411508285039249, -0.6276275581243719, -0.5127434127323376, -0.44710629582025563, -0.3853137917199486, -0.3853137917199486, -0.3485233683113338, -0.3008475838732653, -0.22625442301493448, -0.15166126215660358, -0.14303406197096172, -0.14303406197096172, 0.0026039423924754945, 0.053096546820946555, 0.0993039111521562, 0.0993039111521562, 0.0993039111521562, 0.2614257843723851, 0.2976054000388951, 0.3696609138599017, 0.41681589265222346, 0.4810240984571372, 0.5372644961355958, 0.5791825316605546, 0.5791825316605546, 0.7424373299418614, 0.7530096700595077, 0.7530096700595077, 0.906789394649255, 0.9348353467475007, 1.0732676924300075, 1.1096789897786465, 1.2573676114538015, 1.396016513175573, 1.5142180513519905, 1.7248140344744871, 2.0772952939655918, 3.7276747214206036 };
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
        assertThat(s.value(xmin), closeTo(0.0, epsEq));
        assertThat(s.value(xmax), closeTo(1.0, epsEq));
    }

    @Test
    public void testReproData3() {
        double[] x = { 1.00005000333357E-4, 0.019409251359818067, 0.04068123914173336, 0.06123553535520177, 0.08444390461077006, 0.10548305587813296, 0.12805036553141316, 0.15022912395741614, 0.17343564129838546, 0.19945618172583213, 0.22415461197566466, 0.24870841855300307, 0.27466430083894233, 0.30210172980035993, 0.32823392882926855, 0.3554836985273914, 0.3858289402923884, 0.41632088667531236, 0.44693110468213626, 0.4791849849063438, 0.5122592850560248, 0.5457298342288471, 0.5798774432413372, 0.6180240086177533, 0.6537701765373716, 0.6915901126921077, 0.7320916945823159, 0.775624812882162, 0.819289912175568, 0.865797828169946, 0.9138084914304859, 0.9700061938204145, 1.0218370263297414, 1.0794564486050406, 1.1378577173142972, 1.2065319159171217, 1.2747986658568298, 1.3460486046891293, 1.4288468320054415, 1.512114305837601, 1.61132514870369, 1.7160871295537956, 1.8279229962456165, 1.9702377298309495, 2.125429826250494, 2.306241157048467, 2.5295702641685995, 2.807343679326415, 3.2358050235131897, 3.953996926202895, 7.878844401611691 };
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
        assertThat(s.value(xmin), closeTo(0.0, epsEq));
        assertThat(s.value(xmax), closeTo(1.0, epsEq));
    }

    void splineTrials(final RealDistribution dist, final int trials, final int intervals) {
        final double eps = 1e-9;
        double[] x = new double[1 + intervals];
        double[] y = new double[1 + intervals];
        final double step = 1.0 / (double)intervals;
        final RealDistribution u = new UniformRealDistribution(-step / 20.0, step / 20.0);
        for (int j = 0; j < intervals; ++j) y[j] = step * (double)j;
        y[intervals] = 1.0;
        for (int ii = 0; ii < trials; ++ii) {
            for (int j = 0; j <= intervals; ++j) {
                double yy = y[j] + u.sample();
                yy = Math.max(0.0001, yy);
                yy = Math.min(0.9999, yy);
                x[j] = dist.inverseCumulativeProbability(yy);
            }
            final double xmin = x[0];
            final double xmax = x[intervals];
            MonotonicSplineInterpolator interpolator = new MonotonicSplineInterpolator();
            // set the bounds of interpolation to the data range
            interpolator.setBounds(xmin, xmax);
            // set 20 splining intervals
            interpolator.setM(20);
            // these constraints fix cdf(xmin) to be "effectively zero" and also enforce > 0
            interpolator.addEqualityConstraint(xmin, eps);
            interpolator.addGreaterThanConstraint(xmin, 0.0);
            // these constraints fix cdf(xmax) to be "effectively one" and also enforce < 1
            interpolator.addEqualityConstraint(xmax, 1.0 - eps);
            interpolator.addLessThanConstraint(xmax, 1.0);
            interpolator.addInterpolationOptions(new MaxIter(500));
            // get the splines that approximate CDF and the PDF
            PolynomialSplineFunction s = null;
            try {
                s = interpolator.interpolate(x, y);
            } catch (Exception e) {
                System.err.format("Splining trial FAILED:\n");
                System.err.format("x= %s\n", Arrays.toString(x));
                System.err.format("y= %s\n", Arrays.toString(y));
                throw e;
            }
            testMonotone(s);
            assertThat(s.value(xmin + eps), closeTo(0.0, epsEq));
            assertThat(s.value(xmax - eps), closeTo(1.0, epsEq));
        }
    }

    @Test
    public void testSplineTrialsUniform() {
        splineTrials(new UniformRealDistribution(), 100, 50);
    }

    @Test
    public void testSplineTrialsNormal() {
        splineTrials(new NormalDistribution(), 100, 50);
    }

    @Test
    public void testSplineTrialsGamma() {
        splineTrials(new GammaDistribution(1.0, 1.0), 100, 50);
    }
}
