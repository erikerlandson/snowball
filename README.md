# snowball
Monotonic smoothing splines for the JVM ecosystem and Apache Commons Math.

### Documentation

API javadoc is available at:
[https://erikerlandson.github.io/snowball/java/api/](https://erikerlandson.github.io/snowball/java/api/)

A few examples are below.

### Features

* Fit monotonic interpolating splines to data, including data that has noise or is otherwise non-monotonic.
* Enforce equality constraints of the form s(x) = y, where s is the spline function
* Enforce gradient constraints of the form ds(x)/dx = g
* Enforce inequality constraints of the form s(x) < y and s(x) > y

### How to use `snowball` in your project
The `snowball` package is implemented in java, and so it can be used in both java and scala. It is built on, and designed to work with, Apache Commons Math 3.6.

#### using SBT
```scala
libraryDependencies ++= Seq(
  "com.manyangled" % "snowball" % "0.3.0"
  )
```

#### using maven
```xml
<dependency>
  <groupId>com.manyangled</groupId>
  <artifactId>snowball</artifactId>
  <version>0.3.0</version>
  <type>pom</type>
</dependency>
<dependency>
  <groupId>com.manyangled</groupId>
  <artifactId>gibbous</artifactId>
  <version>0.3.0</version>
  <type>pom</type>
</dependency>
```

### Examples

#### Java
```java
import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction;
import com.manyangled.snowball.analysis.interpolation.MonotonicSplineInterpolator;

double[] x = { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0 };
double[] y = { 0.0, 0.05, 0.02, 0.3, 0.5, 0.7, 0.99, 0.95, 1.0 };
MonotonicSplineInterpolator interpolator = new MonotonicSplineInterpolator();
PolynomialSplineFunction s = interpolator.interpolate(x, y);
```

#### Scala REPL
```scala
scala> import com.manyangled.snowball.analysis.interpolation._, com.manyangled.gnuplot4s._
import com.manyangled.snowball.analysis.interpolation._
import com.manyangled.gnuplot4s._

scala> val interpolator = new MonotonicSplineInterpolator()
interpolator: com.manyangled.snowball.analysis.interpolation.MonotonicSplineInterpolator = com.manyangled.snowball.analysis.interpolation.MonotonicSplineInterpolator@6834fd1b

scala> val xdata = Array(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0)
xdata: Array[Double] = Array(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0)

scala> val ydata = Array(0.0, 0.2, 0.05, 0.3, 0.5, 0.7, 0.95, 0.8, 1.0)
ydata: Array[Double] = Array(0.0, 0.2, 0.05, 0.3, 0.5, 0.7, 0.95, 0.8, 1.0)

scala> val s = interpolator.interpolate(xdata, ydata)
s: org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction = org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction@5852d898

scala> Session().block("data", xdata.zip(ydata)).block("spline", (1.0 to 9.0 by 0.1).map { x => (x, s.value(x)) }).plot(Plot().block("data").using(1,2).style(PlotStyle.Points)).plot(Plot().block("spline").using(1,2).style(PlotStyle.Lines)).term(Dumb().size(80,40)).render

scala> 
                                                                                       
                                                                                
    1 +-+------+--------+-------+--------+--------+--------+-------+------+-A   
      +        +        +       +        +        +        +       +     ####   
      |                                             $data uAing 1:2   #A#   |   
      |                                           $spline using 1:######### |   
      |                                                        ###          |   
      |                                                     ###             |   
      |                                                   ##                |   
  0.8 +-+                                                ##        A      +-+   
      |                                               ###                   |   
      |                                             ##                      |   
      |                                           A#                        |   
      |                                           #                         |   
      |                                         ##                          |   
      |                                       ##                            |   
  0.6 +-+                                    #                            +-+   
      |                                     #                               |   
      |                                    #                                |   
      |                                  A#                                 |   
      |                                 ##                                  |   
      |                                #                                    |   
      |                               #                                     |   
  0.4 +-+                            #                                    +-+   
      |                            ##                                       |   
      |                          ##                                         |   
      |                         A                                           |   
      |                        #                                            |   
      |                      ##                                             |   
      |                   ###                                               |   
  0.2 +-+      A        ##                                                +-+   
      |                ##                                                   |   
      |             ###                                                     |   
      |          ###                                                        |   
      |      ####                                                           |   
      |   ###           A                                                   |   
      ####     +        +       +        +        +        +       +        +   
    0 A-+------+--------+-------+--------+--------+--------+-------+------+-+   
      1        2        3       4        5        6        7       8        9   
                                                                                
```

### References:
1. H. Fujioka and H. Kano: [Monotone smoothing spline curves using normalized uniform cubic B-splines](/monotone-cubic-B-splines-2013.pdf), Trans. Institute of Systems, Control and Information Engineers, Vol. 26, No. 11, pp. 389–397, 2013

1. Hiroyuki KANO, Hiroyuki FUJIOKA, and Clyde F. MARTIN, [Optimal Smoothing Spline with Constraints on Its Derivatives](https://www.jstage.jst.go.jp/article/jcmsi/7/2/7_104/_pdf), SICE Journal of Control, Measurement, and System Integration, Vol.7, No. 2, pp. 104–111, March 2014

1. M. Nagahara, Y. Yamamoto, C. Martin, [Quadratic Programming for Monotone Control Theoretic Splines](https://www.researchgate.net/profile/Clyde_Martin/publication/224182849_Quadratic_programming_for_monotone_control_theoretic_splines/links/00b7d52da8b1e52d6c000000/Quadratic-programming-for-monotone-control-theoretic-splines.pdf), SICE, 2010.

1. M. Egerstedt and C. Martin. [Monotone Smoothing Splines](http://magnus.ece.gatech.edu/Papers/MonoSplines.pdf). Mathematical Theory of Networks and Systems. Perpignan, France, June 2000.
