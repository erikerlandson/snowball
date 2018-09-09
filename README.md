# snowball
Monotonic smoothing splines for the JVM ecosystem and Apache Commons Math.

### Documentation

API javadoc is available at:
[https://erikerlandson.github.io/snowball/java/api/](https://erikerlandson.github.io/snowball/java/api/)

### How to use `snowball` in your project
The `snowball` package is implemented in java, and so it can be used in both java and scala. It is built on, and designed to work with, Apache Commons Math 3.6

`snowball` expects you to provide `commons-math3` and `gibbous` dependencies, as shown here:

```scala
resolvers += "manyangled" at "https://dl.bintray.com/manyangled/maven/"

libraryDependencies ++= Seq(
  "com.manyangled" % "snowball" % "0.1.0",
  "com.manyangled" % "gibbous" % "0.1.0",
  "org.apache.commons" % "commons-math3" % "3.6.1")
```

### References:
1. H. Fujioka and H. Kano: [Monotone smoothing spline curves using normalized uniform cubic B-splines](/monotone-cubic-B-splines-2013.pdf), Trans. Institute of Systems, Control and Information Engineers, Vol. 26, No. 11, pp. 389–397, 2013

1. Hiroyuki KANO, Hiroyuki FUJIOKA, and Clyde F. MARTIN, [Optimal Smoothing Spline with Constraints on Its Derivatives](https://www.jstage.jst.go.jp/article/jcmsi/7/2/7_104/_pdf), SICE Journal of Control, Measurement, and System Integration, Vol.7, No. 2, pp. 104–111, March 2014

1. M. Nagahara, Y. Yamamoto, C. Martin, [Quadratic Programming for Monotone Control Theoretic Splines](https://www.researchgate.net/profile/Clyde_Martin/publication/224182849_Quadratic_programming_for_monotone_control_theoretic_splines/links/00b7d52da8b1e52d6c000000/Quadratic-programming-for-monotone-control-theoretic-splines.pdf), SICE, 2010.

1. M. Egerstedt and C. Martin. [Monotone Smoothing Splines](http://magnus.ece.gatech.edu/Papers/MonoSplines.pdf). Mathematical Theory of Networks and Systems. Perpignan, France, June 2000.
