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

package com.manyangled.snowball

import com.joptimizer.functions.PDQuadraticMultivariateRealFunction
import com.joptimizer.functions.PSDQuadraticMultivariateRealFunction
import com.joptimizer.functions.ConvexMultivariateRealFunction
import com.joptimizer.functions.LinearMultivariateRealFunction
import com.joptimizer.optimizers.OptimizationRequest
import com.joptimizer.optimizers.JOptimizer

import org.apache.commons.math3.linear.Array2DRowRealMatrix
import org.apache.commons.math3.linear.DiagonalMatrix

/**
 * Implements various components from:
 *
 * H. Fujioka and H. Kano: Monotone smoothing spline curves using normalized uniform cubic B-splines,
 * Trans. Institute of Systems, Control and Information Engineers, Vol. 26, No. 11, pp. 389–397, 2013
 *
 * https://github.com/erikerlandson/snowball/blob/master/monotone-cubic-B-splines-2013.pdf
 * 
 */
object infra {
  // conventions
  // n: the number of data points to interpolate (N in the paper)
  // (u1, d1), (u2, d2)... (un, dn): data points to interpolate
  // w1, w2, ... wn: interpolation weights for data; larger --> fit closer
  // umin, umax: bounds of the interpolation range; e.g. could be (u1,un)
  // m = number of knot intervals covering interpolation interval [umin, umax]
  // mm = the number of knot intervals (M in the paper); mm = m + 3
  // In the paper, knot points ti are indexed from -3, -2, ... m-1
  // knot points t(-3), t(-2), t(-1) lie to the left of umin
  // t(0) = umin.  knots t(0) ... t(m-1) lie on the interval [umin, umax]
  // alpha: 1/alpha is the spacing between knot points; alpha*(t(j+1) - t(j)) = 1; see Eq(2)
  // alpha = m / (umax - umin)
  // lambda: smoothing weight parameter; larger --> smoother, i.e. lower curvature
  // tau: vector of control points.  tau(i) is ith control point.

  case class MonotoneSplineSpec(
    u: Array[Double],
    d: Array[Double],
    m: Int,
    umin: Double,
    umax: Double,
    lambda: Double,
    w: Array[Double],
    tk: Array[Double],
    alpha: Double
  ) {
    def mm = m + 3

    def solve = infra.solve(this)
  }

  def solve(spec: MonotoneSplineSpec) = {
    val (gg, g, r) = qpObjectiveMatrices(spec)
    // Matrix G (aka gg) is basically positive semi-def, and JOptimizer appears to choke if G isn't full rank.
    // Adding 1 along the diagonal restores full rank, and also is equivalently
    // adding (tau.tau) to the cost function, so is effectively requesting "smallest solution" from
    // the underdetermined original matrix.
    for { j <- 0 until gg.length } { gg(j)(j) += 1.0 }
    val objective = new PDQuadraticMultivariateRealFunction(gg, g, r)
    val (hh, t) = qpMonotoneConstraints(spec)
    val ineq: Array[ConvexMultivariateRealFunction] = hh.zip(t).map { case (h, x) =>
      new LinearMultivariateRealFunction(h, x)
    }
    // JOptimizer requires an initial guess that is strictly interior to the constraints
    val ip = interiorPoint(ineq)
    val oreq = new OptimizationRequest()
    oreq.setInitialPoint(ip)
    oreq.setF0(objective)
    oreq.setFi(ineq)
    oreq.setToleranceFeas(1e-10)
    oreq.setTolerance(1e-10)
    val opt = new JOptimizer
    opt.setOptimizationRequest(oreq)
    opt.optimize()
    opt.getOptimizationResponse().getSolution()
  }

  // This is what I've come to; generate and test.
  def interiorPoint(ineq: Array[ConvexMultivariateRealFunction], maxTries: Int = 1000000): Array[Double] = {
    var found = false
    var tries = 0
    var ip: Array[Double] = null
    val m = ineq(0).asInstanceOf[LinearMultivariateRealFunction].getQ().size().toInt
    while (!found) {
      tries += 1
      if (tries > maxTries) {
        throw new Exception(s"interiorPoint failed after $maxTries attempts")
      }
      val t = Array.fill(m) { scala.util.Random.nextDouble() }
      val v = new cern.colt.matrix.tdouble.impl.DenseDoubleMatrix1D(t)
      val sat = ineq.forall { q => q.value(v) < 0.0 }
      if (sat) {
        println(s"found interior point after $tries attempts")
        ip = t
        found = true
      }
    }
    ip
  }

  // returns matrix Rinf Eq(26) for size MxM, where M = mm
  def matrixRinf(mm: Int): Array[Array[Double]] = {
    // unsure if this is hard minimum but makes logic safe for now
    require(mm >= 7)
    val d0 = 16.0 / 6.0
    val d1 = -9.0 / 6.0
    val d2 = 0.0 / 6.0
    val d3 = 1.0 / 6.0
    val band = Array(d3, d2, d1, d0, d1, d2, d3)
    val b = band.length
    val rf = Array.tabulate(3) { jj =>
      val j = 3 - jj
      band.drop(j) ++ Array.fill(mm - (b - j))(0.0)
    }
    val rl = Array.tabulate(3) { jj =>
      val j = jj + 1
      Array.fill(mm - (b - j))(0.0) ++ band.dropRight(j)
    }
    val rm = Array.tabulate(mm - (2 * 3)) { j =>
      Array.fill(j)(0.0) ++ band ++ Array.fill(mm - (b + j))(0.0)
    }
    rf ++ rm ++ rl
  }

  def matrixR(mm: Int): Array[Array[Double]] = {
    // From Eq(26)
    val r = matrixRinf(mm)
    // From Eq(27), upper left corner only, delaying 1/6 factor
    val rm = Array(
      Array(14.0,  -6.0,   0.0),
      Array(-6.0,   8.0,  -3.0),
      Array( 0.0,  -3.0,   2.0)
    )
    // From Eq(28), lower right corner only
    val rp = Array(
      Array( 2.0,  -3.0,   0.0),
      Array(-3.0,   8.0,  -6.0),
      Array( 0.0,  -6.0,   14.0)
    )
    // From Eq(25): subtract rm and rp
    for {
      j <- 0 until 3;
      k <- 0 until 3
    } {
      r(j)(k) -= (rm(j)(k) / 6.0)
      r(mm - 3 + j)(mm - 3 + k) -= (rp(j)(k) / 6.0)
    }
    r
  }

  // basis B3(t) from Eq(3) and Table 1
  // Here, t is treated as normalized to unit distance between knot points
  def b3(t: Double) = t match {
    case t if t < 0.0 => 0.0
    case t if t < 1 => math.pow(1.0 - t, 3) / 6.0
    case t if t < 2 => (4.0 - (6.0 * t * t) + (3.0 * math.pow(t, 3))) / 6.0
    case t if t < 3 => (1.0 + (3.0 * t) + (3.0 * t * t) - (3.0 * math.pow(t, 3))) / 6.0
    case t if t < 4 => math.pow(t, 3) / 6.0
    case _ => 0.0
  }

  // Eq(11)
  // tk is array of knot points, dimension mm, aka M
  def vectorb(alpha: Double, tk: Array[Double])(t: Double): Array[Double] =
    tk.map { ti => b3(alpha * (t - ti)) }

  // tk is array of knot points, dimension mm, aka M
  // note this returns B-transposed relative to Eq(12).
  // u is vector of data "x" points; in same axis as knot points tk
  def matrixBT(alpha: Double, tk: Array[Double])(u: Array[Double]): Array[Array[Double]] = {
    val b = vectorb(alpha, tk)(_)
    u.map(b)
  }

  // constructs G, g and r from Eq(14) and Eq(68)
  def qpObjectiveMatrices(spec: MonotoneSplineSpec) = {
    val r = 0.0  // per Eq(68) - constant does not change location of minimum
    val BT = new Array2DRowRealMatrix(matrixBT(spec.alpha, spec.tk)(spec.u), false)
    val B = BT.transpose()
    val R = new Array2DRowRealMatrix(matrixR(spec.mm), false)
    val lambdaQ = R.scalarMultiply(spec.lambda * math.pow(spec.alpha, 3)) // Eq(15) & Eq(22)
    val W = new DiagonalMatrix(spec.w, false)
    val G = lambdaQ.add(B.multiply(W).multiply(BT)) // Eq(15)
    val g = B.multiply(W).operate(spec.d) // Eq(16); 'operate' method is "multiply by vector"
    (G.getData(), g, r)
  }

  def qpMonotoneConstraints(spec: MonotoneSplineSpec) = {
    // negative of Eq(64), because JOptimizer expects constraints of form Hx <= 0 instead of Hx >= 0
    val Fm3 = Array(1.0, 0.0, -1.0,  0.0)
    val Fm2 = Array(1.0, 2.0, -3.0,  0.0)
    val Fm1 = Array(0.0, 1.0,  0.0, -1.0)
    val Fm0 = Array(0.0, 3.0, -2.0, -1.0)

    // apply constraints over the interpolation interval (starting at knot interval [t0, t1])
    // Eq(66), replicated for each interval starting at t0, t1, ... t(m-1)
    val H = scala.collection.mutable.ArrayBuffer.empty[Array[Double]]
    for { j <- 0 until spec.m } {
      H += Array.fill(j)(0.0) ++ Fm3 ++ Array.fill(spec.mm - (j + 4))(0.0)
      H += Array.fill(j)(0.0) ++ Fm2 ++ Array.fill(spec.mm - (j + 4))(0.0)
      H += Array.fill(j)(0.0) ++ Fm1 ++ Array.fill(spec.mm - (j + 4))(0.0)
      H += Array.fill(j)(0.0) ++ Fm0 ++ Array.fill(spec.mm - (j + 4))(0.0)
    }

    (H.toArray, Array.fill(H.length)(0.0))
  }
}
