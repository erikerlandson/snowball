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

package com.manyangled

package object snowball {

  implicit class AAShow(m: Array[Array[Double]]) {
    def show = m.map(_.mkString(", ")).mkString("\n")
  }
  implicit class AShow(v: Array[Double]) {
    def show = v.mkString(", ")
  }

  import infra.MonotoneSplineSpec

  object MonotoneSplineSetup {
    def apply(
      data: Seq[(Double, Double)],
      m: Int,
      bounds: Option[(Double, Double)] = None,
      lambda: Option[Double] = None,
      w: Option[Seq[Double]] = None
    ): MonotoneSplineSpec = {
      require(m >= 4)
      require(data.length >= m)
      val u = data.map(_._1).toArray
      val d = data.map(_._2).toArray
      val ll = lambda.getOrElse(1.0)
      require(ll > 0.0)
      val (umin, umax) = bounds.getOrElse((u.min, u.max))
      require(umin < umax)
      val ww = w.map(_.toArray).getOrElse(Array.fill(u.length)(1.0))
      require(ww.length == u.length)
      require(ww.forall(_ > 0.0))
      val alpha = m.toDouble / (umax - umin)
      val tk = (-3 to (m - 1)).toArray.map { j => umin + (j.toDouble / alpha) }
      MonotoneSplineSpec(u, d, m, umin, umax, ll, ww, tk, alpha)
    }
  }

  object test {
    val data = Vector(
      (0.0, 0.0), (1.0, 0.1), (2.0, 0.4), (3.0, 0.5), (4.0, 0.6), (5.0, 0.9), (6.0, 1.0)
    )
  }

  object testJO {
    import com.joptimizer.functions.PDQuadraticMultivariateRealFunction
    import com.joptimizer.functions.PSDQuadraticMultivariateRealFunction
    import com.joptimizer.functions.ConvexMultivariateRealFunction
    import com.joptimizer.functions.LinearMultivariateRealFunction
    import com.joptimizer.optimizers.OptimizationRequest
    import com.joptimizer.optimizers.JOptimizer

    // solution space is dimension n; in this example n = 2
    // matrix P is n X n
    val P = Array(
      Array(1.0, 0.4),
      Array(0.4, 1.0)
    )
    // vector q is length n
    val q = Array(0.0, 0.0)
    // r is a scalar
    val r = 0.0
    // My P is positive definite but I'm testing if it will take positive semi-def
    val objective = new PSDQuadraticMultivariateRealFunction(P, q, r)

    // A is k X n, where k is number of equalities. Here k = 1
    // vector b has dimension k
    // conceptually, setting up a matrix equation Ax = b
    // This constraint corresponds to x[0] + x[1] = 1
    val A = Array(Array(1.0, 1.0))
    val b = Array(1.0)

    // Convex inequalities. These are typical linear inequalities.
    // Each linear inequality corresponds to form: c.x <= t , where c is
    // a vector of length n and t is a scalar threshold
    // The following corresponds to x[0] >= 0.6, x[1] >= 0
    val ineq: Array[ConvexMultivariateRealFunction] = Array(
      new LinearMultivariateRealFunction(Array(-1.0, 0.0), 0.6),
      new LinearMultivariateRealFunction(Array(0.0, -1.0), 0.0)
    )

    val oreq = new OptimizationRequest()
    oreq.setF0(objective)
    // set this to a strictly feasible point
    // here, we have an equality constraint x[0] + x[1] = 1, pick a point satisfying that
    // Also we have x[0] >= 0.6, so make sure x[0] in this initial point is >= 0.6
    oreq.setInitialPoint(Array(0.9, 0.1))
    oreq.setFi(ineq)
    oreq.setA(A)
    oreq.setB(b)
    // these seem pretty aggressive, consider loosening a bit in real life?
    oreq.setToleranceFeas(1e-12)
    oreq.setTolerance(1e-12)

    // given the problem definition above, this should return x = [0.6, 0.4]
    def run = {
      val opt = new JOptimizer()
      opt.setOptimizationRequest(oreq)
      opt.optimize()
      opt.getOptimizationResponse().getSolution()
    }
  }
}
