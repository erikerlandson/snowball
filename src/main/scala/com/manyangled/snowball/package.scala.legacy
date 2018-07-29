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
}
