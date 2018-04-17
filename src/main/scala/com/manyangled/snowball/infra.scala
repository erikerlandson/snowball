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
  // m: the number of piecewise intervals between knots
  // (t0, tm): min and maximum of the interpolation range
  // alpha: (tm - t0) / m, the spacing between knot points
  // (t0, t1, ... tm): the knot points (equally spaced by alpha)
  // n: the number of data points to interpolate (N in the paper)
  // (u1, d1), (u2, d2)... (un, dn): data points to interpolate
  // mm = m + 3 (M in the paper); 3 is from degree-3 (cubic) B-spline
  // lambda: smoothing weight parameter
  // w1, w2, ... wn: interpolation weights for data (larger -> fit closer)

  // returns matrix Rinf (26) for size MxM, where M = m+3
  def matrixRinf(m: Int): Array[Array[Double]] = {
    val mm = m + 3
    // unsure if this is hard minimum but makes logic safe for now, and
    // corresponds to >= 4 knot intervals which doesn't seem overly restrictive
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

  def matrixR(m: Int): Array[Array[Double]] = {
    val mm = m + 3
    // From (26)
    val r = matrixRinf(m)
    // From (27), upper left corner only, delaying 1/6 factor
    val rm = Array(
      Array(14.0,  -6.0,   0.0),
      Array(-6.0,   8.0,  -3.0),
      Array( 0.0,  -3.0,   2.0)
    )
    // From (28), lower right corner only
    val rp = Array(
      Array( 2.0,  -3.0,   0.0),
      Array(-3.0,   8.0,  -6.0),
      Array( 0.0,  -6.0,   14.0)
    )
    // From (25): subtract rm and rp
    for {
      j <- 0 until 3;
      k <- 0 until 3
    } {
      r(j)(k) -= (rm(j)(k) / 6.0)
      r(mm - 3 + j)(mm - 3 + k) -= (rp(j)(k) / 6.0)
    }
    r
  }
}
