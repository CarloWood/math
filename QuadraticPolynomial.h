#pragma once

#include "ConceptMathOperations.h"
#include "utils/square.h"
#include "utils/macros.h"
#include <array>
#include "debug.h"
#if CW_DEBUG
#include "utils/almost_equal.h"
#endif
#ifdef CWDEBUG
#include "utils/has_print_on.h"
#endif

namespace math {
#ifdef CWDEBUG
using utils::has_print_on::operator<<;
#endif

template<typename T>
requires ConceptMathOperations<T>
class QuadraticPolynomial
{
 private:
  std::array<T, 3> coefficients_{};

 public:
  // Create a zero polynomial.
  QuadraticPolynomial() { }
  // Create a polynomial  a + b x + c x^2.
  QuadraticPolynomial(T a, T b, T c) : coefficients_{{a, b, c}} { }

  // Evaluation.
  T operator()(T w) const
  {
    return coefficients_[0] + (coefficients_[1] + coefficients_[2] * w) * w;
  };

  // Evaluate derivative.
  T derivative(T w) const
  {
    return coefficients_[1] + 2 * coefficients_[2] * w;
  }

  // Access coefficients.
  T operator[](int i) const { ASSERT(0 <= i && i < coefficients_.size()); return coefficients_[i]; }
  T& operator[](int i) { ASSERT(0 <= i && i < coefficients_.size()); return coefficients_[i]; }

  int get_roots(std::array<T, 2>& roots_out) const
  {
    if (coefficients_[2] == 0)
    {
      roots_out[0] = -coefficients_[0] / coefficients_[1];
      return isfinite(roots_out[0]) ? 1 : 0;
    }

    T const D = utils::square(coefficients_[1]) - 4 * coefficients_[2] * coefficients_[0];
    if (D < 0)
      return 0;

    // If this fails then probably the coefficients are SO small that both utils::square(coefficients_[1])
    // as well as coefficients_[2] * coefficients_[0] are zero. That is "fine" as long as the discriminant
    // is not, in fact, less than zero.
    ASSERT(D > 0 || abs(coefficients_[1]) >= 2 * sqrt(abs(coefficients_[2])) * sqrt(abs(coefficients_[0])) ||
        utils::almost_equal(abs(coefficients_[1]),
          2 * sqrt(abs(coefficients_[2])) * sqrt(abs(coefficients_[0])), 1e-14));

    // Use a sqrt with the same sign as coefficients_[1];
    T const signed_sqrt_D = copysign(sqrt(D), coefficients_[1]);

    // Calculate the root closest to zero.
    roots_out[0] = -2 * coefficients_[0] / (coefficients_[1] + signed_sqrt_D);

    if (AI_UNLIKELY(isnan(roots_out[0])))
    {
      // This means we must have divided by zero, which means that both, coefficients_[1] as well as sqrtD, must be zero.
      // The latter means that coefficients_[0] is zero (coefficients_[2] was already checked not to be zero).
      // Therefore we have: f(x) = c x^2 with one root at x=0.
      roots_out[0] = 0;
      return 1;
    }

    // Calculate the root further away from zero.
    roots_out[1] = (coefficients_[1] + signed_sqrt_D) / (-2 * coefficients_[2]);

    // The second one is larger in absolute value.
    ASSERT(abs(roots_out[1]) > abs(roots_out[0]) || utils::almost_equal(abs(roots_out[0]), abs(roots_out[1]), 1e-14));

    return 2;
  }

  // Returns the x coordinate of the vertex of the parabola.
  T vertex_x() const
  {
    // f(x) = a + bx + cx^2
    // f'(x) = b + 2c x
    // f'(x) = 0 --> x = -b / 2c
    return coefficients_[1] / (-2 * coefficients_[2]);
  }

  // Returns the y coordinate of the vertex of the parabola.
  T vertex_y() const
  {
    // f(vertex_x()) = a + b(-b / 2c) + c(-b / 2c)^2 = a - b^2 / 2c + b^2 / 4c = a - b^2 / 4c.
    return coefficients_[0] - utils::square(coefficients_[1]) / (4 * coefficients_[2]);
  }

  QuadraticPolynomial height() const
  {
    // f(x) - vertex_y()
    QuadraticPolynomial result(*this);
    result[0] = utils::square(coefficients_[1]) / (4 * coefficients_[2]);
    return result;
  }

  QuadraticPolynomial& operator-=(QuadraticPolynomial const& rhs)
  {
    for (int i = 0; i < coefficients_.size(); ++i)
      coefficients_[i] -= rhs.coefficients_[i];
    return *this;
  }

  QuadraticPolynomial& operator+=(QuadraticPolynomial const& rhs)
  {
    for (int i = 0; i < coefficients_.size(); ++i)
      coefficients_[i] += rhs.coefficients_[i];
    return *this;
  }

  QuadraticPolynomial operator+(QuadraticPolynomial const& rhs) const
  {
    QuadraticPolynomial result(*this);
    result += rhs;
    return result;
  }

  QuadraticPolynomial operator-(QuadraticPolynomial const& rhs) const
  {
    QuadraticPolynomial result(*this);
    result -= rhs;
    return result;
  }

  QuadraticPolynomial& operator*=(T factor)
  {
    for (int i = 0; i < coefficients_.size(); ++i)
      coefficients_[i] *= factor;
    return *this;
  }

  friend QuadraticPolynomial operator*(T factor, QuadraticPolynomial const& rhs)
  {
    QuadraticPolynomial result(rhs);
    result *= factor;
    return result;
  }

  // Returns true if the old parabola is close enough to this parabola at infinity.
  // toggles_out is filled with an even number of x-coordinates at which the
  // y-coordinate of the old parabola goes in or out the acceptable range.
  bool equal_intervals(QuadraticPolynomial const& old, std::array<T, 4>& toggles_out, int& count_out) const
  {
    // This object is the most recent parabola, the one that came after old.

    // Create a new parabola, diff = 10 * (new - old).
    QuadraticPolynomial diff{10 * (coefficients_[0] - old[0]), 10 * (coefficients_[1] - old[1]), 10 * (coefficients_[2] - old[2])};

    // If abs(diff[2]) < abs(new[2]) then the old one is close enough at infinity.
    bool close_at_inf = abs(diff[2]) < abs(coefficients_[2]);

    // Set v_y to the y-coordinate of the vertex of the new parabola.
    T v_y = vertex_y();

    std::array<std::array<T, 2>, 2> toggle;
    std::array<int, 2> count{};
    std::array<T, 2> c{coefficients_[2] - diff.coefficients_[2], coefficients_[2] + diff.coefficients_[2]};
    std::array<T, 2> abs_c{abs(c[0]), abs(c[1])};
    for (int i = 0; i < 2; ++i) // i=0: subtract diff, i=1: add diff.
    {
      T sign = i == 0 ? -1 : 1;
      if (abs_c[i] > 1e-30)
      {
        T b = coefficients_[1] + sign * diff.coefficients_[1];
        T D = utils::square(b) - 4 * (coefficients_[0] - v_y + sign * diff.coefficients_[0]) * c[i];
        if (D > 0)
        {
          T avg = b / (-2 * c[i]);
          T delta = sqrt(D) / (2 * abs_c[i]);
          toggle[i][0] = avg - delta;
          toggle[i][1] = avg + delta;
          count[i] = 2;
        }
      }
    }
    count_out = count[0] + count[1];
    int i0 = 0;
    int i1 = 0;
    int i = 0;

    // Sort the result and write it to toggles_out.
    while (i < count_out)
    {
      if (i1 == count[1] || (i0 < count[0] && toggle[0][i0] < toggle[1][i1]))
        toggles_out[i] = toggle[0][i0++];
      else
        toggles_out[i] = toggle[1][i1++];
      ++i;
    }

    ASSERT(count_out == 0 ||
           (count_out == 2 && toggles_out[0] < toggles_out[1]) ||
           (count_out == 4 && toggles_out[0] < toggles_out[1] && toggles_out[1] < toggles_out[2] && toggles_out[2] < toggles_out[3]));

    return close_at_inf;
  }

#ifdef CWDEBUG
  void print_on(std::ostream& os) const
  {
    bool first = true;
    int exponent = 0;
    for (T coefficient : coefficients_)
    {
      if (coefficient != 0)
      {
        if (first)
          os << coefficient;
        else if (coefficient > 0)
          os << " + " << coefficient;
        else
          os << " - " << -coefficient;
        if (exponent > 0)
        {
          os << " x";
          if (exponent > 1)
            os << '^' << exponent;
        }
        first = false;
      }
      ++exponent;
    }
  }
#endif
};

} // namespace math
