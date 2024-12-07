#pragma once

#include "Polynomial.h"
#include "QuadraticPolynomial.h"
#include "utils/square.h"
#include "utils/macros.h"
#include <array>
#include <cmath>
#include <cstring>
#include "debug.h"
#ifdef CWDEBUG
#include "utils/has_print_on.h"

NAMESPACE_DEBUG_CHANNELS_START
extern channel_ct cubic;
NAMESPACE_DEBUG_CHANNELS_END
#endif // CWDEBUG

namespace math {
#ifdef CWDEBUG
using utils::has_print_on::operator<<;
#endif

template<typename T>
requires ConceptMathOperations<T>
class CubicPolynomial
{
 protected:
  std::array<T, 4> coefficients_{};

 public:
  // Create a zero polynomial.
  CubicPolynomial() { }
  // Create a polynomial  a + b x + c x^2 + d x^3.
  CubicPolynomial(T a, T b, T c, T d) : coefficients_{{a, b, c, d}} { }
  // Create a polynomial from an array.
  CubicPolynomial(std::array<T, 4> const& coefficients) : coefficients_(coefficients) { }

  // Construct a cubic polynomial from a Polynomial.
  CubicPolynomial(Polynomial<T> const& polynomial)
  {
    std::vector<T> const& coefficients = polynomial.coefficients();
    // polynomial must have degree three or less.
    ASSERT(coefficients.size() == coefficients_.size());
    std::copy(coefficients.begin(), coefficients.end(), coefficients_.begin());
  }

  void initialize(T x0, T y0, T dxdy0, T x1, T y1, T dxdy1);
  int get_extrema(std::array<T, 2>& extrema_out, bool left_most_first = true) const;

  // Return the derivative of this cubic.
  QuadraticPolynomial<T> derivative() const
  {
    return {coefficients_[1], 2 * coefficients_[2], 3 * coefficients_[3]};
  }

  // Evaluation.
  T operator()(T x) const
  {
    return coefficients_[0] + (coefficients_[1] + (coefficients_[2] + coefficients_[3] * x) * x) * x;
  };

  T evaluate(T x) const
  {
    return coefficients_[0] + (coefficients_[1] + (coefficients_[2] + coefficients_[3] * x) * x) * x;
  };

  // Evaluate derivative.
  T derivative(T x) const
  {
    return coefficients_[1] + (2 * coefficients_[2] + 3 * coefficients_[3] * x) * x;
  }

  // Evaluate half times the second derivative.
  T half_second_derivative(T x) const
  {
    return coefficients_[2] + 3 * coefficients_[3] * x;
  }

  // Evaluate second derivative.
  T second_derivative(T x) const
  {
    return 2 * half_second_derivative(x);
  }

  T inflection_point() const
  {
    if (AI_UNLIKELY(coefficients_[3] == 0))
      return std::numeric_limits<T>::infinity();           // Anything really large would do (or towards negative infinity too).
    return -coefficients_[2] / (3 * coefficients_[3]);
  }

  // Access coefficients.
  T operator[](int i) const { ASSERT(0 <= i && i < coefficients_.size()); return coefficients_[i]; }
  T& operator[](int i) { ASSERT(0 <= i && i < coefficients_.size()); return coefficients_[i]; }

  CubicPolynomial& operator-=(CubicPolynomial const& rhs)
  {
    for (int i = 0; i < coefficients_.size(); ++i)
      coefficients_[i] -= rhs.coefficients_[i];
    return *this;
  }

  CubicPolynomial& operator+=(CubicPolynomial const& rhs)
  {
    for (int i = 0; i < coefficients_.size(); ++i)
      coefficients_[i] += rhs.coefficients_[i];
    return *this;
  }

  CubicPolynomial operator+(CubicPolynomial const& rhs) const
  {
    CubicPolynomial result(*this);
    result += rhs;
    return result;
  }

  CubicPolynomial operator-(CubicPolynomial const& rhs) const
  {
    CubicPolynomial result(*this);
    result -= rhs;
    return result;
  }

  CubicPolynomial& operator*=(T factor)
  {
    for (int i = 0; i < coefficients_.size(); ++i)
      coefficients_[i] *= factor;
    return *this;
  }

  friend CubicPolynomial operator*(T factor, CubicPolynomial const& rhs)
  {
    CubicPolynomial result(rhs);
    result *= factor;
    return result;
  }

  // Return the division of this Polynomial by the factor (x - r).
  QuadraticPolynomial<T> long_division(T r, T& remainder) const
  {
    QuadraticPolynomial<T> result;
    result[2] = coefficients_[3];
    result[1] = coefficients_[2] + r * result[2];
    result[0] = coefficients_[1] + r * result[1];
    remainder = coefficients_[0] + r * result[0];
    return result;
  }

  int get_roots(std::array<T, 3>& roots_out) const;

#if CW_DEBUG
  // Return true if one was assigned from the other.
  friend bool operator==(CubicPolynomial const& lhs, CubicPolynomial const& rhs)
  {
    return lhs.coefficients_ == rhs.coefficients_;
  }

  void print_on(std::ostream& os) const;
#endif
};

template<typename T>
requires ConceptMathOperations<T>
void CubicPolynomial<T>::initialize(T x0, T y0, T dxdy0, T x1, T y1, T dxdy1)
{
  T delta_x_inverse = 1 / (x0 - x1);

  // The theory of this approach is described here:
  // https://math.stackexchange.com/a/4926903/489074
  T d = (dxdy0 + dxdy1 - 2 * (y0 - y1) * delta_x_inverse) * (delta_x_inverse * delta_x_inverse);
  T c = ((dxdy0 - dxdy1) * delta_x_inverse - 3 * d * (x0 + x1)) / 2;
  T b = (x0 * dxdy1 - x1 * dxdy0) * delta_x_inverse + 3 * d * x0 * x1;

  coefficients_[3] = d;
  coefficients_[2] = c;
  coefficients_[1] = b;
  coefficients_[0] = 0;
  coefficients_[0] = y0 - operator()(x0);
}

template<typename T>
requires ConceptMathOperations<T>
int CubicPolynomial<T>::get_extrema(std::array<T, 2>& extrema_out, bool left_most_first) const
{
  DoutEntering(dc::notice, "CubicPolynomial<" << libcwd::type_info_of<T>().demangled_name() <<
      ">::get_extrema(extrema_out, left_most_first = " << std::boolalpha << left_most_first << ")");

  // The cubic is:
  // coefficients_[0] + coefficients_[1] * x + coefficients_[2] * x^2 + coefficients_[3] * x^3.
  // The derivative is:
  // coefficients_[1] + 2 * coefficients_[2] * x + 3 * coefficients_[3] * x^2.

  if (coefficients_[3] == 0)
  {
    // If coefficients_[3] is zero, then the derivative has a single root when
    // coefficients_[1] + 2 * coefficients_[2] * x = 0 -->
    // x = -coefficients_[1] / (2 * coefficients_[2]).
    extrema_out[0] = coefficients_[1] / (-2 * coefficients_[2]);
    return isfinite(extrema_out[0]) ? 1 : 0;
  }

  // The determinant is (2 * coefficients_[2])^2 - 4 * coefficients_[1] * (3 * coefficients_[3]).
  T one_fourth_D = utils::square(coefficients_[2]) - 3 * coefficients_[1] * coefficients_[3];
  Dout(dc::notice, "D = " << (4 * one_fourth_D));

  if (one_fourth_D < 0)
  {
    // The inflection point is
    // -(2 * coefficients_[2]) / (2 * (3 * coefficients_[3])) =
    // -coefficients_[2] / (3 * coefficients_[3]).

    // Write the inflection point to index 0.
    extrema_out[0] = -coefficients_[2] / (3 * coefficients_[3]);
    return 0;
  }

  // Use a sqrt with the same sign as coefficients_[2];
  T const signed_half_sqrt_D = copysign(sqrt(one_fourth_D), coefficients_[2]);
  // The first root that is calculated should go here (see brute_force_index_first_root.cxx).
  int index_first_root = (left_most_first ? (coefficients_[3] > 0) : false) == (copysign(T{1}, coefficients_[2]) > 0);

  // Calculate the root closest to zero.
  extrema_out[index_first_root] = -coefficients_[1] / (coefficients_[2] + signed_half_sqrt_D);

  if (AI_UNLIKELY(isnan(extrema_out[0])))
  {
    // This means we must have divided by zero, which means that both, coefficients_[1] as well as sqrtD, must be zero.
    // The latter means that coefficients_[0] is zero (coefficients_[2] was already checked not to be zero).
    // Therefore we have: f(x) = c x^2 with one root at x=0.
    extrema_out[0] = 0;
    return 1;
  }

  // Calculate the root further away from zero.
  extrema_out[1 - index_first_root] = -(coefficients_[2] + signed_half_sqrt_D) / (3 * coefficients_[3]);

  Dout(dc::notice, "extrema_out = " << std::setprecision(std::numeric_limits<T>::digits10) << extrema_out);

  // The smallest one must be in index 0 if left_most_first is true.
  ASSERT(one_fourth_D == 0 || !left_most_first || extrema_out[0] <= extrema_out[1]);
  // The minimum must be in index 0 if left_most_first is false.
  ASSERT(one_fourth_D == 0 || left_most_first || (extrema_out[0] < extrema_out[1] == coefficients_[3] < 0));

  return (one_fourth_D == 0) ? 1 : 2;
}

template<typename T>
requires ConceptMathOperations<T>
int CubicPolynomial<T>::get_roots(std::array<T, 3>& roots_out) const
{
  DoutEntering(dc::cubic, "CubicPolynomial<" << libcwd::type_info_of<T>().demangled_name() << ">::get_roots() for " << *this);

  // Include the body of the function.
# include "CubicPolynomial_get_roots.cpp"
}

#ifdef CWDEBUG
template<typename T>
requires ConceptMathOperations<T>
void CubicPolynomial<T>::print_on(std::ostream& os) const
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

} // namespace math
