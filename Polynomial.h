#pragma once

#include "ConceptMathOperations.h"
#include "utils/macros.h"
#include <vector>
#include <array>
#include <ranges>
#include <complex>
#include "debug.h"
#ifdef CWDEBUG
#include "utils/has_print_on.h"
#include <string>
#endif

// This should become part of machine-learning in the end.
namespace math {
#ifdef CWDEBUG
using utils::has_print_on::operator<<;
#endif

template<typename T>
requires ConceptMathOperations<T>
class Polynomial
{
 protected:
  std::vector<T> coefficients_;
#ifdef CWDEBUG
  std::string symbol_name_;
#endif

 public:
  // Create a polynomial of at most degree `number_of_coefficients - 1`, starting
  // with all coefficients set to zero.
  Polynomial(int number_of_coefficients COMMA_CWDEBUG_ONLY(std::string const& symbol_name)) :
    coefficients_(number_of_coefficients) COMMA_CWDEBUG_ONLY(symbol_name_(symbol_name)) { }

  T operator[](int i) const { ASSERT(0 <= i && i < coefficients_.size()); return coefficients_[i]; }
  T& operator[](int i) { ASSERT(0 <= i && i < coefficients_.size()); return coefficients_[i]; }

  // Accessor.
  std::vector<T> const& coefficients() const { return coefficients_; }

  int actual_degree() const
  {
    int degree = coefficients_.size() - 1;
    while (degree > 0 && coefficients_[degree] == 0.0)
      --degree;
    return degree;
  }

  T operator()(T w) const
  {
    T result = 0.0;
    for (T coefficient : std::ranges::reverse_view(coefficients_))
      result = w * result + coefficient;
    return result;
  };

  Polynomial& operator-=(T rhs)
  {
    coefficients_[0] -= rhs;
    return *this;
  }

  Polynomial& operator+=(T rhs)
  {
    coefficients_[0] += rhs;
    return *this;
  }

  Polynomial operator-(T rhs) const
  {
    Polynomial result(*this);
    result -= rhs;
    return result;
  }

  Polynomial operator+(T rhs) const
  {
    Polynomial result(*this);
    result += rhs;
    return result;
  }

  Polynomial& operator-=(Polynomial const& rhs)
  {
    for (int i = 0; i < coefficients_.size(); ++i)
      coefficients_[i] -= rhs.coefficients_[i];
    return *this;
  }

  Polynomial& operator+=(Polynomial const& rhs)
  {
    for (int i = 0; i < coefficients_.size(); ++i)
      coefficients_[i] += rhs.coefficients_[i];
    return *this;
  }

  Polynomial operator+(Polynomial const& rhs) const
  {
    Polynomial result(*this);
    result += rhs;
    return result;
  }

  Polynomial operator-(Polynomial const& rhs) const
  {
    Polynomial result(*this);
    result -= rhs;
    return result;
  }

  Polynomial& operator*=(T factor)
  {
    for (int i = 0; i < coefficients_.size(); ++i)
      coefficients_[i] *= factor;
    return *this;
  }

  friend Polynomial operator*(T factor, Polynomial const& rhs)
  {
    Polynomial result(rhs);
    result *= factor;
    return result;
  }

  Polynomial operator*=(Polynomial const& rhs);
  PRAGMA_DIAGNOSTIC_PUSH_IGNORE("-Wnon-template-friend")
  friend Polynomial operator*(Polynomial const& lhs, Polynomial const& rhs);
  PRAGMA_DIAGNOSTIC_POP

  Polynomial derivative() const
  {
    Polynomial result(coefficients_.size() - 1 COMMA_CWDEBUG_ONLY(symbol_name_));
    for (int i = 1; i < coefficients_.size(); ++i)
      result[i - 1] = coefficients_[i] * i;
    return result;
  }

  Polynomial integrate() const
  {
    Polynomial result(coefficients_.size() + 1 COMMA_CWDEBUG_ONLY(symbol_name_));
    for (int i = 0; i < coefficients_.size(); ++i)
      result[i + 1] = coefficients_[i] / (i + 1);
    return result;
  }

  // Return the division of this Polynomial by the factor (w - z).
  Polynomial long_division(T z, T& remainder) const
  {
    // f(w) = 3 * w^3 +      5 * w^2 - 4 * w + 10.
    //        3 * w^3 + (-2)*3 * w^2
    //      - ------------------------------------
    //                      11 * w^2 -     4 * w + 10.
    //                      11 * w^2 + (-2)*11 w
    //                    - --------------------------
    //                                    18 * w +     10.
    //                                    18 * w + (-2)18
    //                                  - ---------------
    //                                                 46
    // Divide by (w - 2)
    // 3 * w^2 + 11 * w + 18

    // NOTICE        : 10 + -4 w + 5 w^2 + 3 w^3
    // (w - 2)(3 w^2 + 11 w + 18) = 3 w^3 + 5 w^2 - 4 w - 36

    if (coefficients_.size() < 2)
    {
      ASSERT(coefficients_.size() == 1);
      remainder = coefficients_[0];
      return {1 COMMA_CWDEBUG_ONLY(symbol_name_)};
    }
    Polynomial result(coefficients_.size() - 1 COMMA_CWDEBUG_ONLY(symbol_name_));
    result[coefficients_.size() - 2] = coefficients_[coefficients_.size() - 1];
    for (int i  = coefficients_.size() - 2; i > 0; --i)
      result[i - 1] = coefficients_[i] + z * result[i];
    remainder = coefficients_[0] + z * result[0];
    return result;
  }

  // Get the real roots of a polynomial, up till degree two.
  int get_roots(std::array<T, 2>& roots_out) const;

  // Get the complex roots of a polynomial, up till degree five.
  // Returns the number of roots (equal to the degree of the Polynomial).
  int get_roots(std::array<std::complex<T>, 5>& roots_out) const;

#ifdef CWDEBUG
  void print_on(std::ostream& os) const;
#endif
};

} // namespace math

#include <Eigen/Dense>
#include "utils/macros.h"
#include "utils/square.h"

namespace math {

template<typename T>
requires ConceptMathOperations<T>
Polynomial<T> Polynomial<T>::operator*=(Polynomial<T> const& rhs)
{
  int const ad_lhs = actual_degree();
  int const ad_rhs = rhs.actual_degree();
  size_t const number_of_coefficients = 1 + ad_lhs + ad_rhs;
  std::vector<T> coefficients_lhs(std::max({number_of_coefficients, coefficients_.size(), rhs.coefficients().size()}), 0.0);
  std::swap(coefficients_, coefficients_lhs);
  for (int d_lhs = 0; d_lhs <= ad_lhs; ++d_lhs)
    for (int d_rhs = 0; d_rhs <= ad_rhs; ++d_rhs)
      coefficients_[d_lhs + d_rhs] += coefficients_lhs[d_lhs] * rhs[d_rhs];
  return *this;
}

template<typename T>
requires ConceptMathOperations<T>
Polynomial<T> operator*(Polynomial<T> const& lhs, Polynomial<T> const& rhs)
{
  int const ad_lhs = lhs.actual_degree();
  int const ad_rhs = rhs.actual_degree();
  size_t const number_of_coefficients = 1 + ad_lhs + ad_rhs;
  Polynomial<T> result(std::max({number_of_coefficients, lhs.coefficients().size(), rhs.coefficients().size()})
      COMMA_CWDEBUG_ONLY(lhs.symbol_name_));

  for (int d_lhs = 0; d_lhs <= ad_lhs; ++d_lhs)
    for (int d_rhs = 0; d_rhs <= ad_rhs; ++d_rhs)
      result[d_lhs + d_rhs] += lhs[d_lhs] * rhs[d_rhs];

  return result;
}

template<typename T>
requires ConceptMathOperations<T>
int Polynomial<T>::get_roots(std::array<T, 2>& roots_out) const
{
  // This can be at most a parabola.
  ASSERT(1 <= coefficients_.size() && coefficients_.size() <= 3);
  if (coefficients_.size() < 3 || coefficients_[2] == 0.0)
  {
    if (coefficients_.size() < 2)
      return 0;
    roots_out[0] = -coefficients_[0] / coefficients_[1];
    return isfinite(roots_out[0]) ? 1 : 0;
  }

  T const D = utils::square(coefficients_[1]) - 4.0 * coefficients_[2] * coefficients_[0];
  if (D < 0.0)
    return 0;
  // Use a sqrt with the same sign as coefficients_[1];
  T const signed_sqrt_D = copysign(sqrt(D), coefficients_[1]);

  // Calculate the root closest to zero.
  roots_out[0] = -2.0 * coefficients_[0] / (coefficients_[1] + signed_sqrt_D);

  if (AI_UNLIKELY(isnan(roots_out[0])))
  {
    // This means we must have divided by zero, which means that both, coefficients_[1] as well as sqrtD, must be zero.
    // The latter means that coefficients_[0] is zero (coefficients_[2] was already checked not to be zero).
    // Therefore we have: f(x) = c x^2 with one root at x=0.
    roots_out[0] = 0.0;
    return 1;
  }

  // Calculate the root further away from zero.
  roots_out[1] = -0.5 * (coefficients_[1] + signed_sqrt_D) / coefficients_[2];

  // The second one is larger in absolute value.
  ASSERT(abs(roots_out[1]) > abs(roots_out[0]));

  return 2;
}

template<typename T>
requires ConceptMathOperations<T>
int Polynomial<T>::get_roots(std::array<std::complex<T>, 5>& roots_out) const
{
  int degree = coefficients_.size() - 1;

  // This function only works up till degree 5.
  ASSERT(degree <= 5);

  // Abbreviate the coefficients as a[0] etc.
  auto const& a = coefficients_;

  if (degree == 5 && a[5] != 0.0)
  {
    // Construct the 5x5 companion matrix.
    Eigen::Matrix<T, 5, 5> C;
    C <<  0, 0, 0, 0, -a[0] / a[5],
          1, 0, 0, 0, -a[1] / a[5],
          0, 1, 0, 0, -a[2] / a[5],
          0, 0, 1, 0, -a[3] / a[5],
          0, 0, 0, 1, -a[4] / a[5];

    // Compute the eigenvalues of the companion matrix.
    Eigen::EigenSolver<Eigen::Matrix<T, 5, 5>> solver(C);
    Eigen::Matrix<std::complex<T>, 5, 1> roots = solver.eigenvalues();

    // Store the roots in the output array.
    for (int i = 0; i < roots.size(); ++i)
        roots_out[i] = roots[i];
  }
  else if (degree >= 4 && a[4] != 0.0)
  {
    // Construct the 4x4 companion matrix.
    Eigen::Matrix4d C;
    C <<  0, 0, 0, -a[0] / a[4],
          1, 0, 0, -a[1] / a[4],
          0, 1, 0, -a[2] / a[4],
          0, 0, 1, -a[3] / a[4];

    // Compute the eigenvalues of the companion matrix.
    Eigen::EigenSolver<Eigen::Matrix4d> solver(C);
    Eigen::Vector4cd roots = solver.eigenvalues();

    // Store the roots in the output array.
    for (int i = 0; i < roots.size(); ++i)
      roots_out[i] = roots[i];
  }
  else if (degree >= 3 && a[3] != 0.0)
  {
    degree = 3;

    // Construct the 3x3 companion matrix.
    Eigen::Matrix3d C;
    C <<  0, 0, -a[0] / a[3],
          1, 0, -a[1] / a[3],
          0, 1, -a[2] / a[3];

    // Compute the eigenvalues of the companion matrix.
    Eigen::EigenSolver<Eigen::Matrix3d> solver(C);
    Eigen::Vector3cd roots = solver.eigenvalues();

    // Store the roots in the output array.
    for (int i = 0; i < roots.size(); ++i)
      roots_out[i] = roots[i];
  }
  else if (degree >= 2 && a[2] != 0.0)
  {
    degree = 2;

    // Construct the 2x2 companion matrix.
    Eigen::Matrix2d C;
    C <<  0, -a[0] / a[2],
          1, -a[1] / a[2];

    // Compute the eigenvalues of the companion matrix.
    Eigen::EigenSolver<Eigen::Matrix2d> solver(C);
    Eigen::Vector2cd roots = solver.eigenvalues();

    // Store the roots in the output array.
    for (int i = 0; i < roots.size(); ++i)
      roots_out[i] = roots[i];
  }
  else if (degree >= 1 && a[1] != 0.0)
  {
    degree = 1;

    roots_out[0] = -a[0] / a[1];
  }
  else
    degree = 0;

  return degree;
}

#ifdef CWDEBUG
template<typename T>
requires ConceptMathOperations<T>
void Polynomial<T>::print_on(std::ostream& os) const
{
  bool first = true;
  int exponent = 0;
  for (T coefficient : coefficients_)
  {
    if (coefficient != 0.0)
    {
      if (first)
        os << coefficient;
      else if (coefficient > 0.0)
        os << " + " << coefficient;
      else
        os << " - " << -coefficient;
      if (exponent > 0)
      {
        os << ' ' << symbol_name_;
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
