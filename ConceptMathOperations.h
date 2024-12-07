#pragma once

#include <cmath>
#include <concepts>

namespace math {

// Define the used mathematical operators for builtin type arguments in namespace math.
template<std::floating_point T>
[[gnu::always_inline]] inline T abs(T a)
{
  return std::abs(a);
}

template<std::floating_point T>
[[gnu::always_inline]] inline T sqrt(T a)
{
  return std::sqrt(a);
}

template<std::floating_point T>
[[gnu::always_inline]] inline T copysign(T mag, T sgn)
{
  return std::copysign(mag, sgn);
}

template<std::floating_point T>
[[gnu::always_inline]] inline bool isfinite(T a)
{
  return std::isfinite(a);
}

template<std::floating_point T>
[[gnu::always_inline]] inline bool isnan(T a)
{
  return std::isnan(a);
}

template<typename T, typename U>
concept ConceptSameAsOrConst = std::same_as<T, U> || std::same_as<T, std::add_const_t<U>>;

// Concept to check for required operations.
template<typename T>
concept ConceptMathOperations = requires(T a) {
  { abs(a) } -> ConceptSameAsOrConst<T>;
  { sqrt(a) } -> ConceptSameAsOrConst<T>;
  { copysign(a, a) } -> ConceptSameAsOrConst<T>;
  { isfinite(a) } -> std::convertible_to<bool>;
  { isnan(a) } -> std::convertible_to<bool>;
};

} // namespace math
