#pragma once

#include "utils/uint_leastN_t.h"
#include "utils/BitSet.h"
#include "utils/macros.h"
#include "utils/create_mask.h"
#include <Eigen/Core>
#include <type_traits>
#include <bit>

namespace math {

// Forward declare Universe.
template<typename ID, int MAX_N, typename T>
struct Universe;

//-----------------------------------------------------------------------------
// ConceptUniverse
//

namespace detail {
// Type trait for Universe.

template<typename>
struct is_universe : std::false_type { };

template<typename ID, int MAX_N, typename T>
struct is_universe<Universe<ID, MAX_N, T>> : std::true_type { };

template<typename U>
constexpr bool is_universe_v = is_universe<U>::value;

} // namespace detail

template<typename U>
concept ConceptUniverse = detail::is_universe_v<std::remove_cvref_t<U>>;

// End of ConceptUniverse
//-----------------------------------------------------------------------------

//=============================================================================
// Basis
//
// An affine subspace of n dimensions.
//
// used_dimensions : a bit mask (one bit for each dimension of U) where a set
//                   bit means that that dimension is part of this basis.

namespace detail {

// Trait to extract MAX_N.
template<typename> struct universe_max_n;

template<typename ID, int MAX_N, typename T>
struct universe_max_n<Universe<ID, MAX_N, T>> : std::integral_constant<int, MAX_N> { };

// Optional helper variable:
template<typename U>
constexpr int universe_max_n_v = universe_max_n<U>::value;

// Trait to extract T.
template<typename> struct universe_float_type;

template<typename ID, int MAX_N, typename T>
struct universe_float_type<Universe<ID, MAX_N, T>> { using type = T; };

template<typename U>
using universe_float_type_t = typename universe_float_type<U>::type;

} // namespace detail

template<ConceptUniverse U, int N = detail::universe_max_n_v<U>>
class Basis
{
 public:
  static constexpr int n = N;
  using basis_type = Basis;
  using float_type = detail::universe_float_type_t<U>;
  using rotation_matrix_type = Eigen::Matrix<float_type, detail::universe_max_n_v<U>, detail::universe_max_n_v<U>>;
//  using axes_type = utils::BitSet<utils::uint_leastN_t<detail::universe_max_n_v<U>>>;
//  using subset_type = utils::BitSetPOD<typename axes_type::mask_type>;

 private:
  float_type scale_factor_;                     // A unit vector in this Basis is scale_factor_ times the unit vector in Universe coordinates pointing the same way.
                                                // That means that the same abstract vector v has a norm that is scale_factor_ smaller expressed in this Basis than in Universe coordinates.
  rotation_matrix_type rotation_matrix_;        // Proper rotation that aligns the Universe standard basis with this Basis (first n axes retained).

 public:
  Basis() : scale_factor_(1), rotation_matrix_(rotation_matrix_type::Identity()) { }

  Basis(float_type scale_factor) : scale_factor_(scale_factor), rotation_matrix_(rotation_matrix_type::Identity()) { }
};

// End of Basis
//=============================================================================

//=============================================================================
// Universe
//
// An affine space of MAX_N dimensions, that is not a subspace of any other space.
//
// Two Universe objects with different ID are completely independent.
// They have each a basis that can not be compared: there exists no transformation
// between the two basis, not even if they have the same dimension.
//
template<typename ID, int MAX_N, typename T = double>
struct Universe
{
  using float_type = T;
  static constexpr int max_n = MAX_N;
  using basis_type = Basis<Universe>;
  using rotation_matrix_type = basis_type::rotation_matrix_type;

//  using axes_mask_type = utils::uint_leastN_t<max_n>;
//  using axes_type = basis_type::axes_type;
//  using subset_type = basis_type::subset_type;

  static basis_type standard_basis;

  template<int N>
  struct CoordinateSubspace
  {
    static_assert(N <= MAX_N, "The CoordinateSubspace can't have more dimensions than its Universe.");
//    static constexpr auto axes = BitSetPOD.m_bitmask;
    using basis_type = Basis<Universe, N /*, axes*/>;
  };
};

//static
template<typename ID, int MAX_N, typename T>
Basis<Universe<ID, MAX_N, T>> Universe<ID, MAX_N, T>::standard_basis;

// End of Universe
//=============================================================================

} // namespace math
