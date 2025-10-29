#pragma once

#include "utils/uint_leastN_t.h"
#include "utils/BitSet.h"
#include "utils/macros.h"
#include "utils/create_mask.h"
#include <type_traits>
#include  <bit>

namespace math {

// Forward declare Universe.
template<typename ID, int MAX_N>
struct Universe;

//-----------------------------------------------------------------------------
// ConceptUniverse
//

namespace detail {
// Type trait for Universe.

template<typename>
struct is_universe : std::false_type { };

template<typename ID, int MAX_N>
struct is_universe<Universe<ID, MAX_N>> : std::true_type { };

template<typename T>
constexpr bool is_universe_v = is_universe<T>::value;

} // namespace detail

template<typename T>
concept ConceptUniverse = detail::is_universe_v<std::remove_cvref_t<T>>;

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

template<typename ID, int MAX_N>
struct universe_max_n<Universe<ID, MAX_N>> : std::integral_constant<int, MAX_N> { };

// Optional helper variable:
template<typename T>
constexpr int universe_max_n_v = universe_max_n<T>::value;

} // namespace detail

template<ConceptUniverse U, auto used_dimensions = utils::create_mask<utils::uint_leastN_t<detail::universe_max_n_v<U>>, detail::universe_max_n_v<U>>()>
class Basis
{
 public:
  using basis_type = Basis;
  static constexpr int n = std::popcount(used_dimensions);
  using axes_type = utils::BitSet<utils::uint_leastN_t<detail::universe_max_n_v<U>>>;
  using subset_type = utils::BitSetPOD<typename axes_type::mask_type>;

 private:
  float_type dilation;
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
template<typename ID, int MAX_N>
struct Universe
{
  static constexpr int max_n = MAX_N;
  using axes_mask_type = utils::uint_leastN_t<max_n>;
  using basis_type = Basis<Universe>;
  using axes_type = basis_type::axes_type;
  using subset_type = basis_type::subset_type;

  static basis_type standard_basis;

  template<auto BitSetPOD>
  struct CoordinateSubspace
  {
    static constexpr auto axes = BitSetPOD.m_bitmask;
PRAGMA_DIAGNOSTIC_PUSH_IGNORE("-Wchanges-meaning")
    using Basis = Basis<Universe, axes>;
PRAGMA_DIAGNOSTIC_POP
  };
};

//static
template<typename ID, int MAX_N>
Basis<Universe<ID, MAX_N>> Universe<ID, MAX_N>::standard_basis;

// End of Universe
//=============================================================================

} // namespace math
