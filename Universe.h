#pragma once

#include "utils/uint_leastN_t.h"
#include "utils/BitSet.h"
#include "utils/macros.h"
#include "utils/create_mask.h"
#include <Eigen/Core>
#include <type_traits>
#include <bit>
#ifdef CWDEBUG
#include <utils/has_print_on.h>
#endif
#include "debug.h"

namespace math {
#ifdef CWDEBUG
using utils::has_print_on::operator<<;
#endif

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

  // A BitSet with at least as many bits as the number of axis of the universe.
  using axes_type = utils::BitSet<utils::uint_leastN_t<detail::universe_max_n_v<U>>>;
  // The POD type for that BitSet.
  using subset_type = utils::BitSetPOD<typename axes_type::mask_type>;

 private:
  float_type scale_factor_;                     // A unit vector in this Basis is scale_factor_ times the unit vector in Universe coordinates pointing the same way.
                                                // That means that the same abstract vector v has a norm that is scale_factor_ smaller expressed in this Basis than in Universe coordinates.
  rotation_matrix_type rotation_matrix_;        // Proper rotation that aligns the Universe standard basis with this Basis (first n axes retained).

 public:
  Basis() : scale_factor_(1), rotation_matrix_(rotation_matrix_type::Identity()) { }

  Basis(float_type scale_factor) : scale_factor_(scale_factor), rotation_matrix_(rotation_matrix_type::Identity()) { }

 private:
  friend U;
  // Construct a Basis with a rotation matrix that is all zeroes.
  Basis(nullptr_t) : scale_factor_(1), rotation_matrix_(rotation_matrix_type::Zero()) { }

#if CWDEBUG
 public:
  void print_on(std::ostream& os) const;
#endif
};

// End of Basis
//=============================================================================

template<int N>
struct Permutation
{
  std::array<uint8_t, N> seq_{};

  Permutation(std::array<uint8_t, N> const& seq) : seq_(seq) { }

  template<std::integral Index>
  Permutation(std::initializer_list<Index> perm)
  {
    ASSERT(static_cast<int>(perm.size()) == N);
    auto it = perm.begin();
    for (int i = 0; i < N; ++i, ++it)
      seq_[i] = static_cast<uint8_t>(*it);
  }

  bool is_odd() const
  {
    int count = 0;
    uint64_t mask = 0;
    for (uint8_t e : seq_)
    {
      uint64_t em = uint64_t{1} << e;
      if ((mask & em))
        ++count;
      mask ^= (em - 1);
    }
    return (count & 1);
  }

  auto begin() const { return seq_.begin(); }
  auto end() const { return seq_.end(); }

#if CW_DEBUG
  // Check if the Permutation is valid for the given range.
  bool in_range(uint8_t range_start, uint8_t range_end) const
  {
    for (uint8_t e : seq_)
      if (!(range_start <= e && e < range_end))
        return false;
    return true;
  }
#endif
};

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
  using axes_mask_type = utils::uint_leastN_t<max_n>;
  using axes_type = basis_type::axes_type;
  using subset_type = basis_type::subset_type;

  static basis_type standard_basis;

  template<int N>
  struct CoordinateSubspace
  {
    static_assert(N <= MAX_N, "The CoordinateSubspace can't have more dimensions than its Universe.");
    using basis_type = Basis<Universe, N>;

    // Return a basis from a permutation.
    // If the permutation is odd, flip the sign of the basis axis at flip_output_axis_if_odd (default 0) to restore orientation.
    static basis_type from_permutation(Permutation<N> permutation, T scale = 1, int flip_output_axis_if_odd = 0);
  };
};

//static
template<typename ID, int MAX_N, typename T>
Basis<Universe<ID, MAX_N, T>> Universe<ID, MAX_N, T>::standard_basis;

//static
template<typename ID, int MAX_N, typename T>
template<int N>
Universe<ID, MAX_N, T>::CoordinateSubspace<N>::basis_type
Universe<ID, MAX_N, T>::CoordinateSubspace<N>::from_permutation(Permutation<N> permutation, T scale, int flip_output_axis_if_odd)
{
  Universe<ID, MAX_N, T>::CoordinateSubspace<N>::basis_type basis{nullptr};
  basis.scale_factor_ = scale;

  // Track which Universe axes are still free.
  axes_mask_type all_axes_mask = utils::create_mask<axes_mask_type, max_n>();
  axes_type remaining{all_axes_mask};
  using Index = typename axes_type::Index;

  ASSERT(permutation.in_range(0, MAX_N));
  ASSERT(flip_output_axis_if_odd >= 0 && flip_output_axis_if_odd < N);

  bool const permutation_is_odd = permutation.is_odd();
  auto col = permutation.begin();
  for (int row = 0; row < N; ++row, ++col)
  {
    ASSERT(*col < max_n);
    Index axis{utils::bitset::IndexPOD{static_cast<int8_t>(*col)}};
    ASSERT(remaining.test(axis));          // Must be unused so far.
    basis.rotation_matrix_(row, *col) = (permutation_is_odd && row == flip_output_axis_if_odd) ? -1 : 1;
    remaining.reset(axis);                 // Mark axis as used.
  }

  // Fill remaining rows with the unused Universe axes, ascending.
  int row = N;
  for (auto bit : remaining)
  {
    ASSERT(row < max_n);
    Index axis = axes_type::mask2index(bit());
    basis.rotation_matrix_(row++, axis()) = 1;
  }

  return basis;
}

// End of Universe
//=============================================================================

} // namespace math

#if CWDEBUG

template<math::ConceptUniverse U, int N>
void math::Basis<U, N>::print_on(std::ostream& os) const
{
  os << "{scale_factor:" << scale_factor_ << ", rotation_matrix:\n" << rotation_matrix_ << '}';
}
#endif
