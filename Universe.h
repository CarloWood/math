#pragma once

#include "math/Vector.h"
#include "utils/uint_leastN_t.h"
#include "utils/BitSet.h"
#include "utils/macros.h"
#include "utils/create_mask.h"
#include "utils/almost_equal.h"
#include "utils/AIAlert.h"
#include "utils/print_range.h"
#include <Eigen/Core>
#include <type_traits>
#include <bit>
#include <limits>
#ifdef CWDEBUG
#include <utils/has_print_on.h>
#endif
#include "debug.h"

namespace math {
#ifdef CWDEBUG
using utils::has_print_on::operator<<;
#endif

// Forward declare Universe.
template<typename ID, size_t MAX_N, typename T>
struct Universe;

//-----------------------------------------------------------------------------
// ConceptUniverse
//

namespace detail {
// Type trait for Universe.

template<typename>
struct is_universe : std::false_type { };

template<typename ID, size_t MAX_N, typename T>
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
// An affine subspace of N dimensions.
//

namespace detail {

// Trait to extract MAX_N.
template<typename> struct universe_max_n;

template<typename ID, size_t MAX_N, typename T>
struct universe_max_n<Universe<ID, MAX_N, T>> : std::integral_constant<size_t, MAX_N> { };

// Optional helper variable:
template<typename U>
constexpr size_t universe_max_n_v = universe_max_n<U>::value;

// Trait to extract T.
template<typename> struct universe_float_type;

template<typename ID, size_t MAX_N, typename T>
struct universe_float_type<Universe<ID, MAX_N, T>> { using type = T; };

template<typename U>
using universe_float_type_t = typename universe_float_type<U>::type;

template<ConceptUniverse U, int N>
class BasisBuilder;

} // namespace detail

template<ConceptUniverse U, size_t N = detail::universe_max_n_v<U>>
class Basis
{
 public:
  static constexpr int n = N;
  static constexpr int max_n = detail::universe_max_n_v<U>;
  using basis_type = Basis;
  using float_type = detail::universe_float_type_t<U>;
  using vector_type = Vector<max_n, float_type>;
  using rotation_matrix_type = Eigen::Matrix<float_type, max_n, max_n>;

  // A BitSet with at least as many bits as the number of axis of the universe.
  using axes_type = utils::BitSet<utils::uint_leastN_t<max_n>>;
  // The POD type for that BitSet.
  using subset_type = utils::BitSetPOD<typename axes_type::mask_type>;

 private:
  // A unit vector in this Basis is scale_factor_ times the unit vector in Universe
  // coordinates pointing the same way. That means that the same abstract vector v
  // has a norm that is scale_factor_ smaller expressed in this Basis than in
  // Universe coordinates.
  float_type scale_factor_;

  // Proper rotation that converts Universe standard coordinates to this Basis:
  // for any Universe vector v (in standard coordinates), rotation_matrix_ * v
  // yields the coordinates of v in the Basis.
  //
  // For example, in a 5D universe containing a 3D subspace that is spanned by
  // the (5D) orthogonal unit vectors x̂, ŷ and ẑ (in universe coordinates),
  // rotation_matrix_ would be:
  //
  //                          ⎛ <⋯ row unit vector x̂ ⋯> ⎞       ⎫
  //                          ⎜ <⋯ row unit vector ŷ ⋯> ⎟       ⎬ The first N rows are the N basis axes expressed as unit vectors in Universe coordinates.
  //   rotation_matrix_ (R) = ⎜ <⋯ row unit vector ẑ ⋯> ⎟       ⎭
  //                          ⎜       orthonormal       ⎟       ⎞ The remaining rows form an orthonormal completion of the Universe.
  //                          ⎝       completion        ⎠       ⎠
  //
  // Such that `Rx̂ = (1, 0, 0, 0, 0)`, `Rŷ = (0, 1, 0, 0, 0)` and `Rẑ = (0, 0, 1, 0, 0)`.
  // In other words, it expresses vectors in Universe coordinates as vectors in Basis coordinates.
  rotation_matrix_type rotation_matrix_;

 public:
  Basis() : scale_factor_(1), rotation_matrix_(rotation_matrix_type::Identity()) { }

  Basis(float_type scale_factor) : scale_factor_(scale_factor), rotation_matrix_(rotation_matrix_type::Identity()) { }

  Basis(detail::BasisBuilder<U, N> builder);

  // Apply the unique proper rotation that maps v1 to v2 by a rotation in the plane
  // spanned by v1 and v2. The orthogonal complement is left invariant.
  // If v1·v2 < 0 then v2 is first negated so that the rotation angle is at most 90 degrees.
  void rotate_from_to(vector_type v1, vector_type v2);

 private:
  friend U;
  // Construct a Basis with a rotation matrix that is all zeroes.
  Basis(std::nullptr_t) : scale_factor_(1), rotation_matrix_(rotation_matrix_type::Zero()) { }

#if CWDEBUG
 public:
  void print_on(std::ostream& os) const;
#endif
};

// End of Basis
//=============================================================================

template<ConceptUniverse U, size_t N = detail::universe_max_n_v<U>>
class StandardBasis : public Basis<U, N>
{
 public:
   Vector<N, typename Basis<U, N>::float_type> e(int i) const
   {
     return {Vector<N, typename Basis<U, N>::float_type>::eigen_type::Unit(i)};
   }
};

template<size_t N>
struct Permutation
{
  std::array<uint8_t, N> seq_{};

  Permutation(std::array<uint8_t, N> const& seq) : seq_(seq) { }

  template <std::input_iterator It>
  Permutation(It begin, It end)
  {
    ASSERT(std::distance(begin, end) == N);
    std::copy(begin, end, seq_.begin());
  }

  template<std::integral Index>
  Permutation(std::array<Index, N> const& perm) : Permutation(perm.begin(), perm.end()) { }

  template<std::integral Index>
  Permutation(Index const (&perm)[N]) : Permutation(&perm[0], &perm[N]) { }

  template<std::integral Index>
  Permutation(std::initializer_list<Index> perm) : Permutation(perm.begin(), perm.end()) { }

  template<std::integral... Index>
  requires (sizeof...(Index) == N)
  explicit Permutation(Index... indices) : seq_{static_cast<uint8_t>(indices)...} { }

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
  // Return true iff every element of seq_ falls in the range [range_start, range_end).
  bool in_range(uint8_t range_start, uint8_t range_end) const
  {
    for (uint8_t e : seq_)
      if (!(range_start <= e && e < range_end))
        return false;
    return true;
  }
#endif

#if CWDEBUG
 public:
  void print_on(std::ostream& os) const;
#endif
};

// CTAD guide: deduce N from the number of arguments.
template<std::integral... Index>
Permutation(Index...) -> Permutation<sizeof...(Index)>;

//=============================================================================
// Universe
//
// An affine space of MAX_N dimensions, that is not a subspace of any other space.
//
// Two Universe objects with different ID are completely independent.
// They have each a basis that can not be compared: there exists no transformation
// between the two basis, not even if they have the same dimension.
//
template<typename ID, size_t MAX_N, typename T = double>
struct Universe
{
  static_assert(MAX_N <= 64, "Universe::MAX_N must be <= 64 (integral masks are used).");
  using float_type = T;
  static constexpr int max_n = MAX_N;
  using basis_type = Basis<Universe>;
  using rotation_matrix_type = basis_type::rotation_matrix_type;
  using axes_mask_type = utils::uint_leastN_t<max_n>;
  using axes_type = basis_type::axes_type;
  using subset_type = basis_type::subset_type;
  using vector_type = Vector<max_n, T>;

  static StandardBasis<Universe> standard_basis;

  struct CoordinateSubspace
  {
    template<size_t N>
    using basis = Basis<Universe<ID, MAX_N, T>, N>;

    // Return a basis from a permutation.
    // If the permutation is odd, flip the sign of the basis axis at flip_output_axis_if_odd (default 0) to restore orientation.
    template<size_t N>
    static Basis<Universe<ID, MAX_N, T>, N> from_permutation(Permutation<N> permutation, T scale = 1, int flip_output_axis_if_odd = 0);
  };
};

//static
template<typename ID, size_t MAX_N, typename T>
StandardBasis<Universe<ID, MAX_N, T>> Universe<ID, MAX_N, T>::standard_basis;

//static
template<typename ID, size_t MAX_N, typename T>
template<size_t N>
Basis<Universe<ID, MAX_N, T>, N>
Universe<ID, MAX_N, T>::CoordinateSubspace::from_permutation(Permutation<N> permutation, T scale, int flip_output_axis_if_odd)
{
  DoutEntering(dc::notice, libcwd::type_info_of<Universe<ID, MAX_N, T>>().demangled_name() <<
      "::CoordinateSubspace::from_permutation(" << permutation << ", " << scale << ", " << flip_output_axis_if_odd << ")");

  Basis<Universe<ID, MAX_N, T>, N> basis{nullptr};
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
template<ConceptUniverse U, size_t N>
void Basis<U, N>::rotate_from_to(vector_type v1, vector_type v2)
{
  DoutEntering(dc::notice, "Basis<" << libcwd::type_info_of<U>().demangled_name() << ", " << N << ">::rotate_from_to(" << v1 << ", " << v2 << ")");

  constexpr float_type zero = 0;
  constexpr float_type one  = 1;

  // Normalize input vectors.
  if (!v1.normalize() || !v2.normalize())
    return;

  // Make sure the angle between v1 and v2 is at most 90 degrees by possibly flipping v2.
  float_type cos_theta = v1.dot(v2);
  if (cos_theta < zero)
  {
    v2.negate();
    cos_theta = -cos_theta;
  }

  // Can be constexpr as of C++26.
  const/*expr*/ float_type abs_relative_error = std::sqrt(std::numeric_limits<float_type>::epsilon());

  // If v1 and v2 are (almost) identical then there is no rotation.
  if (utils::almost_equal(cos_theta, one, abs_relative_error))
    return;

  // Now dot <= 1 - eps such that (see comment with utils::almost_equal):
  //
  //                           |    1 - dot    |   |     eps       |
  //      abs_relative_error = | ------------- | = | ------------- | ~ eps
  //                           | (1 + dot) / 2 |   | (2 - eps) / 2 |
  //
  // Thus dot = a·b = cos(θ) <= 1 - abs_relative_error
  // and sin(θ) >= √(1 - (1 - abs_relative_error)²) ~ √(1 - (1 - 2 abs_relative_error + ⋯)) = √(2 abs_relative_error).

  // Component of v2 orthogonal to v1.
  auto vo = v2 - cos_theta * v1;

  //                vn ^       v2=Rv1
  //                vo ^‾‾‾‾‾‾^
  //                   |     /⋮
  //     R vn          |    / ⋮
  //       ·.          |  1/  ⋮ ← sin(θ)
  //       ⋮θ ·.       |  /   ⋮
  // cos(θ)⋮     ·.    | /    ⋮
  //       ⋮        ·.θ|/θ    ⋮
  //      -+-----------o------+-------> v1
  //            ↑          ↑    1
  //         sin(θ)     cos(θ)

  float_type sin_theta = vo.norm();
  ASSERT(sin_theta > std::sqrt(abs_relative_error));

  // Normalize vo.
  auto vn = vo / sin_theta;

  // Now {v1, vn} is an orthonormal basis for the rotation plane.

  // Build the rotation matrix R that acts as:
  //   R v1 = v2 = cos_theta * v1 + sin_theta * vn
  //   R vn =     -sin_theta * v1 + cos_theta * vn
  //   R w  = w  for all w perpendicular to v1 and vn.

  auto const& e1 = v1.eigen();
  auto const& en = vn.eigen();
  auto const& eo = vo.eigen();

  rotation_matrix_type const A =
    (cos_theta - one) * (e1 * e1.transpose() + en * en.transpose()) +
                        (eo * e1.transpose() - e1 * eo.transpose());

  // Define and apply incremental rotation R.
  rotation_matrix_type const R = rotation_matrix_type::Identity() + A;
  rotation_matrix_ = R * rotation_matrix_;
}

// An n-dimensional subspace, defined by orthonormal vectors in universe coordinates.
template<ConceptUniverse U, int N>
class SubSpace
{
 public:
  static constexpr int n = N;
  using vector_type = typename U::vector_type;

 private:
  // Orthonormal vectors spanning the subspace, expressed in Universe coordinates.
  std::array<vector_type, N> orthonormal_basis_;

 public:
  SubSpace(vector_type const& normal) requires (N == 1) : orthonormal_basis_{normal}
  {
    if (!orthonormal_basis_[0].normalize())
      THROW_LALERT("Failed to normalize vector [VECTOR]", AIArgs("[VECTOR]", normal));
  }

  // Accessor for the i-th orthonormal basis vector (0 <= i < N).
  vector_type const& basis_vector(int i) const
  {
    ASSERT(i >= 0 && i < N);
    return orthonormal_basis_[i];
  }

#if CWDEBUG
 public:
  void print_on(std::ostream& os) const;
#endif
};

namespace detail {

template<ConceptUniverse U, int N>
class BasisBuilder
{
 public:
  static constexpr int n     = N;
  static constexpr int max_n = U::max_n;
  using float_type           = typename U::float_type;
  using vector_type          = typename U::vector_type;
  using rotation_matrix_type = typename U::rotation_matrix_type;
  using axes_type            = typename U::axes_type;
  using Index                = typename axes_type::Index;
  using orthogonal_subspace_type = SubSpace<U, max_n - n>;

 private:
  // The universe subspace orthogonal to the basis.
  orthogonal_subspace_type orthogonal_subspace_;

  // Orthonormal vectors spanning the processed subspace, expressed in Universe coordinates.
  std::array<vector_type, n> processed_subspace_;
  size_t processed_subspace_size_{0};            // The number of valid vectors in processed_subspace_ so far.

  // Universe standard axes that have already been “claimed” as closest to some processed_subspace_ vector.
  axes_type processed_axes_{0};

 public:
  BasisBuilder(orthogonal_subspace_type orthogonal_subspace) : orthogonal_subspace_(orthogonal_subspace) { }

  rotation_matrix_type build_rotation_matrix();

 private:
  // Project v onto the complement of the orthogonal subspace and the already selected basis vectors.
  vector_type project_onto_complement(vector_type v) const
  {
    // Remove components along the orthogonal subspace (orthonormal).
    for (int j = 0; j < orthogonal_subspace_type::n; ++j)
    {
      auto const& w = orthogonal_subspace_.basis_vector(j);
      float_type const dot = v.dot(w);
      v -= dot * w;
    }

    // Remove components along already constructed basis vectors.
    for (size_t i = 0; i < processed_subspace_size_; ++i)
    {
      auto const& u = processed_subspace_[i];
      float_type const dot = v.dot(u);
      v -= dot * u;
    }

    return v;
  }

  // Add the Universe standard axis whose projection onto the complement is longest.
  // Returns true on success; false if no suitable axis was found.
  bool add_best_axis()
  {
    using namespace utils::bitset;
    constexpr IndexPOD index_end = { max_n };

    float_type  best_norm_squared = 0;
    Index       best_axis;
    vector_type best_vector;

    for (Index axis = index_begin; axis != index_end; ++axis)
    {
      if (processed_axes_.test(axis))
        continue;

      vector_type e_axis(axis());
      vector_type projected = project_onto_complement(e_axis);
      float_type const norm_squared = projected.norm_squared();

      if (norm_squared > best_norm_squared)
      {
        best_norm_squared = norm_squared;
        best_axis         = axis;
        best_vector       = projected;
      }
    }

    // Did any remaining axis have a non-zero projection?
    if (best_norm_squared == 0)
      return false;

    float_type const len = std::sqrt(best_norm_squared);
    ASSERT(len > 0);
    best_vector /= len;

    ASSERT(processed_subspace_size_ < n);
    processed_subspace_[processed_subspace_size_] = best_vector;
    ++processed_subspace_size_;

    processed_axes_.set(best_axis);
    return true;
  }

#if CWDEBUG
 public:
  void print_on(std::ostream& os) const;
#endif
};

template<ConceptUniverse U, int N>
BasisBuilder<U, N>::rotation_matrix_type BasisBuilder<U, N>::build_rotation_matrix()
{
  // Iteratively construct N orthonormal basis vectors for the complement of
  // the orthogonal subspace by selecting the Universe standard axis whose
  // projection onto that complement is longest.
  for (int basis_row = 0; basis_row < n; ++basis_row)
  {
    if (!add_best_axis())
      THROW_LALERT("Failed to construct basis vector [INDEX].", AIArgs("[INDEX]", basis_row));
  }
  ASSERT(processed_subspace_size_ == n);

  // Fill rotation_matrix with the constructed basis vectors and the orthogonal subspace basis.
  rotation_matrix_type rotation_matrix;

  // First N rows: basis axes in Universe coordinates.
  for (int row = 0; row < n; ++row)
  {
    vector_type const& u = processed_subspace_[row];
    for (int col = 0; col < max_n; ++col)
      rotation_matrix(row, col) = u[col];
  }

  // Remaining rows: orthogonal subspace basis vectors.
  for (int row = n; row < max_n; ++row)
  {
    vector_type const& w = orthogonal_subspace_.basis_vector(row - n);
    for (int col = 0; col < max_n; ++col)
      rotation_matrix(row, col) = w[col];
  }

  return rotation_matrix;
}

} // namespace detail

// Construct a Basis from a BasisBuilder.
//
// The BasisBuilder describes the subspace orthogonal to the Basis (its normal
// subspace) and keeps track of Universe axes that were already used to define
// basis directions. Here we construct an orthonormal basis of the complement
// of that orthogonal subspace by selecting Universe standard axes whose
// projections onto that complement are as long as possible.
//
// The resulting rotation_matrix_ converts Universe coordinates to coordinates
// in this Basis: the first N rows are the N basis axes in Universe
// coordinates, and the remaining rows are an orthonormal completion given by
// the orthogonal subspace.
template<ConceptUniverse U, size_t N>
Basis<U, N>::Basis(detail::BasisBuilder<U, N> builder) : scale_factor_{1}, rotation_matrix_{builder.build_rotation_matrix()} { }

#if CWDEBUG
template<ConceptUniverse U, size_t N>
void Basis<U, N>::print_on(std::ostream& os) const
{
  os << "{scale_factor:" << scale_factor_ << ", rotation_matrix:\n" << rotation_matrix_ << '}';
}

template<size_t N>
void Permutation<N>::print_on(std::ostream& os) const
{
  os << '(';
  char const* separator = "";
  for (int element : seq_)
  {
    os << separator << element;
    separator = " ";
  }
  os << ')';
}

template<ConceptUniverse U, int N>
void SubSpace<U, N>::print_on(std::ostream& os) const
{
  os << "{orthonormal_basis:" << utils::print_range(orthonormal_basis_.begin(), orthonormal_basis_.end()) << "}";
}

template<ConceptUniverse U, int N>
void detail::BasisBuilder<U, N>::print_on(std::ostream& os) const
{
  std::string const processed_axes = processed_axes_.to_string();
  os << "{orthogonal_subspace_:" << orthogonal_subspace_ <<
    ", processed_subspace:" << utils::print_range(processed_subspace_.begin(), processed_subspace_.begin() + processed_subspace_size_) <<
    ", processed_subspace_size:" << processed_subspace_size_ << ", processed_axes:" << processed_axes.substr(processed_axes.length() - N) << "}";
}
#endif

} // namespace math
