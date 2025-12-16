#pragma once

#include "utils/BitSet.h"
#include "utils/create_mask.h"
#include "utils/deposit_extract.h"
#include "utils/uint_leastN_t.h"
#include "utils/VectorIndex.h"
#include "math/binomial.h"
#include "debug.h"
#ifdef CWDEBUG
#include "utils/has_print_on.h"
#endif

namespace math {
#ifdef CWDEBUG
using utils::has_print_on::operator<<;
#endif

template<int n, int k>
struct kFaceData;

// n : number of dimensions of the hypercube.
// k : number of dimensions of the encoded k-face.
//
// This encodes a k-face as a bitmask that can be used as an index into a Vector.
// The bitmask exists of a 'rank' for the k axes that the k-face aligns with,
// followed by n-k 'fixed bits' for the other axes:
//
//     <rrr><fffffff>
//      ^^^  ^^^^^^^
//
// The fixed bits (<fffffff>) are those bits from the coordinate mask of each
// corner of the k-face that do not change: the bits that belong the one of the
// k axis (that are thus not fixed) are removed.
//
// The rank bits (<rrr>) encode the lexiographic rank of the k axes that the k-face
// aligns with, as returned by the member function `rank`.
//
template<int n, int k>
class kFaceIndex;

// Note that a 0-face is a corner (k=0), a 1-face is an edge (k=1).
template<int n>
using CornerIndex = kFaceIndex<n, 0>;

template<int n>
using EdgeIndex = kFaceIndex<n, 1>;

template<int n, int k>
struct kFaceData
{
  using axes_type = utils::BitSet<utils::uint_leastN_t<n>>;

  axes_type k_axes;             // The k axes that are aligned with this k-face.
  CornerIndex<n> zero_corner;   // The corner index, that is part of the k-face, with zeroes for each of the k axes.
};

template<int n, int k>
class kFace;

template<int n, int k>
class kFaceIndex : public utils::VectorIndex<kFaceData<n, k>>
{
 public:
  // A bitset that can hold at least n bits, where each bit represents an axis.
  using axes_type = utils::BitSet<utils::uint_leastN_t<n>>;
  static constexpr uint32_t fixed_mask = utils::create_mask<uint32_t, n - k>();
  static constexpr uint32_t corner_mask = utils::create_mask<uint32_t, n>();

  // The number of k-faces of an n-cube.
  static constexpr size_t size = static_cast<size_t>(binomial(n, k)) << (n - k);
  // The number of (k-1)-faces of a given k-face.
  static constexpr size_t number_of_facets = 2 * k;

  static uint32_t rank(axes_type k_axes);
  static axes_type unrank(uint32_t r);

  // Construct an undefined kface.
  constexpr kFaceIndex() = default;

  // Construct a corner.
  // Alow constructing a k-face index from ibegin and iend (for use with for loops).
  constexpr kFaceIndex(size_t value) : utils::VectorIndex<kFaceData<n, k>>(value)
  {
#if CW_DEBUG
    if constexpr (k == 0)
      // All unused bits must be zero, except if value is possibly `end` by being one larger than the maximum allowed index.
      ASSERT((value & corner_mask) == value || value == size_t{1} << n);
    else
      // This should only be used as part of calling ibegin() and/or iend().
      ASSERT(value == size_t{0} || static_cast<size_t>(value) == size);
#endif
  }

  // Construct a kFaceIndex from some given kFace.
  kFaceIndex(kFaceData<n, k> const& kface) : utils::VectorIndex<kFaceData<n, k>>(static_cast<size_t>(rank(kface.k_axes) << (n - k) | kface.zero_corner.remove_bits(kface.k_axes))) { }

  // Construct a zero-corner by inserting zeroes into `fixed_bits` at the positions of `k_axes`.
  kFaceIndex(axes_type k_axes, uint32_t fixed_bits) requires (k == 0) : utils::VectorIndex<kFaceData<n, k>>(fixed_bits) { insert_bits(k_axes); }

  std::array<kFaceIndex<n, k - 1>, number_of_facets> facet_indexes() const requires (k > 0);

  kFace<n, k> as_kface() const;

 private:
  template<int, int> friend class kFaceIndex;

  uint32_t remove_bits(axes_type k_axes) const requires (k == 0)
  {
    ASSERT(!this->undefined());
    static_assert(n <= 32, "remove_bits returns at most 32 bits.");

    using mask_type = typename axes_type::mask_type;

    mask_type const value = static_cast<mask_type>(this->get_value());
    mask_type const mask = ~k_axes();

    return static_cast<uint32_t>(utils::extract_bits(value, mask));
  }

  void insert_bits(axes_type k_axes)
  {
    ASSERT(!this->undefined());
    static_assert(n <= 32, "insert_bits returns at most 32 bits.");

    using mask_type = typename axes_type::mask_type;

    mask_type const value = static_cast<mask_type>(this->get_value());
    mask_type const mask = ~k_axes();

    this->m_value = utils::deposit_bits(value, mask);
  }

 public:
  static constexpr kFaceIndex ibegin() { return kFaceIndex{size_t{0}}; }
  static constexpr kFaceIndex iend() { return kFaceIndex{size}; }
};

template<int n, int k>
class kFace : public kFaceData<n, k>
{
 public:
  using axes_type = kFaceData<n, k>::axes_type;

  kFace(kFaceData<n, k> const& data) : kFaceData<n, k>(data) { }

 private:
  friend class kFaceIndex<n, k>;
  kFace<n, k - 1> get_facet(axes_type axis, int c) const
  {
    // `axis` must be an axis from k_axes.
    ASSERT(axis.is_single_bit() && !(this->k_axes & axis).none());

    size_t zero_corner = this->zero_corner.get_value();
    if (c)
      zero_corner |= axis();
    else
      zero_corner &= ~axis();

    return {{this->k_axes & ~axis, zero_corner}};
  }

#ifdef CWDEBUG
 public:
  void print_on(std::ostream& os) const
  {
    os << "{k_axes:" << this->k_axes << ", zero_corner:" << this->zero_corner << "}";
  }
#endif
};

template<int n, int k>
kFace<n, k> kFaceIndex<n, k>::as_kface() const
{
  uint32_t const mask = this->get_value();
  axes_type k_axes = unrank(mask >> (n - k));
  CornerIndex<n> zero_corner(k_axes, mask & fixed_mask);
  return {{k_axes, zero_corner}};
}

template<int n, int k>
std::array<kFaceIndex<n, k - 1>, kFaceIndex<n, k>::number_of_facets>
kFaceIndex<n, k>::facet_indexes() const requires (k > 0)
{
  std::array<kFaceIndex<n, k - 1>, kFaceIndex<n, k>::number_of_facets> result;

  ASSERT(!this->undefined());

  // Fill the result array with all facets.
  kFace<n, k> kface = as_kface();
  int i = 0;
  for (auto axis : kface.k_axes)
    for (int c = 0; c <= 1; ++c)
      result[i++] = kface.get_facet(axis, c);

  return result;
}

// Returns the rank for this choice of axes for the set of all possible ways one can choose k axes out of n.
//static
template<int n, int k>
uint32_t kFaceIndex<n, k>::rank(axes_type k_axes)
{
  // There should be one bit set for each of the k axes.
  ASSERT(k_axes.count() == k);

  using namespace utils::bitset;

  // Calculate the rank according to:
  //
  //                         k  dᵢ-1 ⎛n-1-x⎞
  //   rank(d₁, d₂, …, dₖ) = ∑   ∑   ⎜     ⎟
  //                        i=1 x=mᵢ ⎝ k-i ⎠
  //
  // where m₁=0, mᵢ = dᵢ₋₁ + 1 for i>1.
  uint32_t sum = 0;

  int mi = 0;
  Index d_i = index_pre_begin;
  for (int i = 1; i <= k; ++i)
  {
    // Advance d_i to the first/next axis index.
    d_i.next_bit_in(k_axes);
    for (int x = mi; x <= d_i() - 1; ++x)
      sum += binomial(n - 1 - x, k - i);
    mi = d_i() + 1;
  }
  return sum;
}

// Given a rank in [0, C(n,k)) return the corresponding axes BitSet.
//static
template<int n, int k>
typename kFaceIndex<n, k>::axes_type kFaceIndex<n, k>::unrank(uint32_t r)
{
  using namespace utils::bitset;
  // Index of the least significant bit.
  constexpr IndexPOD ilsb = { 0 };
  // Initialize the result variable.
  axes_type axes{axes_type::zero};

  // Our binomials are never negative (n and k are positive).
  ASSERT(r < (uint32_t)binomial(n, k));

  Index mi = ilsb;
  for (int i = 1; i <= k; ++i)
  {
    Index d = mi;
    // Find smallest d such that the block size falls past r.
    for (;; ++d)
    {
      uint32_t const block = binomial(n - 1 - d(), k - i);
      if (r < block)
        break;
      r -= block;
    }
    axes.set(d);
    mi = d + 1;
  }

  return axes;
}

} // namespace math
