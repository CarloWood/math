#pragma once

#include "math/Hyperplane.h"
#include "utils/Vector.h"
#include <concepts>
#include <vector>
#include "debug.h"

namespace math {

namespace category {
struct Hyperblock;
} // namespace category

using CornerIndex = utils::VectorIndex<category::Hyperblock>;

template<int n , std::floating_point T = double>
class Hyperblock
{
 public:
  using VectorType = math::Vector<n, T>;
  static constexpr int number_of_corners = 1 << n;

 private:
  utils::Vector<VectorType, CornerIndex> C_;          // The 2^n corners of the hyperblock.

 public:
  // Construct an axis-aligned hyperblock from two opposite corner vectors.
  Hyperblock(VectorType const& c1, VectorType const& c2) : C_(number_of_corners)
  {
    VectorType base = c2 - c1;
    for (CornerIndex ci = C_.ibegin(); ci != C_.iend(); ++ci)
    {
      VectorType c;
      c.eigen().setZero();
      for (int d = 0; d < n; ++d)
      {
        int bit = 1 << d;
        if ((ci.get_value() & bit))
          c[d] += base[d];
      }
      C_[ci] = c1 + c;
    }
  }

  // Return the coordinates of corner ci as a vector.
  VectorType const& operator[](CornerIndex ci) const
  {
    return C_[ci];
  }

  CornerIndex ibegin() const { return C_.ibegin(); }
  CornerIndex iend() const { return C_.iend(); }

  std::vector<VectorType> intersection_points(Hyperplane<n, T> const& plane);

#ifdef CWDEBUG
  void print_on(std::ostream& os) const
  {
    LIBCWD_USING_OSTREAM_PRELUDE
    os << "{C:" << C_ << "}";
  }
#endif
};

template<int n, std::floating_point T>
std::vector<typename Hyperblock<n, T>::VectorType> Hyperblock<n, T>::intersection_points(Hyperplane<n, T> const& plane)
{
  DoutEntering(dc::notice, "Hyperblock<" << n << ", " << type_info_of<T>().demangled_name() << ">::intersection_points(" << plane << ")");
  std::vector<VectorType> intersections;

  utils::Vector<Sign, CornerIndex> side(number_of_corners);
  for (CornerIndex ci = C_.ibegin(); ci != C_.iend(); ++ci)
    side[ci] = plane.side(C_[ci]);

  Dout(dc::notice|continued_cf, "side = ");
  for (CornerIndex ci = C_.ibegin(); ci != C_.iend(); ++ci)
    Dout(dc::continued, (side[ci] == positive ? "+" : side[ci] == in_plane ? "0" : "-"));
  Dout(dc::finish, "");

  for (int e = 0; e < n; ++e)         // Consider all edges along axis `e`.
  {
    // Run over all corners that have coordinate 0 for axis `e`.
    for (CornerIndex ci_e0 = C_.ibegin(); ci_e0 != C_.iend(); ++ci_e0)
    {
      size_t const mask_e = size_t{1} << e;
      size_t const bit_e = ci_e0.get_value() & mask_e;
      if (bit_e != 0)
        continue;

      // The corner on the opposite side along axis `e`.
      CornerIndex ci_e1(ci_e0.get_value() | mask_e);

      // The these two corners are on opposite sides of the hyper plane then the edge between them intersects.
      if (side[ci_e0] != side[ci_e1])
      {
        // Run over all 2-faces that this edge is a part of.
        for (int f = 0; f < n; ++f)
        {
          if (f == e)
            continue;
          // Dimensions e and f span a 2-face.
          size_t const mask_f = size_t{1} << f;
          size_t const bit_f = ci_e0.get_value() & mask_f;
          // Determine the clock-wise order in which the four corners of this 2-face occur, looking at it from the outside of the hyper-block.
          //
          // If e < f, then this picture matches the discussion with chatgpt:
          // See https://chatgpt.com/share/692e1f87-b2b0-800d-960e-dea2e1b0b83a
          //
          //                             +1            -1       <-- change in e.
          //                            f=1           f=0
          //                         .---^---.     .---^---.
          // Order A means:   00 --> 01 --> 11 --> 10 --> 00
          //                         ||            ||
          //                      +1 ||         -1 ||           <-- change in e.
          //                     f=0 ||        f=1 ||
          //                  .---^---.     .---^---.
          // Order B means:   00 --> 10 --> 11 --> 01 --> 00
          //                         ^^            ^^
          //                      e_/  \_f      e_/  \_f
          //
          // If f < e, then the picture holds:
          //                      +1            -1              <-- change in e.
          //                     f=0           f=1
          //                  .---^---.     .---^---.
          // Order A means:   00 --> 01 --> 11 --> 10 --> 00
          //                         ||            ||
          //                         ||  +1        ||  -1       <-- change in e.
          //                         || f=1        || f=0
          //                         .---^---.     .---^---.
          // Order B means:   00 --> 10 --> 11 --> 01 --> 00
          //                         ^^            ^^
          //                      f_/  \_e      f_/  \_e        <-- iff f < e.
          //
          bool order_A = (std::popcount(ci_e0.get_value() & ~mask_f) + n + (f - e)) & 1;
          // All other coordinates of this vector are zero.
          int const edge_direction = /* vector along the edge with coordinate for axis `e`: */ (order_A == (bit_f == 0)) != (f < e) ? -1 : 1;
          //
          // Get the (sign of the) component of the normal of the hyper plane along the edge direction.
          int const n_e = plane.normal()[e] < 0.0 ? -1 : 1;
          // Now we can calculate the (sign of the) dot product between the hyperplane normal and the directional vector that goes clockwise around the 2-face.
          // If this sign is positive then this is the entry point of this 2-face.
          bool is_entry_point = edge_direction * n_e > 0;

          // Example:                             end (defining the direction of the line segment SE)
          //            (0,1)                     /               (1,1)
          //               O---------------------E-------->--------O
          //               |                    /                  |
          //               |                   /                   |   clock wise direction
          //               ^         N .      /                    | /
          //               |            '-_  /                     v
          //               |               'O                      |
          //               |edge_direction /                       |
          //               |      \       /                        |
          //               O------<------S-------------------------O
          //      ci_e0=(0,0)           /                   ci_e1=(1,0)
          //                          start                          ^
          //                                                          \ bit_f
          //
          // The point F has corners ci_e0 (0,0) and ci_e1 (1,0), e=0 (the first coordinate).
          // The 2-face is spanned by e=0 and f=1 (the second (and only other) coordinate).
          // bit_f = 0
          // ci_e0.get_value() = b00, n=2, f-e=1 --> order_A = true: going clockwise we encounter in order: (0,0) --> (0,1) --> (1,1) --> (1,0) --> (0,0)
          // Thus edge_direction=-1 (going from x=1 to x=0).
          // The component of the normal along e also -1 (N points "to the left").
          // Therefore is_entry_point=true : F is the point where the line enters the rectangle.
        }

        // Found two corners on opposite sides of the hyperplane.
        intersections.push_back(plane.intersection(C_[ci_e0], C_[ci_e1]));
      }
    }
  }

  return intersections;
}

} // namespace math
