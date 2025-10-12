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

  for (CornerIndex ci = C_.ibegin(); ci != C_.iend(); ++ci)
  {
    for (int d = 0; d < n; ++d)
    {
      int bit = 1 << d;
      CornerIndex ci2(ci.get_value() | bit);
      if (side[ci] != side[ci2])
      {
        // Found two corners on opposite sides of the hyperplane.
        intersections.push_back(plane.intersection(C_[ci], C_[ci2]));
      }
    }
  }

  return intersections;
}

} // namespace math
