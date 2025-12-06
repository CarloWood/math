#pragma once

#include "math/Vector.h"
#include "debug.h"

namespace math {

enum Sign
{
  negative = -1,      // The signed distance is negative: the point is on the -N side of the plane.
  in_plane = 0,       // The signed distance is very small: the point can be considered to be in the plane.
  positive = 1,       // The signed distance is positive: the point is on the +N side of the plane.
};

// Hyperplane
//
// An (n-1)-dimensional hyperplane orthogonal to a given normal vector N,
// located at the signed offset -b/‖N‖ from the origin.
//
// This class follows the convention originally set by Ludwig Otto Hesse, with respect to
// the signed distance. The Hesse form (https://en.wikipedia.org/wiki/Hesse_normal_form) of
// a hyperplane is where N is a unit normal vector (‖N‖ is one), indicated by giving N a hat,
// and reads:
//     ∧
//     N·X = ρ
//
// Here ρ is the signed offset of the plane from the origin: ρ is positive iff
// the plane lies in the direction that the normal points:
//
//              hyperplane
//                 |
//          N      |  N
//        O---->   |---->
//     origin      |
//                 P----------->A
//        -------->|  S(A) > 0
//           ρ > 0 |
//
// However, the signed distance S(A) from the plane to a point A is relative to the plane
// itself and positive if A is on the side of the plane that N is pointing towards
// (e.g. https://courses.csail.mit.edu/6.036/spring_2016/assignments/hw0_final.pdf (2d)).
//
// Let a point A be given by its projection P onto the plane and its signed distance s from
// the plane:
//               ∧
//     A = P + s N
//                                                 ∧
// Then the signed distance of A to the hyperplane N·X - ρ = 0 is:
//            ∧         ∧        ∧         ∧            ∧ ∧
//     S(A) = N·A - ρ = N·(P + s N) - ρ = (N·P - ρ) + s N·N = s
//
// because P is in the plane.
//
// Note that the signed distance of the origin itself is:
//            ∧
//     S(O) = N·O - ρ = -ρ,
//
// and has the opposite sign of the plane offset!
//
//              hyperplane
//                 |
//                 |  N
//     origin      |---->
//        O<-------Q-.         ∧
//            /    |  \_ Q = ρ·N
//   S(O) = -ρ < 0 |
//
// The constructor of this class accepts the normal vector N (not necessarily
// a unit vector) as first argument and the negative of the dot product of N with
// some point P in the plane (b = -N·P).
//
// Note that in the literature it is common to write the equation of a hyperplane as:
//
//     N·X + b = 0,
//
// Hyperplane also stores N and b.
//
// Because P is in the plane, we have N·P + b = 0 --> b = -N·P.
// If we write the plane equation in Hesse form, by dividing by ‖N‖,
//
//     (N·X + b)/‖N‖ = 0 -> (N/‖N‖)·X = -b/‖N‖
//
// we see that b = -ρ ‖N‖, the plane offset multiplied by -‖N‖.
//
// Note:
// * The signed offset of the plane from the origin is ρ = -b/‖N‖.
//              ∧
// * The point ρN = (-b/‖N‖)(N/‖N‖) = (-b/‖N‖²)N is the projection of the origin onto the hyperplane.
//
// * The perpendicular distance between the origin and the plane is |ρ| = |b|/‖N‖.
//
// * The signed distance of the origin from the plane is -ρ = b/‖N‖.
//
template<int N, std::floating_point T = double>          // N is the number of dimensions here.
class Hyperplane
{
 public:
  using scalar_type = T;
  using vector_type = math::Vector<N, T>;                // N is the number of dimensions here; below 'N' in the comments refers to the normal_.

 private:
  vector_type normal_;                                   // The normal of the hyperplane.
  T b_;                                                  // The plane constant where N·X + b = 0.

 public:
  // Create a hyperplane that satisfies N·X + b = 0, where b = -N·P for some P on the plane.
  Hyperplane(vector_type const& normal, typename vector_type::scalar_type b) : normal_(normal), b_(b) { }

  // Accessors.
  vector_type const& normal() const { return normal_; }
  T b() const { return b_; }

  // Return the signed distance in Euclidean units (positive in +N direction).
  T signed_distance(vector_type const& A) const
  {
    //        ∧
    // S(A) = N·A - ρ = (N·A + b) / ‖N‖.
    return (normal_.dot(A) + b_) / normal_.norm();
  }

  // Return h such that A_projected = A - h N lies on the plane.
  T height_along_N(vector_type const& A) const
  {
    // 0 = N·A_projected + b = N·(A - h N) + b = N·A + b - h N·N -->
    // h = (N·A + b) / N·N.
    return (normal_.dot(A) + b_) / normal_.norm_squared();
  }

  // Orthogonal projection of A onto the plane.
  vector_type project(vector_type const& A) const
  {
    return A - height_along_N(A) * normal_;
  }

  // Return which side of the plane A is on.
  Sign side(vector_type const& A) const
  {
    T h = height_along_N(A);
    T ah = std::abs(h);
    constexpr T abs_relative_error = 1e-6;
    if (ah < abs_relative_error) [[unlikely]]
      return in_plane;
    return h > 0.0 ? positive : negative;
  }

  // Return intersection of the line through C1 and C2 with this Hyperplane.
  // Only call this function for points C1 and C2 where their difference is not perpendicular to N,
  // for example, where `distance` returns significantly different values for both points.
  vector_type intersection(vector_type const& C1, vector_type const& C2) const
  {
    // Let E be a line through C1 and C2: E: C1 + g(C2 - C1), where g parameterizes the points on E.
    // Fill that into the line equation to find the intersection:
    // N·(C1 + g(C2 - C1)) + b = 0 --> N·C1 + b + g N·(C2 - C1) = 0 --> g = -(N·C1 + b) / N·(C2 - C1)
    vector_type diff = C2 - C1;
    T g = -(normal_.dot(C1) + b_) / (normal_.dot(diff));
    return C1 + g * diff;
  }

#ifdef CWDEBUG
  void print_on(std::ostream& os) const
  {
    os << "{normal:" << normal_ << ", b:" << b_ << "}";
  }
#endif
};

} // namespace math
