#pragma once

#include "math/Hyperplane.h"
#include "math/kFace.h"
#include "utils/Vector.h"
#include <concepts>
#include <vector>
#include <queue>
#include <set>
#ifdef CWDEBUG
#include "utils/has_print_on.h"
#endif
#include "debug.h"

namespace math {

namespace detail {
#ifdef CWDEBUG
// Classes in this namespace define a print_on method.
using utils::has_print_on::operator<<;
#endif

inline size_t to_mask(int dim)
{
  return size_t{1} << dim;
}

template<int n , std::floating_point T>
class EdgeId
{
 private:
  CornerIndex<n> lo_;   // The lower-valued corner index of the edge (lexicographic first endpoint).
  CornerIndex<n> hi_;   // The higher-valued corner index of the edge (lexicographic second endpoint).

 public:
  EdgeId() = default;   // Construct an EdgeId with undefined corner indexes.

  EdgeId(CornerIndex<n> ci0, CornerIndex<n> ci1) :
    lo_(std::min(ci0, ci1)),
    hi_(std::max(ci0, ci1))
  {
    ASSERT(!lo_.undefined() && !hi_.undefined());
  }

  // Order EdgeId's by lo_ first.
  friend bool operator<(EdgeId const& lhs, EdgeId const& rhs)
  {
    // Don't call operator< on undefined EdgeId's.
    ASSERT(!lhs.undefined() && !rhs.undefined());
    if (lhs.lo_ != rhs.lo_)
      return lhs.lo_ < rhs.lo_;
    return lhs.hi_ < rhs.hi_;
  }

  bool undefined() const
  {
    return lo_.undefined();
  }

#ifdef CWDEBUG
  void print_on(std::ostream& os) const
  {
    os << "{lo:" << lo_ << ", hi:" << hi_ << "}";
  }
#endif
};

template<int n , std::floating_point T>
class FaceId
{
 private:
  CornerIndex<n> base_; // A corner of the face with the varying bits (dim0, dim1) cleared; encodes the fixed coordinates.
  int dim_lo_;          // Smaller varying dimension of the 2-face.
  int dim_hi_;          // Larger varying dimension of the 2-face.

 public:
  FaceId() = default;   // Construct an uninitialized FaceId.

  FaceId(CornerIndex<n> a_corner, int dim0, int dim1) :
    base_(a_corner.get_value() & ~(to_mask(dim0) | to_mask(dim1))),
    dim_lo_(std::min(dim0, dim1)),
    dim_hi_(std::max(dim0, dim1))
  {
  }

  friend bool operator<(FaceId const& lhs, FaceId const& rhs)
  {
    if (lhs.base_ != rhs.base_)
      return lhs.base_ < rhs.base_;
    if (lhs.dim_lo_ != rhs.dim_lo_)
      return lhs.dim_lo_ < rhs.dim_lo_;
    return lhs.dim_hi_ < rhs.dim_hi_;
  }

#ifdef CWDEBUG
  void print_on(std::ostream& os) const
  {
    os << "{base:" << base_ << ", dim_lo:" << dim_lo_ << ", dim_hi:" << dim_hi_ << "}";
  }
#endif
};

template<int n , std::floating_point T>
class FaceSegment
{
 private:
  EdgeId<n, T> entry_point_;
  EdgeId<n, T> exit_point_;

 public:
  FaceSegment() = default;

  void add(EdgeId<n, T> edge_id, bool entry_point)
  {
    if (entry_point)
      entry_point_ = edge_id;
    else
      exit_point_ = edge_id;
  }

  // Accessor.
  EdgeId<n, T> entry_point() const
  {
    // Call add for the entry point.
    ASSERT(!entry_point_.undefined());
    return entry_point_;
  }

  EdgeId<n, T> exit_point() const
  {
    // Call add for the exit point.
    ASSERT(!exit_point_.undefined());
    return exit_point_;
  }
};

} // namespace detail

template<int n , std::floating_point T = double>
class Hyperblock
{
 public:
  struct IntersectionPointCategory;
  using IntersectionPointIndex = utils::VectorIndex<IntersectionPointCategory>;

 public:
  static constexpr int number_of_corners = 1 << n;
  using vector_type = math::Vector<n, T>;
  using IntersectionPoints = utils::Vector<vector_type, IntersectionPointIndex>;

 private:
  utils::Vector<vector_type, CornerIndex<n>> C_;        // The 2^n corners of the hyperblock.

 public:
  // Construct an axis-aligned hyperblock from two opposite corner vectors.
  Hyperblock(vector_type const& c1, vector_type const& c2) : C_(number_of_corners)
  {
    vector_type base = c2 - c1;
    for (CornerIndex<n> ci = C_.ibegin(); ci != C_.iend(); ++ci)
    {
      vector_type c;
      c.eigen().setZero();
      for (int d = 0; d < n; ++d)
      {
        size_t const bit = detail::to_mask(d);
        if ((ci.get_value() & bit))
          c[d] += base[d];
      }
      C_[ci] = c1 + c;
    }
  }

  // Return the coordinates of corner ci as a vector.
  vector_type const& operator[](CornerIndex<n> ci) const
  {
    return C_[ci];
  }

  // Run over all corners of the HyperBlock by index.
  CornerIndex<n> ibegin() const { return C_.ibegin(); }
  CornerIndex<n> iend() const { return C_.iend(); }

  IntersectionPoints intersection_points(Hyperplane<n, T> const& plane) const;

#ifdef CWDEBUG
  void print_on(std::ostream& os) const
  {
    LIBCWD_USING_OSTREAM_PRELUDE
    os << "{C:" << C_ << "}";
  }
#endif

 private:
  struct Edge
  {
    vector_type intersection_point_;            // The intersection point of the HyperPlane with this Edge.
    std::set<detail::FaceId<n, T>> entry_faces; // List of all faces that (the intersection point on) this Edge is an entry point for.
    std::set<detail::FaceId<n, T>> exit_faces;  // List of all faces that (the intersection point on) this Edge is an exit point for.
    mutable bool visited_{false};

    void add_intersection_point(vector_type const& intersection_point);
    void store(bool entry_point, detail::FaceId<n, T> const& face_id);

    // Accessors.

    // Return true if this edge is the entry point for at least one 2-face.
    bool is_entry_point() const { return !entry_faces.empty(); }
    bool is_visited() const { return visited_; }
    auto entry_face_begin() const { return entry_faces.begin(); }
    auto entry_face_end() const { return entry_faces.end(); }

    // Return the lexiographically smallest face id that this edge is an entry point of.
    detail::FaceId<n, T> first_entry_face() const
    {
      // Only call this function if is_entry_point returns true.
      ASSERT(is_entry_point());
      return *entry_faces.begin();
    }

    void visited() const
    {
      // Mark this edge as visited.
      visited_ = true;
    }
  };

  void breadth_first_search(
    Edge const& current_edge,
    std::map<detail::EdgeId<n, T>, Edge> const& edges,
    std::map<detail::FaceId<n, T>, detail::FaceSegment<n, T>> const& face_segments,
    IntersectionPoints& intersection_points) const;
};

template<int n , std::floating_point T>
void Hyperblock<n, T>::Edge::add_intersection_point(vector_type const& intersection_point)
{
  intersection_point_ = intersection_point;
}

template<int n , std::floating_point T>
void Hyperblock<n, T>::Edge::store(bool entry_point, detail::FaceId<n, T> const& face_id)
{
  if (entry_point)
    entry_faces.insert(face_id);
  else
    exit_faces.insert(face_id);
}

template<int n, std::floating_point T>
typename Hyperblock<n, T>::IntersectionPoints Hyperblock<n, T>::intersection_points(Hyperplane<n, T> const& plane) const
{
  IntersectionPoints intersection_points;

  utils::Vector<Sign, CornerIndex<n>> side(number_of_corners);
  for (CornerIndex ci = C_.ibegin(); ci != C_.iend(); ++ci)
    side[ci] = plane.side(C_[ci]);

  using namespace detail;

  // Push in_plane corners to one side, so we consistently count 'cut' edges such that
  // the respective corner is counted the minimum number of times (zero or one time in the 3D case):
  // push it to the side that has the most adjacent corners.
  //
  // For example,
  // consider the case where two adjacent corners, A and B, are in_plane (0).
  // Let the corner C be adjacent to A and on the 'negative' side of the plane,
  // and the corner D also adjacent to A but on the 'positive' side of the plane.
  // Likewise for E and F:
  //
  //               D               F
  //                +               +        n
  //                 \               \       ^
  //                  \      edge     \      │
  //          ⋯⋯⋯⋯⋯⋯⋯A⋯0───────────────0⋯B⋯⋯⋯⋯⋯⋯ (the plane)
  //                  /               /
  //                 /               /
  //                -               -
  //               C               E
  //
  // Then either A or B is encountered first. Lets say A.
  // In that case we count an equal number of adjacent corners to be on the
  // positive side as on the negative side (one each), in which the tie-breaker
  // is used to put A on the positive side (the plane is assumed to cut AC, not AD).
  //
  // Next B is encountered which now sees one adjacent corner (E) to be
  // on the negative side and two (A and F) to be on the positive side.
  // Therefore B is also put on the positive side.
  //
  // This is consistent with moving the whole plane an infinitesimal bit
  // in the direction opposite to its normal:
  //
  //               D               F
  //                +               +
  //                 \               \       n
  //                  \      edge     \      ^
  //                 A +───────────────+ B   │
  //          ⋯⋯⋯⋯⋯⋯⋯⋯/⋯⋯⋯⋯⋯⋯⋯⋯⋯⋯⋯⋯⋯⋯⋯/⋯⋯⋯⋯⋯⋯⋯⋯⋯ (the plane, assumed to cut AC and BE).
  //                 /               /
  //                -               -
  //               C               E
  //
  for (CornerIndex ci = C_.ibegin(); ci != C_.iend(); ++ci)
  {
    if (side[ci] != in_plane)
      continue;

    int pos = 0, neg = 0;
    for (int d = 0; d < n; ++d)
    {
      size_t const bit = to_mask(d);
      CornerIndex<n> const adjacent_ci{ci.get_value() ^ bit};
      Sign const adjacent_side = side[adjacent_ci];
      if (adjacent_side == positive)
        ++pos;
      else if (adjacent_side == negative)
        ++neg;
    }

    // Any side (above or below the plane) will do, but this is a good choice.
    side[ci] = (pos >= neg) ? positive : negative;
  }

  // A map from EdgeId to Edge data.
  std::map<EdgeId<n, T>, Edge> edges;
  // A map from FaceId to the two intersected edges of that 2-face.
  std::map<FaceId<n, T>, FaceSegment<n, T>> face_segments;

  for (int e = 0; e < n; ++e)         // Consider all edges along axis `e`.
  {
    // Run over all corners that have coordinate 0 for axis `e`.
    for (CornerIndex ci_e0 = C_.ibegin(); ci_e0 != C_.iend(); ++ci_e0)
    {
      size_t const mask_e = to_mask(e);
      size_t const bit_e = ci_e0.get_value() & mask_e;

      // Skip the corners where coordinate `e` isn't zero.
      if (bit_e != 0)
        continue;

      // The corner on the opposite side along axis `e`.
      CornerIndex<n> ci_e1(ci_e0.get_value() | mask_e);

      // Skip edges that don't intersect with the hyper plane.
      if (side[ci_e0] == side[ci_e1])
        continue;

      // Found two corners on opposite sides of the hyperplane.

      // Construct a unique ID for the current edge between the two corners.
      EdgeId<n, T> const edge_id(ci_e0, ci_e1);
      Edge& current_edge = edges[edge_id];
      // Calculate the intersection point on this edge and store it in edges.
      current_edge.add_intersection_point(plane.intersection(C_[ci_e0], C_[ci_e1]));

      // Run over all 2-faces that this edge is a part of.
      for (int f = 0; f < n; ++f)
      {
        if (f == e)
          continue;

        // Dimensions e and f span a 2-face.
        size_t const mask_f = to_mask(f);
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
        // If this sign is positive then this is the entry point (S) of this 2-face, otherwise it is the exit point (E).
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
        // The point S has corners ci_e0 (0,0) and ci_e1 (1,0), e=0 (the first coordinate).
        // The 2-face is spanned by e=0 and f=1 (the second (and only other) coordinate).
        // bit_f = 0
        // ci_e0.get_value() = b00, n=2, f-e=1 --> order_A = true: going clockwise we encounter in order: (0,0) --> (0,1) --> (1,1) --> (1,0) --> (0,0)
        // Thus edge_direction=-1 (going from x=1 to x=0).
        // The component of the normal along e also -1 (N points "to the left").
        // Therefore is_entry_point=true : S is the point where the line enters the rectangle.

        // Construct a unique ID for the current 2-face.
        FaceId<n, T> const face_id(ci_e0, e, f);
        // Remember if the intersection point on the current edge is the entry point for the current face or not.
        current_edge.store(is_entry_point, face_id);

        // Remember per face which edges are cut and if that is an entry or exit point.
        face_segments[face_id].add(edge_id, is_entry_point);
      }
    }
  }

  // Does the hyper plane intersect with any edge at all?
  if (!edges.empty())
  {
    // Begin with the lexiographically smallest EdgeId that is an entry point for at least one face.
    auto first_edge_iter = edges.begin();
    while (!first_edge_iter->second.is_entry_point())
    {
      ++first_edge_iter;
      ASSERT(first_edge_iter != edges.end());
    }
    Edge const& first_edge{first_edge_iter->second};

    // Run over all edges from here, following face segments for which we are the entry point in a breath-first manner and store the intersection point of each edge.
    breadth_first_search(first_edge, edges, face_segments, intersection_points);
  }

  return intersection_points;
}

template<int n, std::floating_point T>
void Hyperblock<n, T>::breadth_first_search(
    Edge const& current_edge,
    std::map<detail::EdgeId<n, T>, Edge> const& edges,
    std::map<detail::FaceId<n, T>, detail::FaceSegment<n, T>> const& face_segments,
    IntersectionPoints& intersection_points) const
{
  using namespace detail;

  std::queue<Edge const*> queue;
  queue.push(&current_edge);

  while (!queue.empty())
  {
    Edge const* edge = queue.front();
    queue.pop();

    if (edge->is_visited())
      continue;

    intersection_points.emplace_back(edge->intersection_point_);
    edge->visited();

    // Enqueue all exit edges reachable from faces for which this edge is an entry point.
    for (auto entry_face_iter = edge->entry_face_begin(); entry_face_iter != edge->entry_face_end(); ++entry_face_iter)
    {
      FaceId<n, T> const entry_face_id = *entry_face_iter;
      auto face_segment = face_segments.find(entry_face_id);
      EdgeId<n, T> const exit_point = face_segment->second.exit_point();
      auto exit_edge = edges.find(exit_point);
      ASSERT(exit_edge != edges.end());
      if (!exit_edge->second.is_visited())
        queue.push(&exit_edge->second);
    }
  }
}

} // namespace math
