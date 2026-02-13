#pragma once

#include "Point.h"
#include "Size.h"
#include "math/Transform.h"

namespace math::cs {

template<CS from_cs, CS to_cs, bool inverted, math::AffineTransformConcept AffineTransformBackend>
Point<to_cs> operator*(Point<from_cs> const& from, math::Transform<from_cs, to_cs, inverted, AffineTransformBackend> const& transform)
{
  static_assert(!inverted, "Do not multiply with a Transform<..., inverted = true, ...> type. Pre-calculate the inverse first by calling Transform::inverted().");
  return transform.map_point(from);
}

template<CS from_cs, CS to_cs, bool inverted, math::AffineTransformConcept AffineTransformBackend>
Vector<to_cs> operator*(Vector<from_cs> const& from, math::Transform<from_cs, to_cs, inverted, AffineTransformBackend> const& transform)
{
  static_assert(!inverted, "Do not multiply with a Transform<..., inverted = true, ...> type. Pre-calculate the inverse first by calling Transform::inverted().");
  return transform.map_vector(from);
}

template<CS from_cs, CS to_cs, bool inverted, math::AffineTransformConcept AffineTransformBackend>
Direction<to_cs> operator*(Direction<from_cs> const& from, math::Transform<from_cs, to_cs, inverted, AffineTransformBackend> const& transform)
{
  static_assert(!inverted, "Do not multiply with a Transform<..., inverted = true, ...> type. Pre-calculate the inverse first by calling Transform::inverted().");
  return transform.map_direction(from);
}

template<CS from_cs, CS to_cs, bool inverted, math::AffineTransformConcept AffineTransformBackend>
Size<to_cs> operator*(Size<from_cs> const& from, math::Transform<from_cs, to_cs, inverted, AffineTransformBackend> const& transform)
{
  static_assert(!inverted, "Do not multiply with a Transform<..., inverted = true, ...> type. Pre-calculate the inverse first by calling Transform::inverted().");
  return transform.map_size(from);
}

} // namespace math::cs
