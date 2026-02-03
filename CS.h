#pragma once

#include <cstdint>

namespace math {

// An extensible structural value type to be used as Coordinate System NTTP tag.
struct CS
{
  std::uint32_t id;

  constexpr bool operator==(CS const&) const = default;

  // Required in order for utils::to_string to work.
  char const* to_string() const;
};

namespace csid {
  inline constexpr CS painter{0};       // The space defined by the painter’s Current Transformation Matrix (CTM) at the instant you issue a painter->draw*() call.
  inline constexpr CS pixels{1};        // The final coordinate system of the window in pixels.
  inline constexpr CS plot{2};          // Used for cairowindow::Plot.
} // namespace csid

inline char const* CS::to_string() const
{
  if (id == csid::painter.id)
    return "cs::painter";
  else if (id == csid::pixels.id)
    return "cs::pixels";
  else if (id == csid::plot.id)
    return "cs::plot";
  return "cs:<unknown>";
}

} // namespace math
