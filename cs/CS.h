#pragma once

#include <cstdint>
#include "utils/has_print_on.h"

namespace math {
// This class defines a print_on method.
using utils::has_print_on::operator<<;

// An extensible structural value type to be used as Coordinate System NTTP tag.
struct CS
{
  char const* id;

  constexpr bool operator==(CS const&) const = default;

  // Required in order for utils::to_string to work.
  char const* to_string() const { return id; }

  // Required in order for utils::has_print_on to work.
  void print_on(std::ostream& os) const;
};

#define DECLARE_CSID(x) \
  inline constexpr char x##_name[] = "csid::" #x; \
  inline constexpr math::CS x{x##_name};

namespace csid {
  DECLARE_CSID(painter);        // The space defined by the painter’s Current Transformation Matrix (CTM) at the instant you issue a painter->draw*() call.
  DECLARE_CSID(pixels);         // The final coordinate system of the window in pixels.
  DECLARE_CSID(plot);           // Used for cairowindow::Plot.
} // namespace csid

} // namespace math
