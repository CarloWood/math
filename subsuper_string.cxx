#include "sys.h"
#include "subsuper_string.h"
#include <array>

namespace math {

static std::array<char const*, 10> u8_subscripts = {
  "₀", "₁", "₂", "₃", "₄", "₅", "₆", "₇", "₈", "₉"
};

static std::array<char const*, 10> u8_superscripts = {
  "⁰", "¹", "²", "³", "⁴", "⁵", "⁶", "⁷", "⁸", "⁹"
};

enum SubOrSuperScript {
  Subscript, Superscript
};

char const* get_subsuper(int digit, SubOrSuperScript what)
{
  return what == Subscript ? u8_subscripts[digit] : u8_superscripts[digit];
}

std::string to_subsuperscript(int value, SubOrSuperScript what)
{
  std::string result;
  if (value < 0)
  {
    result = what == Subscript ? "₋" : "⁻";
    value = -value;
  }
  int const value_div_10 = value / 10;
  int digit = 1;
  while (digit <= value_div_10)
    digit *= 10;
  do
  {
    result += get_subsuper(value / digit, what);
    digit /= 10;
  }
  while (digit);
  return result;
}

std::string to_subscript(int value)
{
  return to_subsuperscript(value, Subscript);
}

std::string to_superscript(int value)
{
  return to_subsuperscript(value, Superscript);
}

} // namespace math
