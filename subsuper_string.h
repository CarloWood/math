#pragma once

#include <string>

// Convert integers to sub- and/or superscript strings.
//
namespace math {

std::string to_subscript(int value);
std::string to_superscript(int value);

} // namespace math
