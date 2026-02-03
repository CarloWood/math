#include "sys.h"
#include "CS.h"
#include <iostream>

namespace math {

void CS::print_on(std::ostream& os) const
{
  std::string name = to_string();
  os << name;
}

} // namespace math
