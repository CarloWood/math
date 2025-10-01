#include "sys.h"
#include "LinePiece.h"
#include "Direction.h"

namespace math {

Direction LinePiece::direction() const
{
  return Direction{from_, to_};
}

#ifdef CWDEBUG
void LinePiece::print_on(std::ostream& os) const
{
  os << "{from:" << from_ << ", to:" << to_ << "}";
}
#endif

} // namespace math
