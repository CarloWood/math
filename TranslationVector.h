#pragma once

#include "CS.h"
#include "Vector.h"

namespace math {

// TranslationVector
//
// Type of the argument of Transform::tranlate.
//
// Classes that can be used as translation type should add an implicit conversion operator.
// For example:
//
// template<CS cs>
// class Foo {
//  private:
//   /* state variables */
//  public:
//   // Implicit converstion to math::TranslationVector<cs>, so we can pass a Foo<cs> to math::Transform<>::translate.
//   operator math::TranslationVector<cs>() const
//   {
//     return math::TranslationVector<cs>::create_from_cs_values(/* state variables */);
//   }
// };
//
template<CS cs>
class TranslationVector
{
 private:
  math::Vector<2> translation_;

 private:
  TranslationVector(Vector<2> const& translation) : translation_(translation) { }
  TranslationVector(double dx, double dy) : translation_(dx, dy) { }

 public:
  double dx() const { return translation_.x(); }
  double dy() const { return translation_.y(); }

  // Allow explicit conversion of non-coordinate-system values to be used to construct a TranslationVector.
  // This should really only be used in the implicit conversion operator of a templated type<cs>.
  static TranslationVector create_from_cs_values(Vector<2> const& vector_cs) { return {vector_cs}; }
  static TranslationVector create_from_cs_values(double dx_cs, double dy_cs) { return {dx_cs, dy_cs}; }

  friend TranslationVector operator*(double scale, TranslationVector const& tv)
  {
    return {scale * tv.translation_};
  }
};

} // namespace math
