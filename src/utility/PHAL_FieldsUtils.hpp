#ifndef PHAL_FIELDS_UTILS_HPP
#define PHAL_FIELDS_UTILS_HPP

#include <type_traits>
#include <string>

namespace PHAL {

/*
 * Enums and utils to automatize MDFields specification
 * 
 * The enums allow to specify the field scalar type, rank, and location.
 * These information can then be used to construct the proper PHX::MDField,
 * with the proper scalar type and MDA layout
 */

// ================== Enums ================ //

// Define an invalid string in one place (to avoid typos-related bugs).
constexpr const char* INVALID_STR = "__INVALID__";

// The type of the field scalar
enum class FieldScalarType : int {
  Real        = 0,
  MeshScalar  = 1,
  ParamScalar = 2,
  Scalar      = 3
};

// Mesh entity where a field is located
enum class FieldLocation : int {
  Cell,
  Node,
  QuadPoint
};

// Type of field (scalar, vector, gradient, tensor)
// Note: gradient is just a vector with length equal to the mesh dimension
enum class FieldRankType : int {
  Scalar,
  Vector,
  Gradient,
  Tensor
};

// Teuchos requires Teuchos::is_printable<T>::type and Teuchos::is_comparable<T>::type
// to be std::true_type in order for T to be storable in a parameter list.
// However, strong enums are not printable by default, so we need an overload for op<<
template<typename T>
typename
std::enable_if<std::is_same<T,FieldLocation>::value ||
               std::is_same<T,FieldRankType>::value,
               std::ostream&>::type
operator<< (std::ostream& out, const T t) {
  out << e2str(t);
  return out;
}

// ================= Utilities ================= //

// Extract underlying integer value from an enum
template<typename EnumType>
typename std::underlying_type<EnumType>::type etoi (const EnumType e) {
  return static_cast<typename std::underlying_type<EnumType>::type>(e);
}

// Convert enum values to strings
inline std::string e2str (const FieldLocation e) {
  switch (e) {
    case FieldLocation::Node:       return "Node";
    case FieldLocation::QuadPoint:  return "QuadPoint";
    case FieldLocation::Cell:       return "Cell";
  }
  return INVALID_STR;
}

inline std::string e2str (const FieldRankType rank) {
  switch (rank) {
    case FieldRankType::Scalar:     return "Scalar";
    case FieldRankType::Vector:     return "Vector";
    case FieldRankType::Gradient:   return "Gradient";
    case FieldRankType::Tensor:     return "Tensor";
  }

  return INVALID_STR;
}

inline std::string e2str (const FieldScalarType e) {
  switch (e) {
    case FieldScalarType::Scalar:       return "Scalar";      break;
    case FieldScalarType::MeshScalar:   return "MeshScalar";  break;
    case FieldScalarType::ParamScalar:  return "ParamScalar"; break;
    case FieldScalarType::Real:         return "Real";        break;
  }
  return INVALID_STR;
}


// Operators to compute the 'strongest' scalar type, that is,
// the scalar type resulting from something like st1+st2.
inline FieldScalarType& operator|= (FieldScalarType& st1,
                                   const FieldScalarType& st2)
{
  // Return the 'strongest' scalar type. In the enum above, they are ordered per 'strength'.
  // The idea is that the assignment of a scalar type A from a scalar type B is legal if
  // A is 'stronger' than B.

  auto st1_int = etoi(st1);
  auto st2_int = etoi(st2);

  if (st2_int>st1_int) {
    st1 = st2;
  }

  return st1;
}

inline FieldScalarType operator| (const FieldScalarType& st1,
                                  const FieldScalarType& st2)
{
  FieldScalarType st3 = st1;
  st3 |= st2;
  return st3;
}

} // namespace PHAL

#endif // PHAL_FIELDS_UTILS_HPP
