#ifndef PHAL_MD_FIELD_HPP
#define PHAL_MD_FIELD_HPP

#include "Albany_MPL.hpp"
#include "Albany_ScalarOrdinalTypes.hpp"
#include "PHAL_FieldsUtils.hpp"

#include <Phalanx_MDField.hpp>

namespace PHAL {

template<typename Scalar1, typename... ScalarT>
struct FieldTypes {
  template<typename S>
  using FT = PHX::MDField<S>;

  using tail = typename FieldTypes<ScalarT...>::type;
  using head = std::tuple<FT<Scalar1>>;

  using type = decltype(std::tuple_cat(std::declval<head>(),std::declval<tail>()));
};

template<typename Scalar>
struct FieldTypes<Scalar> {
  template<typename S>
  using FT = PHX::MDField<S>;
  using type = std::tuple<FT<Scalar>>;
};

template<typename... S>
struct FieldTypes<std::tuple<S...>>{
  using type = typename FieldTypes<S...>::type;
};

template<typename EvalT>
class MDField {
public:
  using RT = RealType;
  using MT = typename EvalT::MeshScalarT;
  using PT = typename EvalT::ParamScalarT;
  using ST = typename EvalT::ScalarT;

  MDField (const FieldScalarType fst, const std::string& name,
           const Teuchos::RCP<PHX::DataLayout>& dl);

  Teuchos::RCP<const PHX::FieldTag> fieldTag () const;

private:

  // A tuple containing all *different* types in the sequence RT,MT,PT,ST.
  using ScalarTypes = typename Albany::TypeSet<RT,MT,PT,ST>::type;
  // A tuple containing all MDFields corresponding to all Scalar types
  using FieldsT = typename FieldTypes<ScalarTypes>::type;


  FieldsT fields;
};

} // namespace PHAL

#endif // PHAL_MD_FIELD_HPP
