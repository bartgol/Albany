#ifndef ALBANY_MPL_HPP
#define ALBANY_MPL_HPP

#include <Phalanx_MDField.hpp>

#include <tuple>

namespace Albany {

// A type to represent nothing. Notice that using 'void'
// *might* be a problem if 'void' is a legitimate type
// to be stored in some MPL type.
struct NullType {};

// template<typename... T>
// struct TypeSet;

// template<typename T>
// struct TypeSet<T> {
//   template<typename S>
//   static constexpr bool has () {
//     return std::is_same<T,S>::value;
//   }

//   using type = std::tuple<T>;

//   static constexpr int size = 1;
// };

// template<typename T, typename... Ts>
// struct TypeSet<T,Ts...> {
//   using tail = TypeSet<Ts...>;

//   template<typename S>
//   static constexpr bool has () {
//     return std::is_same<T,S>::value || tail::template has<S>();
//   }

//   using type = typename std::conditional<tail::template has<T>(),
//         typename tail::type,
//         std::tuple<T,Ts...>
//     >::type;

//   static constexpr int size = 1 + tail::size;
// };

// A compile-time map of types, with each entry being
// an std::pair<K,V>. The same Key may be present
// multiple times.
template<typename... T>
struct CT_Map;

template<>
struct CT_Map<> {
  template<typename KT>
  using value_type = NullType;

  template<typename KT>
  static constexpr bool has_key () { return false; }
};

template<typename Key, typename Val, typename... T>
struct CT_Map<std::pair<Key,Val>,T...> {
  using tail = CT_Map<T...>;

  template<typename KT>
  static constexpr bool has_key () {
    return std::is_same<Key,KT>::value || tail::template has_key<KT>();
  }

  // Find the type corresponding to a particular key
  template<typename KT>
  using value_type = typename std::conditional<
                        std::is_same<KT,Key>::value,
                        Val,
                        typename tail::template value_type<KT>>::type;
                        
};

// ScalarT -> MDField map
template<typename EvalT, typename... Tags>
struct FieldMap {
  // The scalar types provided by EvalT
  using ST = typename EvalT::ScalarT;
  using PT = typename EvalT::ParamScalarT;
  using MT = typename EvalT::MeshScalarT;

  template<typename Scalar>
  using FT = PHX::MDField<Scalar,Tags...>;

  template<typename Scalar>
  using entry = std::pair<Scalar,FT<Scalar>>;

  using type = CT_Map<entry<ST>,entry<PT>,entry<MT>>;

  // using FieldST = FT<ST>;
};

} // namespace Albany

#endif // ALBANY_MPL_HPP
