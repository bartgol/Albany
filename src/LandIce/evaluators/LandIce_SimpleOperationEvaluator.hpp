//*****************************************************************//
//    Albany 3.0:  Copyright 2016 Sandia Corporation               //
//    This Software is released under the BSD license detailed     //
//    in the file "license.txt" in the top-level Albany directory  //
//*****************************************************************//

#ifndef LANDICE_SIMPLE_OPERATION_EVALUATOR_HPP
#define LANDICE_SIMPLE_OPERATION_EVALUATOR_HPP 1

#include "Phalanx_config.hpp"
#include "Phalanx_Evaluator_WithBaseImpl.hpp"
#include "Phalanx_Evaluator_Derived.hpp"
#include "Phalanx_MDField.hpp"

#include "Albany_MPL.hpp"
#include "Albany_Layouts.hpp"
#include "Albany_ScalarOrdinalTypes.hpp"

#include "PHAL_FieldsUtils.hpp"

#include "LandIce_SimpleOperation.hpp"

namespace LandIce
{

template<typename EvalT, typename Traits, typename Operation>
class SimpleOperationBase: public PHX::EvaluatorWithBaseImpl<Traits>,
                           public PHX::EvaluatorDerived<EvalT, Traits>
{
public:
  SimpleOperationBase (const Teuchos::ParameterList& p,
                       const Teuchos::RCP<Albany::Layouts>& dl);

  void postRegistrationSetup (typename Traits::SetupData,
                              PHX::FieldManager<Traits>&) {}

protected:
  using RT = RealType;
  using ST = typename EvalT::ScalarT;
  using PT = typename EvalT::ParamScalarT;
  using MT = typename EvalT::MeshScalarT;
 
  // Input:
  PHX::MDField<const RT> in_rt;
  PHX::MDField<const ST> in_st;
  PHX::MDField<const PT> in_pt;
  PHX::MDField<const MT> in_mt;

  // Output:
  PHX::MDField<const RT> out_rt;
  PHX::MDField<const ST> out_st;
  PHX::MDField<const PT> out_pt;
  PHX::MDField<const MT> out_mt;

  // The scalar type of the input/output field
  PHAL::FieldScalarType inout_fst;

  // The operation
  Operation             op;
};

// =================== Specializations For Unary Operations ================= //

template<typename EvalT, typename Traits, typename UnaryOperation>
class SimpleUnaryOperation : public SimpleOperationBase<EvalT,Traits,UnaryOperation>
{
public:
  SimpleUnaryOperation (const Teuchos::ParameterList& p,
                        const Teuchos::RCP<Albany::Layouts>& dl);

  void evaluateFields(typename Traits::EvalData d);
};

template<typename EvalT, typename Traits>
class UnaryScaleOp : public SimpleUnaryOperation<EvalT,Traits,UnaryOps::Scale>
{
public:
  UnaryScaleOp (const Teuchos::ParameterList& p,
                const Teuchos::RCP<Albany::Layouts>& dl);
};

template<typename EvalT, typename Traits>
class UnaryLogOp : public SimpleUnaryOperation<EvalT,Traits,UnaryOps::Log>
{
public:
  UnaryLogOp (const Teuchos::ParameterList& p,
              const Teuchos::RCP<Albany::Layouts>& dl);
};

template<typename EvalT, typename Traits>
class UnaryExpOp : public SimpleUnaryOperation<EvalT,Traits,UnaryOps::Exp>
{
public:
  UnaryExpOp (const Teuchos::ParameterList& p,
              const Teuchos::RCP<Albany::Layouts>& dl);
};

template<typename EvalT, typename Traits>
class UnaryLowPassOp : public SimpleUnaryOperation<EvalT,Traits,UnaryOps::LowPass>
{
public:
  UnaryLowPassOp (const Teuchos::ParameterList& p,
                  const Teuchos::RCP<Albany::Layouts>& dl);
};

template<typename EvalT, typename Traits>
class UnaryHighPassOp : public SimpleUnaryOperation<EvalT,Traits,UnaryOps::HighPass>
{
public:
  UnaryHighPassOp (const Teuchos::ParameterList& p,
                   const Teuchos::RCP<Albany::Layouts>& dl);
};

template<typename EvalT, typename Traits>
class UnaryBandPassOp : public SimpleUnaryOperation<EvalT,Traits,UnaryOps::BandPass>
{
public:
  UnaryBandPassOp (const Teuchos::ParameterList& p,
                   const Teuchos::RCP<Albany::Layouts>& dl);
};

// =================== Specializations For Binary Operations ================= //

template<typename EvalT, typename Traits, typename BinaryOperation>
class SimpleBinaryOperation : public SimpleOperationBase<EvalT,Traits,BinaryOperation>
{
public:
  SimpleBinaryOperation (const Teuchos::ParameterList& p,
                         const Teuchos::RCP<Albany::Layouts>& dl);

  void postRegistrationSetup (typename Traits::SetupData d,
                              PHX::FieldManager<Traits>& fm);

  void evaluateFields(typename Traits::EvalData d);

private:
  template<typename F1, typename F2, typename F3>
  void iterate (const F1& f1, const F2& f2, F3& f3) {
    auto in1 = mdFieldIterator(f1);
    auto in2 = mdFieldIterator(f2);
    auto out = mdFieldIterator(f3);
    for (; !out.done(); ++in1, ++in2, ++out) {
      *out = this->op(*in1,*in2);
    }
  }
  using RT = RealType;
  using ST = typename EvalT::ScalarT;
  using PT = typename EvalT::ParamScalarT;
  using MT = typename EvalT::MeshScalarT;

  PHX::MDField<const RT> field_rt;
  PHX::MDField<const RT> field_pt;
  PHX::MDField<const RT> field_mt;
  PHX::MDField<const RT> field_st;

  PHAL::FieldScalarType field_fst;
};

template<typename EvalT, typename Traits>
class BinaryScaleOp : public SimpleBinaryOperation<EvalT,Traits,BinaryOps::Scale>
{
public:
  BinaryScaleOp (const Teuchos::ParameterList& p,
                 const Teuchos::RCP<Albany::Layouts>& dl);
};

template<typename EvalT, typename Traits>
class BinarySumOp : public SimpleBinaryOperation<EvalT,Traits,BinaryOps::Sum>
{
public:
  BinarySumOp (const Teuchos::ParameterList& p,
               const Teuchos::RCP<Albany::Layouts>& dl);
};

template<typename EvalT, typename Traits>
class BinaryLogOp : public SimpleBinaryOperation<EvalT,Traits,BinaryOps::Log>
{
public:
  BinaryLogOp (const Teuchos::ParameterList& p,
               const Teuchos::RCP<Albany::Layouts>& dl);
};

template<typename EvalT, typename Traits>
class BinaryExpOp : public SimpleBinaryOperation<EvalT,Traits,BinaryOps::Exp>
{
public:
  BinaryExpOp (const Teuchos::ParameterList& p,
               const Teuchos::RCP<Albany::Layouts>& dl);
};

template<typename EvalT, typename Traits>
class BinaryLowPassOp : public SimpleBinaryOperation<EvalT,Traits,BinaryOps::LowPass>
{
public:
  BinaryLowPassOp (const Teuchos::ParameterList& p,
                   const Teuchos::RCP<Albany::Layouts>& dl);
};

template<typename EvalT, typename Traits>
class BinaryHighPassOp : public SimpleBinaryOperation<EvalT,Traits,BinaryOps::HighPass>
{
public:
  BinaryHighPassOp (const Teuchos::ParameterList& p,
                   const Teuchos::RCP<Albany::Layouts>& dl);
};

template<typename EvalT, typename Traits>
class BinaryBandPassFixedLowerOp : public SimpleBinaryOperation<EvalT,Traits,BinaryOps::BandPassFixedLower>
{
public:
  BinaryBandPassFixedLowerOp (const Teuchos::ParameterList& p,
                              const Teuchos::RCP<Albany::Layouts>& dl);
};

template<typename EvalT, typename Traits>
class BinaryBandPassFixedUpperOp : public SimpleBinaryOperation<EvalT,Traits,BinaryOps::BandPassFixedUpper>
{
public:
  BinaryBandPassFixedUpperOp (const Teuchos::ParameterList& p,
                              const Teuchos::RCP<Albany::Layouts>& dl);
};

// =================== Specializations For Ternary Operations ================= //

template<typename EvalT, typename Traits, typename TernaryOperation>
class SimpleTernaryOperation : public SimpleOperationBase<EvalT,Traits,TernaryOperation>
{
public:
  SimpleTernaryOperation (const Teuchos::ParameterList& p,
                          const Teuchos::RCP<Albany::Layouts>& dl);

  void postRegistrationSetup (typename Traits::SetupData d,
                              PHX::FieldManager<Traits>& fm);

  void evaluateFields(typename Traits::EvalData d);

private:
  template<typename F1, typename F2, typename F3, typename F4>
  void iterate (const F1& f1, const F2& f2, const F3& f3, F4& f4) {
    auto in1 = mdFieldIterator(f1);
    auto in2 = mdFieldIterator(f2);
    auto in3 = mdFieldIterator(f3);
    auto out = mdFieldIterator(f4);
    for (; !out.done(); ++in1, ++in2, ++in3, ++out) {
      *out = this->op(*in1,*in2,*in3);
    }
  }

  using RT = RealType;
  using ST = typename EvalT::ScalarT;
  using PT = typename EvalT::ParamScalarT;
  using MT = typename EvalT::MeshScalarT;

  PHX::MDField<const RT> field1_rt;
  PHX::MDField<const RT> field1_pt;
  PHX::MDField<const RT> field1_mt;
  PHX::MDField<const RT> field1_st;
  PHX::MDField<const RT> field2_rt;
  PHX::MDField<const RT> field2_pt;
  PHX::MDField<const RT> field2_mt;
  PHX::MDField<const RT> field2_st;

  PHAL::FieldScalarType field_fst;
};

template<typename EvalT, typename Traits>
class TernaryBandPassOp : public SimpleTernaryOperation<EvalT,Traits,TernaryOps::BandPass>
{
public:
  TernaryBandPassOp (const Teuchos::ParameterList& p,
                     const Teuchos::RCP<Albany::Layouts>& dl);
};

} // Namespace LandIce

#endif // LANDICE_SIMPLE_OPERATION_EVALUATOR_HPP
