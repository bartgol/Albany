//*****************************************************************//
//    Albany 3.0:  Copyright 2016 Sandia Corporation               //
//    This Software is released under the BSD license detailed     //
//    in the file "license.txt" in the top-level Albany directory  //
//*****************************************************************//

#include "Phalanx_DataLayout.hpp"
#include "Phalanx_Print.hpp"

#include "LandIce_SimpleOperationEvaluator.hpp"
#include "LandIce_ProblemUtils.hpp"

#include "PHAL_Utilities.hpp"

namespace LandIce {

//**********************************************************************
template<typename EvalT, typename Traits, typename Operation>
SimpleOperationBase<EvalT, Traits, Operation>::
SimpleOperationBase (const Teuchos::ParameterList& p,
                     const Teuchos::RCP<Albany::Layouts>& /* dl */)
{
  using FST = PHAL::FieldScalarType;

  std::string fieldInName  = p.get<std::string> ("Input Field Name");
  std::string fieldOutName = p.get<std::string> ("Output Field Name");
  inout_fst = p.get<FST> ("Input Field Scalar Type");

  Teuchos::RCP<PHX::DataLayout> layout = p.get<Teuchos::RCP<PHX::DataLayout>>("Field Layout");
  TEUCHOS_TEST_FOR_EXCEPTION (layout.is_null(), std::runtime_error, "Error! Input layout is null.\n");

  switch (inout_fst) {
    case FST::Real:
      in_rt = decltype(in_rt)(fieldInName, layout);
      out_rt = decltype(out_rt)(fieldOutName, layout);
      this->addDependentField(in_rt);
      this->addEvaluatedField(out_rt);
      break;
    case FST::MeshScalar:
      in_mt = decltype(in_mt)(fieldInName, layout);
      out_mt = decltype(out_mt)(fieldOutName, layout);
      this->addDependentField(in_mt);
      this->addEvaluatedField(out_mt);
      break;
    case FST::ParamScalar:
      in_pt = decltype(in_pt)(fieldInName, layout);
      out_pt = decltype(out_pt)(fieldOutName, layout);
      this->addDependentField(in_pt);
      this->addEvaluatedField(out_pt);
      break;
    case FST::Scalar:
      in_st = decltype(in_st)(fieldInName, layout);
      out_st = decltype(out_st)(fieldOutName, layout);
      this->addDependentField(in_st);
      this->addEvaluatedField(out_st);
      break;
  }

  this->setName("SimpleOperationBase"+PHX::print<EvalT>());
}

//**********************************************************************
template<typename EvalT, typename Traits, typename UnaryOperation>
SimpleUnaryOperation<EvalT, Traits, UnaryOperation>::
SimpleUnaryOperation (const Teuchos::ParameterList& p,
                      const Teuchos::RCP<Albany::Layouts>& dl) :
  SimpleOperationBase<EvalT,Traits,UnaryOperation> (p,dl)
{
  this->op.setup(p);
}

//**********************************************************************
template<typename EvalT, typename Traits, typename UnaryOperation>
void SimpleUnaryOperation<EvalT, Traits, UnaryOperation>::
evaluateFields (typename Traits::EvalData /* workset */)
{
  using RT = RealType;
  using ST = typename EvalT::ScalarT;
  using PT = typename EvalT::ParamScalarT;
  using MT = typename EvalT::MeshScalarT;
  using FST = PHAL::FieldScalarType;

  switch (this->inout_fst) {
    case FST::Real:
    {
      auto in  = mdFieldIterator(this->in_rt);
      auto out = mdFieldIterator(this->out_rt);
      for (; !in.done(); ++in, ++out) {
        *out = this->op(*in);
      }
    }
    case FST::MeshScalar:
    {
      auto in  = mdFieldIterator(this->in_mt);
      auto out = mdFieldIterator(this->out_mt);
      for (; !in.done(); ++in, ++out) {
        *out = this->op(*in);
      }
    }
    case FST::ParamScalar:
    {
      auto in  = mdFieldIterator(this->in_pt);
      auto out = mdFieldIterator(this->out_pt);
      for (; !in.done(); ++in, ++out) {
        *out = this->op(*in);
      }
    }
    case FST::Scalar:
    {
      auto in  = mdFieldIterator(this->in_st);
      auto out = mdFieldIterator(this->out_st);
      for (; !in.done(); ++in, ++out) {
        *out = this->op(*in);
      }
    }
  }
}

//**********************************************************************
template<typename EvalT, typename Traits>
UnaryScaleOp<EvalT, Traits>::
UnaryScaleOp (const Teuchos::ParameterList& p,
              const Teuchos::RCP<Albany::Layouts>& dl) :
  SimpleUnaryOperation<EvalT,Traits,UnaryOps::Scale>(p,dl) {}

template<typename EvalT, typename Traits>
UnaryLogOp<EvalT, Traits>::
UnaryLogOp (const Teuchos::ParameterList& p,
            const Teuchos::RCP<Albany::Layouts>& dl) :
  SimpleUnaryOperation<EvalT,Traits,UnaryOps::Log>(p,dl) {}

template<typename EvalT, typename Traits>
UnaryExpOp<EvalT, Traits>::
UnaryExpOp (const Teuchos::ParameterList& p,
            const Teuchos::RCP<Albany::Layouts>& dl) :
  SimpleUnaryOperation<EvalT,Traits,UnaryOps::Exp>(p,dl) {}

template<typename EvalT, typename Traits>
UnaryLowPassOp<EvalT, Traits>::
UnaryLowPassOp (const Teuchos::ParameterList& p,
                const Teuchos::RCP<Albany::Layouts>& dl) :
  SimpleUnaryOperation<EvalT,Traits,UnaryOps::LowPass>(p,dl) {}

template<typename EvalT, typename Traits>
UnaryHighPassOp<EvalT, Traits>::
UnaryHighPassOp (const Teuchos::ParameterList& p,
                 const Teuchos::RCP<Albany::Layouts>& dl) :
  SimpleUnaryOperation<EvalT,Traits,UnaryOps::HighPass>(p,dl) {}

template<typename EvalT, typename Traits>
UnaryBandPassOp<EvalT, Traits>::
UnaryBandPassOp (const Teuchos::ParameterList& p,
                 const Teuchos::RCP<Albany::Layouts>& dl) :
  SimpleUnaryOperation<EvalT,Traits,UnaryOps::BandPass>(p,dl) {}

//**********************************************************************
template<typename EvalT, typename Traits, typename BinaryOperation>
SimpleBinaryOperation<EvalT, Traits, BinaryOperation>::
SimpleBinaryOperation (const Teuchos::ParameterList& p,
                       const Teuchos::RCP<Albany::Layouts>& dl) :
  SimpleOperationBase<EvalT,Traits,BinaryOperation> (p,dl)
{
  using FST = PHAL::FieldScalarType;

  field_fst = p.get<FST>("Parameter Field 1 Scalar Type");
  auto name   = p.get<std::string>("Parameter Field 1");
  auto layout = p.get<Teuchos::RCP<PHX::DataLayout>>("Field Layout");
  switch (field_fst) {
    case FST::Real:
      field_rt = decltype(field_rt) (name,layout);
      this->addDependentField(field_rt);
      break;
    case FST::MeshScalar:
      field_mt = decltype(field_mt) (name, layout);
      this->addDependentField(field_mt);
      break;
    case FST::ParamScalar:
      field_pt = decltype(field_pt) (name, layout);
      this->addDependentField(field_pt);
      break;
    case FST::Scalar:
      field_st = decltype(field_st) (name, layout);
      this->addDependentField(field_st);
      break;
  }
}

//**********************************************************************
template<typename EvalT, typename Traits, typename BinaryOperation>
void SimpleBinaryOperation<EvalT, Traits, BinaryOperation>::
postRegistrationSetup (typename Traits::SetupData d,
                       PHX::FieldManager<Traits>& fm)
{
  SimpleOperationBase<EvalT,Traits,BinaryOperation>::postRegistrationSetup(d,fm);
}

//**********************************************************************
template<typename EvalT, typename Traits, typename BinaryOperation>
void SimpleBinaryOperation<EvalT, Traits, BinaryOperation>::
evaluateFields (typename Traits::EvalData /* workset */)
{
  using FST = PHAL::FieldScalarType;
  switch (this->inout_fst) {
    case FST::Real:
      switch(field_fst) {
        case FST::Real:
          iterate (this->in_rt, field_rt, this->out_rt);
          break;
        case FST::MeshScalar:
          iterate (this->in_rt, field_mt, this->out_rt);
          break;
        case FST::ParamScalar:
          iterate (this->in_rt, field_pt, this->out_rt);
          break;
        case FST::Scalar:
          iterate (this->in_rt, field_st, this->out_rt);
          break;
      } break;
    case FST::MeshScalar:
      switch(field_fst) {
        case FST::Real:
          iterate (this->in_mt, field_rt, this->out_mt);
          break;
        case FST::MeshScalar:
          iterate (this->in_mt, field_mt, this->out_mt);
          break;
        case FST::ParamScalar:
          iterate (this->in_mt, field_pt, this->out_mt);
          break;
        case FST::Scalar:
          iterate (this->in_mt, field_st, this->out_mt);
          break;
      } break;
    case FST::ParamScalar:
      switch(field_fst) {
        case FST::Real:
          iterate (this->in_pt, field_rt, this->out_pt);
          break;
        case FST::MeshScalar:
          iterate (this->in_pt, field_mt, this->out_pt);
          break;
        case FST::ParamScalar:
          iterate (this->in_pt, field_pt, this->out_pt);
          break;
        case FST::Scalar:
          iterate (this->in_pt, field_st, this->out_pt);
          break;
      } break;
    case FST::Scalar:
      switch(field_fst) {
        case FST::Real:
          iterate (this->in_st, field_rt, this->out_st);
          break;
        case FST::MeshScalar:
          iterate (this->in_st, field_mt, this->out_st);
          break;
        case FST::ParamScalar:
          iterate (this->in_st, field_pt, this->out_st);
          break;
        case FST::Scalar:
          iterate (this->in_st, field_st, this->out_st);
          break;
      } break;
  }
}

//**********************************************************************
template<typename EvalT, typename Traits>
BinaryScaleOp<EvalT, Traits>::
BinaryScaleOp (const Teuchos::ParameterList& p,
               const Teuchos::RCP<Albany::Layouts>& dl) :
  SimpleBinaryOperation<EvalT,Traits,BinaryOps::Scale>(p,dl) {}

template<typename EvalT, typename Traits>
BinarySumOp<EvalT, Traits>::
BinarySumOp (const Teuchos::ParameterList& p,
             const Teuchos::RCP<Albany::Layouts>& dl) :
  SimpleBinaryOperation<EvalT,Traits,BinaryOps::Sum>(p,dl) {}

template<typename EvalT, typename Traits>
BinaryLogOp<EvalT, Traits>::
BinaryLogOp (const Teuchos::ParameterList& p,
             const Teuchos::RCP<Albany::Layouts>& dl) :
  SimpleBinaryOperation<EvalT,Traits,BinaryOps::Log>(p,dl) {}

template<typename EvalT, typename Traits>
BinaryExpOp<EvalT, Traits>::
BinaryExpOp (const Teuchos::ParameterList& p,
             const Teuchos::RCP<Albany::Layouts>& dl) :
  SimpleBinaryOperation<EvalT,Traits,BinaryOps::Exp>(p,dl) {}

template<typename EvalT, typename Traits>
BinaryLowPassOp<EvalT, Traits>::
BinaryLowPassOp (const Teuchos::ParameterList& p,
                 const Teuchos::RCP<Albany::Layouts>& dl) :
  SimpleBinaryOperation<EvalT,Traits,BinaryOps::LowPass>(p,dl) {}

template<typename EvalT, typename Traits>
BinaryHighPassOp<EvalT, Traits>::
BinaryHighPassOp (const Teuchos::ParameterList& p,
                 const Teuchos::RCP<Albany::Layouts>& dl) :
  SimpleBinaryOperation<EvalT,Traits,BinaryOps::HighPass>(p,dl) {}

template<typename EvalT, typename Traits>
BinaryBandPassFixedLowerOp<EvalT, Traits>::
BinaryBandPassFixedLowerOp (const Teuchos::ParameterList& p,
                            const Teuchos::RCP<Albany::Layouts>& dl) :
  SimpleBinaryOperation<EvalT,Traits,BinaryOps::BandPassFixedLower>(p,dl) {}

template<typename EvalT, typename Traits>
BinaryBandPassFixedUpperOp<EvalT, Traits>::
BinaryBandPassFixedUpperOp (const Teuchos::ParameterList& p,
                            const Teuchos::RCP<Albany::Layouts>& dl) :
  SimpleBinaryOperation<EvalT,Traits,BinaryOps::BandPassFixedUpper>(p,dl) {}

//**********************************************************************
template<typename EvalT, typename Traits, typename TernaryOperation>
SimpleTernaryOperation<EvalT, Traits, TernaryOperation>::
SimpleTernaryOperation (const Teuchos::ParameterList& p,
                        const Teuchos::RCP<Albany::Layouts>& dl) :
  SimpleOperationBase<EvalT,Traits,TernaryOperation> (p,dl)
{
  using FST = PHAL::FieldScalarType;

  field_fst = p.get<FST>("Parameter Field 1 Scalar Type");
  auto name1   = p.get<std::string>("Parameter Field 1");
  auto name2   = p.get<std::string>("Parameter Field 2");
  auto layout = p.get<Teuchos::RCP<PHX::DataLayout>>("Field Layout");
  switch (field_fst) {
    case FST::Real:
      field1_rt = decltype(field1_rt) (name1,layout);
      field1_rt = decltype(field2_rt) (name1,layout);
      this->addDependentField(field1_rt);
      this->addDependentField(field2_rt);
      break;
    case FST::MeshScalar:
      field1_mt = decltype(field1_mt) (name1, layout);
      field1_mt = decltype(field2_mt) (name1, layout);
      this->addDependentField(field1_mt);
      this->addDependentField(field2_mt);
      break;
    case FST::ParamScalar:
      field1_mt = decltype(field1_mt) (name1, layout);
      field1_mt = decltype(field2_mt) (name1, layout);
      this->addDependentField(field1_mt);
      this->addDependentField(field2_mt);
      break;
    case FST::Scalar:
      field1_st = decltype(field1_st) (name1, layout);
      field1_st = decltype(field2_st) (name1, layout);
      this->addDependentField(field1_st);
      this->addDependentField(field2_st);
      break;
  }
}

//**********************************************************************
template<typename EvalT, typename Traits, typename TernaryOperation>
void SimpleTernaryOperation<EvalT, Traits, TernaryOperation>::
postRegistrationSetup (typename Traits::SetupData d,
                       PHX::FieldManager<Traits>& fm)
{
  SimpleOperationBase<EvalT,Traits,TernaryOperation>::postRegistrationSetup(d,fm);
}

//**********************************************************************
template<typename EvalT, typename Traits, typename TernaryOperation>
void SimpleTernaryOperation<EvalT, Traits, TernaryOperation>::
evaluateFields (typename Traits::EvalData /* workset */)
{
  using FST = PHAL::FieldScalarType;
  switch (this->inout_fst) {
    case FST::Real:
      switch(field_fst) {
        case FST::Real:
          iterate (this->in_rt, field1_rt, field2_rt, this->out_rt);
          break;                           
        case FST::MeshScalar:              
          iterate (this->in_rt, field1_mt, field2_mt, this->out_rt);
          break;                           
        case FST::ParamScalar:             
          iterate (this->in_rt, field1_pt, field2_pt, this->out_rt);
          break;                           
        case FST::Scalar:                  
          iterate (this->in_rt, field1_st, field2_st, this->out_rt);
          break;
      } break;
    case FST::MeshScalar:
      switch(field_fst) {
        case FST::Real:
          iterate (this->in_mt, field1_rt, field2_rt, this->out_mt);
          break;                           
        case FST::MeshScalar:              
          iterate (this->in_mt, field1_mt, field2_mt, this->out_mt);
          break;                           
        case FST::ParamScalar:             
          iterate (this->in_mt, field1_pt, field2_pt, this->out_mt);
          break;                           
        case FST::Scalar:                  
          iterate (this->in_mt, field1_st, field2_st, this->out_mt);
          break;
      } break;
    case FST::ParamScalar:
      switch(field_fst) {
        case FST::Real:
          iterate (this->in_pt, field1_rt, field2_rt, this->out_pt);
          break;                           
        case FST::MeshScalar:              
          iterate (this->in_pt, field1_mt, field2_mt, this->out_pt);
          break;                           
        case FST::ParamScalar:             
          iterate (this->in_pt, field1_pt, field2_pt, this->out_pt);
          break;                           
        case FST::Scalar:                  
          iterate (this->in_pt, field1_st, field2_st, this->out_pt);
          break;
      } break;
    case FST::Scalar:
      switch(field_fst) {
        case FST::Real:
          iterate (this->in_st, field1_rt, field2_rt, this->out_st);
          break;                           
        case FST::MeshScalar:              
          iterate (this->in_st, field1_mt, field2_mt, this->out_st);
          break;                           
        case FST::ParamScalar:             
          iterate (this->in_st, field1_pt, field2_pt, this->out_st);
          break;                           
        case FST::Scalar:                  
          iterate (this->in_st, field1_st, field2_st, this->out_st);
          break;
      } break;
  }
}

//**********************************************************************
template<typename EvalT, typename Traits>
TernaryBandPassOp<EvalT, Traits>::
TernaryBandPassOp (const Teuchos::ParameterList& p,
                   const Teuchos::RCP<Albany::Layouts>& dl) :
  SimpleTernaryOperation<EvalT,Traits,TernaryOps::BandPass>(p,dl) {}

} // Namespace LandIce

