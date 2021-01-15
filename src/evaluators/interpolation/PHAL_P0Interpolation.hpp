//*****************************************************************//
//    Albany 3.0:  Copyright 2016 Sandia Corporation               //
//    This Software is released under the BSD license detailed     //
//    in the file "license.txt" in the top-level Albany directory  //
//*****************************************************************//

#ifndef PHAL_P0_INTERPOLATION_HPP
#define PHAL_P0_INTERPOLATION_HPP

#include "Phalanx_config.hpp"
#include "Phalanx_Evaluator_WithBaseImpl.hpp"
#include "Phalanx_Evaluator_Derived.hpp"
#include "Phalanx_MDField.hpp"

#include "Albany_Layouts.hpp"
#include "PHAL_Utilities.hpp"
#include "Albany_SacadoTypes.hpp"

#include "Intrepid2_CellTools.hpp"

namespace PHAL
{
/** \brief Average from points to cell/side

    This evaluator averages the node/quadpoints values to
    obtain a single value for the whole cell/side

*/

template<typename EvalT, typename Traits>
class P0Interpolation : public PHX::EvaluatorWithBaseImpl<Traits>,
                        public PHX::EvaluatorDerived<EvalT, Traits>
{
public:

  P0Interpolation (const Teuchos::ParameterList& p,
                   const Teuchos::RCP<Albany::Layouts>& dl);

  void postRegistrationSetup (typename Traits::SetupData d,
                              PHX::FieldManager<Traits>& fm);

  void evaluateFields(typename Traits::EvalData d);

private:

  using RT = RealType;
  using MT = typename EvalT::MeshScalarT;
  using PT = typename EvalT::ParamScalarT;
  using ST = typename EvalT::ScalarT;

  template<typename InT>
  using OT = typename Albany::StrongestScalarType<InT,MT>::type;

  void evaluate_on_side (typename Traits::EvalData d);
  void evaluate_on_cell (typename Traits::EvalData d);

  enum InterpolationType {
    ValueAtCellBarycenter,
    CellAverage
  };

  bool eval_on_side;  // Whether the interpolation is on a volume or side field.

  int numQPs; 
  int numNodes;
  int dim0;     // For rank1 and rank2 fields
  int dim1;     // For rank2 fields

  FieldLocation       loc;
  FieldRankType       rank;
  FieldScalarType     fst;
  InterpolationType   itype;

  std::vector<PHX::DataLayout::size_type> dims;

  std::string sideSetName; // Only if eval_on_side=true

  // These are only needed for barycenter interpolation
  Teuchos::RCP<Intrepid2::Basis<PHX::Device, RealType, RealType> > intrepidBasis;
  Kokkos::DynRankView<RealType, PHX::Device>                       basis_at_barycenter;

  MDFieldMemoizer<Traits> memoizer;

  // Input:
  PHX::MDField<const RT>    field_rt;
  PHX::MDField<const MT>    field_mt;
  PHX::MDField<const PT>    field_pt;
  PHX::MDField<const ST>    field_st;

  PHX::MDField<const RT>    BF;
  PHX::MDField<const MT>    w_measure;

  // Output:
  PHX::MDField<OT<RT>>      field_rt_p0;
  PHX::MDField<OT<MT>>      field_mt_p0;
  PHX::MDField<OT<PT>>      field_pt_p0;
  PHX::MDField<OT<ST>>      field_st_p0;

public:

  typedef Kokkos::View<int***, PHX::Device>::execution_space ExecutionSpace;

  struct Cell_Average_Scalar_Field_Tag{};
  struct Cell_Average_Vector_Field_Tag{};
  struct Cell_Average_Tensor_Field_Tag{};
  struct Cell_Barycenter_Scalar_Field_Tag{};
  struct Cell_Barycenter_Vector_Field_Tag{};
  struct Cell_Barycenter_Tensor_Field_Tag{};

  typedef Kokkos::RangePolicy<ExecutionSpace,Cell_Average_Scalar_Field_Tag> Cell_Average_Scalar_Field_Policy;
  typedef Kokkos::RangePolicy<ExecutionSpace,Cell_Average_Vector_Field_Tag> Cell_Average_Vector_Field_Policy;
  typedef Kokkos::RangePolicy<ExecutionSpace,Cell_Average_Tensor_Field_Tag> Cell_Average_Tensor_Field_Policy;
  typedef Kokkos::RangePolicy<ExecutionSpace,Cell_Barycenter_Scalar_Field_Tag> Cell_Barycenter_Scalar_Field_Policy;
  typedef Kokkos::RangePolicy<ExecutionSpace,Cell_Barycenter_Vector_Field_Tag> Cell_Barycenter_Vector_Field_Policy;
  typedef Kokkos::RangePolicy<ExecutionSpace,Cell_Barycenter_Tensor_Field_Tag> Cell_Barycenter_Tensor_Field_Policy;

  KOKKOS_INLINE_FUNCTION
  void operator() (const Cell_Average_Scalar_Field_Tag& tag, const int& cell) const;
  KOKKOS_INLINE_FUNCTION
  void operator() (const Cell_Average_Vector_Field_Tag& tag, const int& cell) const;
  KOKKOS_INLINE_FUNCTION
  void operator() (const Cell_Average_Tensor_Field_Tag& tag, const int& cell) const;
  KOKKOS_INLINE_FUNCTION
  void operator() (const Cell_Barycenter_Scalar_Field_Tag& tag, const int& cell) const;
  KOKKOS_INLINE_FUNCTION
  void operator() (const Cell_Barycenter_Vector_Field_Tag& tag, const int& cell) const;
  KOKKOS_INLINE_FUNCTION
  void operator() (const Cell_Barycenter_Tensor_Field_Tag& tag, const int& cell) const;
};

} // Namespace PHAL

#endif // PHAL_P0_INTERPOLATION_HPP
