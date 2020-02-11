//*****************************************************************//
//    Albany 3.0:  Copyright 2016 Sandia Corporation               //
//    This Software is released under the BSD license detailed     //
//    in the file "license.txt" in the top-level Albany directory  //
//*****************************************************************//

#ifndef LANDICE_FLOW_RATE_HPP
#define LANDICE_FLOW_RATE_HPP 1

#include "Phalanx_config.hpp"
#include "Phalanx_Evaluator_WithBaseImpl.hpp"
#include "Phalanx_Evaluator_Derived.hpp"
#include "Phalanx_MDField.hpp"

#include "Albany_Layouts.hpp"
#include "Albany_ScalarOrdinalTypes.hpp"
#include "PHAL_Dimension.hpp"

namespace LandIce
{

/** \brief Hydrology Residual Evaluator

    This evaluator evaluates the residual of the Hydrology model
*/

template<typename EvalT, typename Traits, typename TemperatureST>
class FlowRate : public PHX::EvaluatorWithBaseImpl<Traits>,
                 public PHX::EvaluatorDerived<EvalT, Traits>
{
public:

  FlowRate (const Teuchos::ParameterList& p,
            const Teuchos::RCP<Albany::Layouts>& dl);

  void postRegistrationSetup (typename Traits::SetupData d,
                              PHX::FieldManager<Traits>& fm);

  void evaluateFields(typename Traits::EvalData d);

private:

  typedef typename EvalT::ParamScalarT ParamScalarT;

  // Input:
  PHX::MDField<const RealType,Cell>       given_flow_rate;
  PHX::MDField<const TemperatureST,Cell>  temperature;

  // Output:
  PHX::MDField<TemperatureST,Cell> flowRate;

  double A;
  enum FlowRateType {UNIFORM, GIVEN_FIELD, TEMPERATURE_BASED};
  FlowRateType flowRate_type;
};

} // Namespace LandIce

#endif // LANDICE_FLOW_RATE_HPP
