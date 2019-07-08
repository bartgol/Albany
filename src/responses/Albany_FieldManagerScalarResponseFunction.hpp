//*****************************************************************//
//    Albany 3.0:  Copyright 2016 Sandia Corporation               //
//    This Software is released under the BSD license detailed     //
//    in the file "license.txt" in the top-level Albany directory  //
//*****************************************************************//

#ifndef ALBANY_FIELD_MANAGER_SCALAR_RESPONSE_FUNCTION_HPP
#define ALBANY_FIELD_MANAGER_SCALAR_RESPONSE_FUNCTION_HPP

#include "Albany_ScalarResponseFunction.hpp"
#include "Albany_StateInfoStruct.hpp" // contains MeshSpecsStuct
#include "PHAL_AlbanyTraits.hpp"
#include "Phalanx_FieldManager.hpp"

namespace Albany {

class AbstractProblem;
class Application;
class StateManager;
struct MeshSpecsStruct;

/*!
 * \brief Reponse function representing the average of the solution values
 */
class FieldManagerScalarResponseFunction : public ScalarResponseFunction {
public:

  //! Constructor
  FieldManagerScalarResponseFunction(
    const Teuchos::RCP<Application>& application,
    const Teuchos::RCP<AbstractProblem>& problem,
    const Teuchos::RCP<MeshSpecsStruct>&  ms,
    const Teuchos::RCP<StateManager>& stateMgr,
    Teuchos::ParameterList& responseParams);

  //! Destructor
  ~FieldManagerScalarResponseFunction() = default;

  std::string name () const { return "FieldManagerScalarResponseFunction"; }

  //! Get the number of responses
  unsigned int numResponses() const { return num_responses; }

  //! Perform post registration setup
  void postRegSetup();

  //! Evaluate responses
  void evaluateResponse(const double current_time,
    const Teuchos::RCP<const Thyra_Vector>& x,
    const Teuchos::RCP<const Thyra_Vector>& xdot,
    const Teuchos::RCP<const Thyra_Vector>& xdotdot,
    const Teuchos::Array<ParamVec>& p,
    const Teuchos::RCP<Thyra_Vector>& g);

  //! Evaluate tangent = dg/dx*dx/dp + dg/dxdot*dxdot/dp + dg/dp
  void evaluateTangent(const double alpha,
    const double beta,
    const double omega,
    const double current_time,
    bool sum_derivs,
    const Teuchos::RCP<const Thyra_Vector>& x,
    const Teuchos::RCP<const Thyra_Vector>& xdot,
    const Teuchos::RCP<const Thyra_Vector>& xdotdot,
    const Teuchos::Array<ParamVec>& p,
    ParamVec* deriv_p,
    const Teuchos::RCP<const Thyra_MultiVector>& Vx,
    const Teuchos::RCP<const Thyra_MultiVector>& Vxdot,
    const Teuchos::RCP<const Thyra_MultiVector>& Vxdotdot,
    const Teuchos::RCP<const Thyra_MultiVector>& Vp,
    const Teuchos::RCP<Thyra_Vector>& g,
    const Teuchos::RCP<Thyra_MultiVector>& gx,
    const Teuchos::RCP<Thyra_MultiVector>& gp);

  void evaluateGradient(const double current_time,
    const Teuchos::RCP<const Thyra_Vector>& x,
    const Teuchos::RCP<const Thyra_Vector>& xdot,
    const Teuchos::RCP<const Thyra_Vector>& xdotdot,
    const Teuchos::Array<ParamVec>& p,
    ParamVec* deriv_p,
    const Teuchos::RCP<Thyra_Vector>& g,
    const Teuchos::RCP<Thyra_MultiVector>& dg_dx,
    const Teuchos::RCP<Thyra_MultiVector>& dg_dxdot,
    const Teuchos::RCP<Thyra_MultiVector>& dg_dxdotdot,
    const Teuchos::RCP<Thyra_MultiVector>& dg_dp);

  //! Evaluate distributed parameter derivative dg/dp, in MultiVector form
  void evaluateDistParamDeriv(
    const double current_time,
    const Teuchos::RCP<const Thyra_Vector>& x,
    const Teuchos::RCP<const Thyra_Vector>& xdot,
    const Teuchos::RCP<const Thyra_Vector>& xdotdot,
    const Teuchos::Array<ParamVec>& param_array,
    const std::string& dist_param_name,
    const Teuchos::RCP<Thyra_MultiVector>& dg_dp);

protected:

  //! Constructor for derived classes
  /*!
   * Derived classes must call setup after using this constructor.
   */
  FieldManagerScalarResponseFunction(
    const Teuchos::RCP<Application>& application,
    const Teuchos::RCP<AbstractProblem>& problem,
    const Teuchos::RCP<MeshSpecsStruct>&  ms,
    const Teuchos::RCP<StateManager>& stateMgr);

  //! Setup method for derived classes
  void setup(Teuchos::ParameterList& responseParams);

  // Do not hide base class setup method
  using ScalarResponseFunction::setup;

  //! Helper function for visualizing response graph 
  template <typename EvalT> 
  void visResponseGraph(const std::string& res_type);

protected:

  //! Application class
  Teuchos::RCP<Application> application;

  //! Problem class
  Teuchos::RCP<AbstractProblem> problem;

  //! Mesh specs
  Teuchos::RCP<MeshSpecsStruct> meshSpecs;

  //! State manager
  Teuchos::RCP<StateManager> stateMgr;

  //! Field manager for Responses 
  Teuchos::RCP<PHX::FieldManager<PHAL::AlbanyTraits> > rfm;

  //! Number of responses we compute
  unsigned int num_responses;

  //! Visualize response graph
  int vis_response_graph;

  //! Response name for visualization file
  std::string vis_response_name;

private:

  template <typename EvalT> void evaluate(PHAL::Workset& workset);

  //! Restrict the field manager to an element block, as is done for fm and
  //! sfm in Application.
  int element_block_index;

  bool performedPostRegSetup;
};

template <typename EvalT> 
void FieldManagerScalarResponseFunction::
visResponseGraph(const std::string& res_type) {
  // Only write out the graph file first time function is called
  static bool first = true;
  if (first && vis_response_graph > 0) {
    bool detail = false; if (vis_response_graph > 1) detail=true;
    Teuchos::RCP<Teuchos::FancyOStream> out = 
      Teuchos::VerboseObjectBase::getDefaultOStream();
    *out << "Phalanx writing graphviz file for graph of Response fill "
         << "(detail = "<< vis_response_graph << ")"<< std::endl;
    std::string detail_name = "responses_graph_" + vis_response_name + res_type;
    *out << "Process using 'dot -Tpng -O ' " << detail_name << "\n" << std::endl;
    rfm->writeGraphvizFile<EvalT>(detail_name,detail,detail);
    first = false;
  }
}

} // namespace Albany

#endif // ALBANY_FIELD_MANAGER_SCALAR_RESPONSE_FUNCTION_HPP
