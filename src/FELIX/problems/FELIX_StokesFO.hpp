//*****************************************************************//
//    Albany 2.0:  Copyright 2012 Sandia Corporation               //
//    This Software is released under the BSD license detailed     //
//    in the file "license.txt" in the top-level Albany directory  //
//*****************************************************************//

#ifndef FELIX_STOKESFOPROBLEM_HPP
#define FELIX_STOKESFOPROBLEM_HPP

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "Albany_AbstractProblem.hpp"

#include "Phalanx.hpp"
#include "PHAL_Workset.hpp"
#include "PHAL_Dimension.hpp"
//#include "PHAL_AlbanyTraits.hpp"
#include "PHAL_SaveStateField.hpp"
#include "PHAL_LoadStateField.hpp"

//uncomment the following line if you want debug output to be printed to screen
//#define OUTPUT_TO_SCREEN

namespace FELIX
{

template<typename EvalT,typename Traits>
class HomotopyParamValue
{
public:
    static typename EvalT::ScalarT* value;
};

template<typename EvalT,typename Traits>
typename EvalT::ScalarT* HomotopyParamValue<EvalT,Traits>::value = NULL;
/*!
 * \brief Abstract interface for representing a 1-D finite element
 * problem.
 */
class StokesFO : public Albany::AbstractProblem
{
public:

  //! Default constructor
  StokesFO (const Teuchos::RCP<Teuchos::ParameterList>& params,
            const Teuchos::RCP<ParamLib>& paramLib,
            const int numDim_);

  //! Destructor
  ~StokesFO();

  //! Return number of spatial dimensions
  virtual int spatialDimension() const { return numDim; }

  //! Build the PDE instantiations, boundary conditions, and initial solution
  virtual void buildProblem (Teuchos::ArrayRCP<Teuchos::RCP<Albany::MeshSpecsStruct> >  meshSpecs,
                             Albany::StateManager& stateMgr);

  // Build evaluators
  virtual Teuchos::Array< Teuchos::RCP<const PHX::FieldTag> >
  buildEvaluators (PHX::FieldManager<PHAL::AlbanyTraits>& fm0,
                   const Albany::MeshSpecsStruct& meshSpecs,
                   Albany::StateManager& stateMgr,
                   Albany::FieldManagerChoice fmchoice,
                   const Teuchos::RCP<Teuchos::ParameterList>& responseList);

  //! Each problem must generate it's list of valide parameters
  Teuchos::RCP<const Teuchos::ParameterList> getValidProblemParameters() const;

private:

  //! Private to prohibit copying
  StokesFO(const StokesFO&);

  //! Private to prohibit copying
  StokesFO& operator=(const StokesFO&);

public:

  //! Main problem setup routine. Not directly called, but indirectly by following functions
  template <typename EvalT> Teuchos::RCP<const PHX::FieldTag>
  constructEvaluators (PHX::FieldManager<PHAL::AlbanyTraits>& fm0,
                       const Albany::MeshSpecsStruct& meshSpecs,
                       Albany::StateManager& stateMgr,
                       Albany::FieldManagerChoice fmchoice,
                       const Teuchos::RCP<Teuchos::ParameterList>& responseList);

  void constructDirichletEvaluators(const Albany::MeshSpecsStruct& meshSpecs);
  void constructNeumannEvaluators(const Teuchos::RCP<Albany::MeshSpecsStruct>& meshSpecs);

protected:

  // Used to build basal friction evaluator for all evaluation types
  struct ConstructBasalEvaluatorOp
  {
      StokesFO& prob_;
      std::vector<Teuchos::RCP<PHX::Evaluator<PHAL::AlbanyTraits> > >& evaluators_;

      ConstructBasalEvaluatorOp (StokesFO& prob,
                                 std::vector<Teuchos::RCP<PHX::Evaluator<PHAL::AlbanyTraits> > >& evaluators) :
          prob_(prob), evaluators_(evaluators) {}
      template<typename T>
      void operator() (T x) {
      evaluators_.push_back(prob_.template buildBasalFrictionCoefficientEvaluator<T>());
      evaluators_.push_back(prob_.template buildSlidingVelocityEvaluator<T>());
      }
  };

  template<typename EvalT>
  Teuchos::RCP<PHX::Evaluator<PHAL::AlbanyTraits> >
  buildBasalFrictionCoefficientEvaluator ();

  template<typename EvalT>
  Teuchos::RCP<PHX::Evaluator<PHAL::AlbanyTraits> >
  buildSlidingVelocityEvaluator();

  int numDim;
  double gravity;  //gravity
  double rho;  //ice density
  double rho_w;  //water density
  Teuchos::RCP<Albany::Layouts> dl;
};

} // Namespace FELIX

#include "Intrepid_FieldContainer.hpp"
#include "Intrepid_DefaultCubatureFactory.hpp"
#include "Shards_CellTopology.hpp"

#include "Albany_Utils.hpp"
#include "Albany_ProblemUtils.hpp"
#include "Albany_EvaluatorUtils.hpp"
#include "Albany_ResponseUtilities.hpp"

#include "FELIX_StokesFOResid.hpp"
#ifdef CISM_HAS_FELIX
#include "FELIX_CismSurfaceGradFO.hpp"
#endif
#include "FELIX_StokesFOBodyForce.hpp"
#include "FELIX_ViscosityFO.hpp"
#include "FELIX_FieldNorm.hpp"
#include "FELIX_SaveSideSetStateField.hpp"
#include "FELIX_QuadPointsToCellInterpolation.hpp"
#include "FELIX_BasalFrictionCoefficient.hpp"
#include "PHAL_Neumann.hpp"
#include "PHAL_Source.hpp"
#include <type_traits>

template <typename EvalT>
Teuchos::RCP<const PHX::FieldTag>
FELIX::StokesFO::constructEvaluators (PHX::FieldManager<PHAL::AlbanyTraits>& fm0,
                                      const Albany::MeshSpecsStruct& meshSpecs,
                                      Albany::StateManager& stateMgr,
                                      Albany::FieldManagerChoice fieldManagerChoice,
                                      const Teuchos::RCP<Teuchos::ParameterList>& responseList)
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::ParameterList;
  using PHX::DataLayout;
  using PHX::MDALayout;
  using std::vector;
  using std::string;
  using std::map;
  using PHAL::AlbanyTraits;

  RCP<Intrepid::Basis<RealType, Intrepid::FieldContainer<RealType> > >
    intrepidBasis = Albany::getIntrepidBasis(meshSpecs.ctd);
  RCP<shards::CellTopology> cellType = rcp(new shards::CellTopology (&meshSpecs.ctd));

  const int numNodes = intrepidBasis->getCardinality();
  const int worksetSize = meshSpecs.worksetSize;

  Intrepid::DefaultCubatureFactory<RealType> cubFactory;
  RCP <Intrepid::Cubature<RealType> > cubature = cubFactory.create(*cellType, meshSpecs.cubatureDegree);

  const int numQPts = cubature->getNumPoints();
  const int numVertices = cellType->getNodeCount();
  int vecDim = neq;
  std::string elementBlockName = meshSpecs.ebName;

#ifdef OUTPUT_TO_SCREEN
  *out << "Field Dimensions: Workset=" << worksetSize
       << ", Vertices= " << numVertices
       << ", Nodes= " << numNodes
       << ", QuadPts= " << numQPts
       << ", Dim= " << numDim
       << ", vecDim= " << vecDim << std::endl;
#endif

   Albany::StateStruct::MeshFieldEntity entity;
   dl = rcp(new Albany::Layouts(worksetSize,numVertices,numNodes,numQPts,numDim, vecDim));
   Albany::EvaluatorUtils<EvalT, PHAL::AlbanyTraits> evalUtils(dl);
   int offset=0;

   entity= Albany::StateStruct::ElemData;

   // Temporary variable used numerous times below
      RCP<PHX::Evaluator<AlbanyTraits> > ev;
   {
     std::string stateName("temperature");
     RCP<ParameterList> p = stateMgr.registerStateVariable(stateName, dl->cell_scalar2, elementBlockName,true, &entity);
     ev = rcp(new PHAL::LoadStateField<EvalT,AlbanyTraits>(*p));
     fm0.template registerEvaluator<EvalT>(ev);
   }

   {
     std::string stateName("flow_factor");
     RCP<ParameterList> p = stateMgr.registerStateVariable(stateName, dl->cell_scalar2, elementBlockName,true, &entity);
     ev = rcp(new PHAL::LoadStateField<EvalT,AlbanyTraits>(*p));
     fm0.template registerEvaluator<EvalT>(ev);
   }

   entity= Albany::StateStruct::NodalDataToElemNode;

   {
     std::string stateName("surface_height");
     RCP<ParameterList> p = stateMgr.registerStateVariable(stateName, dl->node_scalar, elementBlockName,true, &entity);
     ev = rcp(new PHAL::LoadStateField<EvalT,AlbanyTraits>(*p));
     fm0.template registerEvaluator<EvalT>(ev);
   }
#ifdef CISM_HAS_FELIX
   {
     std::string stateName("xgrad_surface_height"); //ds/dx which can be passed from CISM (defined at nodes)
     RCP<ParameterList> p = stateMgr.registerStateVariable(stateName, dl->node_scalar, elementBlockName,true, &entity);
     ev = rcp(new PHAL::LoadStateField<EvalT,AlbanyTraits>(*p));
     fm0.template registerEvaluator<EvalT>(ev);
   }
   {
     std::string stateName("ygrad_surface_height"); //ds/dy which can be passed from CISM (defined at nodes)
     RCP<ParameterList> p = stateMgr.registerStateVariable(stateName, dl->node_scalar, elementBlockName,true, &entity);
     ev = rcp(new PHAL::LoadStateField<EvalT,AlbanyTraits>(*p));
     fm0.template registerEvaluator<EvalT>(ev);
   }
#endif
   {
     std::string stateName("thickness");
     RCP<ParameterList> p = stateMgr.registerStateVariable(stateName, dl->node_scalar, elementBlockName,true, &entity);
     ev = rcp(new PHAL::LoadStateField<EvalT,AlbanyTraits>(*p));
     fm0.template registerEvaluator<EvalT>(ev);
   }

   entity= Albany::StateStruct::NodalDistParameter;

   {
     std::string stateName("basal_friction");
     const std::string& meshPart = this->params->sublist("Distributed Parameters").get("Mesh Part","");
     RCP<ParameterList> p = stateMgr.registerStateVariable(stateName, dl->node_scalar, elementBlockName,true, &entity, meshPart);
     fm0.template registerEvaluator<EvalT>
         (evalUtils.constructGatherScalarNodalParameter(stateName));
    }

#if defined(CISM_HAS_FELIX) || defined(MPAS_HAS_FELIX)
   {
    // Here is how to register the field for dirichlet condition.
    std::string stateName("dirichlet_field");
    // IK, 12/9/14: Changed "false" to "true" from Mauro's initial implementation for outputting to Exodus
    RCP<ParameterList> p = stateMgr.registerStateVariable(stateName, dl->node_vector, elementBlockName, true, &entity);
     ev = rcp(new PHAL::LoadStateField<EvalT,AlbanyTraits>(*p));
     fm0.template registerEvaluator<EvalT>(ev);
   }
#endif


   // Define Field Names

  Teuchos::ArrayRCP<std::string> dof_names(1);
  Teuchos::ArrayRCP<std::string> dof_names_dot(1);
  Teuchos::ArrayRCP<std::string> resid_names(1);
  dof_names[0] = "Velocity";
  //dof_names_dot[0] = dof_names[0]+"_dot";
  resid_names[0] = "Stokes Residual";
  fm0.template registerEvaluator<EvalT>
  (evalUtils.constructGatherSolutionEvaluator_noTransient(true, dof_names, offset));

  fm0.template registerEvaluator<EvalT>
    (evalUtils.constructDOFVecInterpolationEvaluator(dof_names[0]));

  fm0.template registerEvaluator<EvalT>
    (evalUtils.constructDOFVecGradInterpolationEvaluator(dof_names[0]));

  fm0.template registerEvaluator<EvalT>
    (evalUtils.constructScatterResidualEvaluator(true, resid_names,offset, "Scatter Stokes"));
  offset += numDim;

  fm0.template registerEvaluator<EvalT>
    (evalUtils.constructGatherCoordinateVectorEvaluator());

  fm0.template registerEvaluator<EvalT>
    (evalUtils.constructMapToPhysicalFrameEvaluator(cellType, cubature));

  fm0.template registerEvaluator<EvalT>
    (evalUtils.constructComputeBasisFunctionsEvaluator(cellType, intrepidBasis, cubature));

  std::string sh = "surface_height";
  fm0.template registerEvaluator<EvalT>
    (evalUtils.constructDOFGradInterpolationEvaluator_noDeriv(sh));


  { // FO Stokes Resid
    RCP<ParameterList> p = rcp(new ParameterList("Stokes Resid"));

    //Input
    p->set<std::string>("Weighted BF Name", "wBF");
    p->set<std::string>("Weighted Gradient BF Name", "wGrad BF");
    p->set<std::string>("QP Variable Name", "Velocity");
    p->set<std::string>("QP Time Derivative Variable Name", "Velocity_dot");
    p->set<std::string>("Gradient QP Variable Name", "Velocity Gradient");
    p->set<std::string>("Body Force Name", "Body Force");
    p->set<std::string>("FELIX Viscosity QP Variable Name", "FELIX Viscosity");

    Teuchos::ParameterList& paramList = params->sublist("Equation Set");
    p->set<Teuchos::ParameterList*>("Parameter List", &paramList);

    //Output
    p->set<std::string>("Residual Name", "Stokes Residual");

    ev = rcp(new FELIX::StokesFOResid<EvalT,AlbanyTraits>(*p,dl));
    fm0.template registerEvaluator<EvalT>(ev);
  }

  { // FELIX viscosity
    RCP<ParameterList> p = rcp(new ParameterList("FELIX Viscosity"));

    //Input
    p->set<std::string>("Coordinate Vector Name", "Coord Vec");
    p->set<std::string>("Gradient QP Variable Name", "Velocity Gradient");
    p->set<std::string>("temperature Name", "temperature");
    p->set<std::string>("flow_factor Name", "flow_factor");

    p->set<RCP<ParamLib> >("Parameter Library", paramLib);
    Teuchos::ParameterList& paramList = params->sublist("FELIX Viscosity");
    p->set<Teuchos::ParameterList*>("Parameter List", &paramList);

    //Output
    p->set<std::string>("FELIX Viscosity QP Variable Name", "FELIX Viscosity");

    ev = rcp(new FELIX::ViscosityFO<EvalT,AlbanyTraits>(*p,dl));

    typename EvalT::ScalarT** value = &HomotopyParamValue<EvalT,PHAL::AlbanyTraits>::value;
    if (*value==NULL)
    {
        typedef typename Sacado::ParameterAccessor<EvalT, SPL_Traits> sacado_accessor_type;
        sacado_accessor_type* pa_ptr;
        pa_ptr = dynamic_cast<sacado_accessor_type*>(&(*ev));
        if (pa_ptr==0)
        {
            std::cout << "Error! Cannot cast the pointer...\n";
            std::abort();
        }
        *value = &pa_ptr->getValue("Glen's Law Homotopy Parameter");
    }
    fm0.template registerEvaluator<EvalT>(ev);
  }

  // Sliding velocity calculation
  {
    ev = buildSlidingVelocityEvaluator<EvalT>();
    fm0.template registerEvaluator<EvalT>(ev);
  }

  { // FELIX basal friction coefficient
    ev = buildBasalFrictionCoefficientEvaluator<EvalT>();
    fm0.template registerEvaluator<EvalT>(ev);
  }
#ifdef CISM_HAS_FELIX
  { // FELIX surface gradient from CISM
    RCP<ParameterList> p = rcp(new ParameterList("FELIX Surface Gradient"));

    //Input
    p->set<std::string>("xgrad_surface_height Name", "xgrad_surface_height");
    p->set<std::string>("ygrad_surface_height Name", "ygrad_surface_height");
    p->set<std::string>("BF Name", "BF");

    p->set<RCP<ParamLib> >("Parameter Library", paramLib);
    Teuchos::ParameterList& paramList = params->sublist("FELIX Surface Gradient");
    p->set<Teuchos::ParameterList*>("Parameter List", &paramList);

    //Output
    p->set<std::string>("FELIX Surface Gradient QP Name", "FELIX Surface Gradient");

    ev = rcp(new FELIX::CismSurfaceGradFO<EvalT,AlbanyTraits>(*p,dl));
    fm0.template registerEvaluator<EvalT>(ev);

  }
#endif

  { // Body Force
    RCP<ParameterList> p = rcp(new ParameterList("Body Force"));

    //Input
    p->set<std::string>("FELIX Viscosity QP Variable Name", "FELIX Viscosity");
#ifdef CISM_HAS_FELIX
    p->set<std::string>("FELIX Surface Gradient QP Variable Name", "FELIX Surface Gradient");
#endif
    p->set<std::string>("Coordinate Vector Name", "Coord Vec");
    p->set<std::string>("surface_height Gradient Name", "surface_height Gradient");

    Teuchos::ParameterList& paramList = params->sublist("Body Force");
    p->set<Teuchos::ParameterList*>("Parameter List", &paramList);

    Teuchos::ParameterList& physParamList = params->sublist("Physical Parameters");
    p->set<Teuchos::ParameterList*>("Physical Parameter List", &physParamList);

    //Output
    p->set<std::string>("Body Force Name", "Body Force");

    ev = rcp(new FELIX::StokesFOBodyForce<EvalT,AlbanyTraits>(*p,dl));
    fm0.template registerEvaluator<EvalT>(ev);
  }

  RCP<ParameterList> paramList = rcp(new ParameterList("Param List"));
  { // response
    RCP<const Albany::MeshSpecsStruct> meshSpecsPtr = Teuchos::rcpFromRef(meshSpecs);
    paramList->set<RCP<const Albany::MeshSpecsStruct> >("Mesh Specs Struct", meshSpecsPtr);
    paramList->set<RCP<ParamLib> >("Parameter Library", paramLib);
  }
  if (params->get("Ice-Hydrology Coupling",true) && std::is_same<EvalT,AlbanyTraits::Residual>::value)
  {
    // If we are solving the StokesFO problem, in the big picture of the coupling
    // with the hydrology problem, we need to save some data on the 2D mesh

    // Save Friction Coefficient
    {
      RCP<ParameterList> p = stateMgr.registerSideSetStateVariable ("basalside","beta_field","basal_friction",
                                                                    dl->node_scalar, elementBlockName,true);
      p->set<Teuchos::RCP<PHX::DataLayout> >("Dummy Data Layout",dl->dummy);

      ev = rcp(new FELIX::SaveSideSetStateField<EvalT,PHAL::AlbanyTraits>(*p,meshSpecs));
      fm0.template registerEvaluator<EvalT>(ev);
      PHX::Tag<typename EvalT::ScalarT> tag("beta_field", dl->dummy);
      fm0.requireField<EvalT>(tag);
    }

    // Ice viscosity interpolation from QP to Cell
    {
      RCP<ParameterList> p = rcp(new ParameterList("Quad Points To Cell Interpolation"));
      p->set<std::string>("Field QP Name","FELIX Viscosity");
      p->set<std::string>("Field Cell Name","FELIX Viscosity Cell");

      ev = rcp(new FELIX::QuadPointsToCellInterpolation<EvalT,PHAL::AlbanyTraits>(*p,dl));
      fm0.template registerEvaluator<EvalT>(ev);
    }

    // Save Ice Viscosity
    {
      RCP<ParameterList> p = stateMgr.registerSideSetStateVariable ("basalside","FELIX Viscosity Cell","ice_viscosity",
                                                                    dl->cell_scalar2, elementBlockName,true);
      p->set<Teuchos::RCP<PHX::DataLayout> >("Dummy Data Layout",dl->dummy);

      ev = rcp(new FELIX::SaveSideSetStateField<EvalT,PHAL::AlbanyTraits>(*p,meshSpecs));
      fm0.template registerEvaluator<EvalT>(ev);
      PHX::Tag<typename EvalT::ScalarT> tag("FELIX Viscosity Cell", dl->dummy);
      fm0.requireField<EvalT>(tag);
    }

    // Save Sliding Velocity
    {
      RCP<ParameterList> p = stateMgr.registerSideSetStateVariable ("basalside","Velocity Norm","sliding_velocity",
                                                                    dl->node_scalar, elementBlockName,true);
      p->set<Teuchos::RCP<PHX::DataLayout> >("Dummy Data Layout",dl->dummy);

      ev = rcp(new FELIX::SaveSideSetStateField<EvalT,PHAL::AlbanyTraits>(*p,meshSpecs));
      fm0.template registerEvaluator<EvalT>(ev);
      PHX::Tag<typename EvalT::ScalarT> tag("Velocity Norm", dl->dummy);
      fm0.requireField<EvalT>(tag);
    }
  }

  if (fieldManagerChoice == Albany::BUILD_RESID_FM)  {
    PHX::Tag<typename EvalT::ScalarT> res_tag("Scatter Stokes", dl->dummy);
    fm0.requireField<EvalT>(res_tag);
  }
  else if (fieldManagerChoice == Albany::BUILD_RESPONSE_FM)
  {
    entity= Albany::StateStruct::NodalDataToElemNode;

    {
      std::string stateName("surface_velocity");
      RCP<ParameterList> p = stateMgr.registerStateVariable(stateName, dl->node_vector, elementBlockName,true,&entity);
      ev = rcp(new PHAL::LoadStateField<EvalT,AlbanyTraits>(*p));
      fm0.template registerEvaluator<EvalT>(ev);
    }

    {
      std::string stateName("surface_velocity_rms");
      RCP<ParameterList> p = stateMgr.registerStateVariable(stateName, dl->node_vector, elementBlockName,true,&entity);
      ev = rcp(new PHAL::LoadStateField<EvalT,AlbanyTraits>(*p));
      fm0.template registerEvaluator<EvalT>(ev);
    }

    Albany::ResponseUtilities<EvalT, PHAL::AlbanyTraits> respUtils(dl);
    return respUtils.constructResponses(fm0, *responseList, paramList, stateMgr);
  }

  return Teuchos::null;
}

template<typename EvalT>
Teuchos::RCP<PHX::Evaluator<PHAL::AlbanyTraits> >
FELIX::StokesFO::buildSlidingVelocityEvaluator ()
{
  Teuchos::RCP<Teuchos::ParameterList> p = Teuchos::rcp(new Teuchos::ParameterList("FELIX Velocity Norm"));
  p->set<std::string>("Field Name","Velocity");
  p->set<std::string>("Field Norm Name","Velocity Norm");

  // Need a more specific pointer to access the setHomotopyParamPtr method
  Teuchos::RCP<FELIX::FieldNorm<EvalT,PHAL::AlbanyTraits> > ev;
  ev = Teuchos::rcp(new FELIX::FieldNorm<EvalT,PHAL::AlbanyTraits>(*p,dl));
  ev->setHomotopyParamPtr(HomotopyParamValue<EvalT,PHAL::AlbanyTraits>::value);

  return ev;
}

template<typename EvalT>
Teuchos::RCP<PHX::Evaluator<PHAL::AlbanyTraits> >
FELIX::StokesFO::buildBasalFrictionCoefficientEvaluator ()
{
  Teuchos::RCP<Teuchos::ParameterList> p = Teuchos::rcp(new Teuchos::ParameterList("FELIX Basal Friction Coefficient"));

  //Input fields
  p->set<std::string>("Velocity Norm Name", "Velocity Norm");
  p->set<std::string>("Given Beta Field Name", "basal_friction");
  p->set<std::string>("thickness Field Name", "thickness");

  //Input physics parameters
  Teuchos::ParameterList& physics = this->params->sublist("Physical Parameters");

  Teuchos::ParameterList& paramList = this->params->sublist("FELIX Basal Friction Coefficient");
  p->set<Teuchos::ParameterList*>("Parameter List", &paramList);
  p->set<Teuchos::ParameterList*>("Physical Parameters", &physics);

  //Output
  p->set<std::string>("FELIX Basal Friction Coefficient Name", "beta_field");

  Teuchos::RCP<FELIX::BasalFrictionCoefficient<EvalT,PHAL::AlbanyTraits> > ev;
  ev = Teuchos::rcp(new FELIX::BasalFrictionCoefficient<EvalT,PHAL::AlbanyTraits>(*p,dl));
  ev->setHomotopyParamPtr(HomotopyParamValue<EvalT,PHAL::AlbanyTraits>::value);

  return ev;
}

#endif // FELIX_STOKESFOPROBLEM_HPP
