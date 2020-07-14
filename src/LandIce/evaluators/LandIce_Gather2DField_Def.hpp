//*****************************************************************//
//    Albany 3.0:  Copyright 2016 Sandia Corporation               //
//    This Software is released under the BSD license detailed     //
//    in the file "license.txt" in the top-level Albany directory  //
//*****************************************************************//

#include "Teuchos_TestForException.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Phalanx_DataLayout.hpp"
#include "Phalanx_Print.hpp"
#include "Sacado.hpp"

#include "Albany_AbstractDiscretization.hpp"
#include "Albany_ThyraUtils.hpp"
#include "Albany_GlobalLocalIndexer.hpp"

#include "LandIce_Gather2DField.hpp"

//uncomment the following line if you want debug output to be printed to screen
//#define OUTPUT_TO_SCREEN

namespace LandIce {

//**********************************************************************

template<typename EvalT, typename Traits>
Gather2DFieldBase<EvalT, Traits>::
Gather2DFieldBase(const Teuchos::ParameterList& p,
                  const Teuchos::RCP<Albany::Layouts>& dl)
 : field2D(p.get<std::string>("2D Field Name"), dl->node_scalar)
{
  Teuchos::RCP<Teuchos::FancyOStream> out(Teuchos::VerboseObjectBase::getDefaultOStream());

  this->addEvaluatedField(field2D);
  cell_topo = p.get<Teuchos::RCP<const CellTopologyData> >("Cell Topology");

  std::vector<PHX::DataLayout::size_type> dims;

  dl->node_gradient->dimensions(dims);
  numNodes = dims[1];
  vecDim = dims[2];

  this->setName("Gather2DField"+PHX::print<EvalT>());

  if (p.isType<int>("Offset of First DOF")) {
    offset = p.get<int>("Offset of First DOF");
  } else {
    offset = 2;
  }

  fieldLevel = p.get<int>("Field Level");
}

//**********************************************************************
template<typename EvalT, typename Traits>
void Gather2DFieldBase<EvalT, Traits>::
postRegistrationSetup(typename Traits::SetupData /* d */,
                      PHX::FieldManager<Traits>& fm)
{
    this->utils.setFieldData(field2D,fm);
}

//**********************************************************************

template<typename Traits>
Gather2DField<PHAL::AlbanyTraits::Residual, Traits>::
Gather2DField(const Teuchos::ParameterList& p,
              const Teuchos::RCP<Albany::Layouts>& dl)
 : Gather2DFieldBase<PHAL::AlbanyTraits::Residual, Traits>(p,dl)
{
  if (p.isType<const std::string>("Mesh Part")) {
    this->meshPart = p.get<const std::string>("Mesh Part");
  } else {
    this->meshPart = "upperside";
  }
}

template<typename Traits>
void Gather2DField<PHAL::AlbanyTraits::Residual, Traits>::
evaluateFields(typename Traits::EvalData workset)
{
  auto nodeID = workset.wsElNodeEqID;
  Teuchos::ArrayRCP<const ST> x_constView = Albany::getLocalData(workset.x);

  const Albany::SideSetList& ssList = *(workset.sideSets);
  Albany::SideSetList::const_iterator it = ssList.find(this->meshPart);

  if (it != ssList.end()) {
    const std::vector<Albany::SideStruct>& sideSet = it->second;

    for (std::size_t iSide = 0; iSide < sideSet.size(); ++iSide) { // loop over the sides on this ws and name
      // Get the data that corresponds to the side
      const int elem_LID = sideSet[iSide].elem_LID;
      const int elem_side = sideSet[iSide].side_local_id;
      const CellTopologyData_Subcell& side =  this->cell_topo->side[elem_side];
      int numSideNodes = side.topology->node_count;
      for (int i = 0; i < numSideNodes; ++i){
        std::size_t node = side.node[i];
        this->field2D(elem_LID,node) = x_constView[nodeID(elem_LID,node,this->offset)];
      }
    }
  }
}

template<typename Traits>
Gather2DField<PHAL::AlbanyTraits::Jacobian, Traits>::
Gather2DField(const Teuchos::ParameterList& p,
              const Teuchos::RCP<Albany::Layouts>& dl)
 : Gather2DFieldBase<PHAL::AlbanyTraits::Jacobian, Traits>(p,dl)
{
  if (p.isType<const std::string>("Mesh Part")) {
    this->meshPart = p.get<const std::string>("Mesh Part");
  } else {
    this->meshPart = "upperside";
  }
}

template<typename Traits>
void Gather2DField<PHAL::AlbanyTraits::Jacobian, Traits>::
evaluateFields(typename Traits::EvalData workset)
{
  auto nodeID = workset.wsElNodeEqID;
  Teuchos::ArrayRCP<const ST> x_constView = Albany::getLocalData(workset.x);

  TEUCHOS_TEST_FOR_EXCEPTION(workset.sideSets.is_null(), std::logic_error,
                             "Side sets defined in input file but not properly specified on the mesh.\n");

  const Albany::SideSetList& ssList = *(workset.sideSets);
  Albany::SideSetList::const_iterator it = ssList.find(this->meshPart);

  if (it != ssList.end()) {
    const std::vector<Albany::SideStruct>& sideSet = it->second;

    // Loop over the sides that form the boundary condition
    for (std::size_t iSide = 0; iSide < sideSet.size(); ++iSide) { // loop over the sides on this ws and name

      // Get the data that corresponds to the side
      const int elem_LID = sideSet[iSide].elem_LID;
      const int elem_side = sideSet[iSide].side_local_id;
      const CellTopologyData_Subcell& side =  this->cell_topo->side[elem_side];
      int numSideNodes = side.topology->node_count;

      for (int i = 0; i < numSideNodes; ++i){
        std::size_t node = side.node[i];
        typename PHAL::Ref<ScalarT>::type val = (this->field2D)(elem_LID,node);
        val = FadType(val.size(), x_constView[nodeID(elem_LID,node,this->offset)]);
        val.fastAccessDx(numSideNodes*this->vecDim*this->fieldLevel+this->vecDim*i+this->offset) = workset.j_coeff;
      }
    }
  }
}

template<typename Traits>
Gather2DField<PHAL::AlbanyTraits::Tangent, Traits>::
Gather2DField(const Teuchos::ParameterList& p,
              const Teuchos::RCP<Albany::Layouts>& dl)
 : Gather2DFieldBase<PHAL::AlbanyTraits::Tangent, Traits>(p,dl)
{
  // Nothing to do here
}

template<typename Traits>
Gather2DField<PHAL::AlbanyTraits::DistParamDeriv, Traits>::
Gather2DField(const Teuchos::ParameterList& p,
              const Teuchos::RCP<Albany::Layouts>& dl)
 : Gather2DFieldBase<PHAL::AlbanyTraits::DistParamDeriv, Traits>(p,dl)
{
  // Nothing to do here
}

//********************************

template<typename Traits>
GatherExtruded2DField<PHAL::AlbanyTraits::Residual, Traits>::
GatherExtruded2DField(const Teuchos::ParameterList& p,
                      const Teuchos::RCP<Albany::Layouts>& dl)
 : Gather2DFieldBase<PHAL::AlbanyTraits::Residual, Traits>(p,dl)
{
  this->setName("GatherExtruded2DField Residual");
}

template<typename Traits>
void GatherExtruded2DField<PHAL::AlbanyTraits::Residual, Traits>::
evaluateFields(typename Traits::EvalData workset)
{
  Teuchos::ArrayRCP<const ST> x_constView = Albany::getLocalData(workset.x);

  TEUCHOS_TEST_FOR_EXCEPTION (workset.disc->getLayeredMeshNumbering().is_null(),
    std::runtime_error, "Error! No layered numbering in the mesh.\n");

  const Albany::LayeredMeshNumbering<GO>& layeredMeshNumbering = *workset.disc->getLayeredMeshNumbering();
  const Albany::NodalDOFManager& solDOFManager = workset.disc->getOverlapDOFManager(workset.disc->solution_dof_name());
  const Teuchos::ArrayRCP<Teuchos::ArrayRCP<GO> >& wsElNodeID  = workset.disc->getWsElNodeID()[workset.wsIndex];
  const auto& indexer = *workset.disc->getOverlapNodeGlobalLocalIndexer();

  for (std::size_t cell=0; cell < workset.numCells; ++cell ) {
    const Teuchos::ArrayRCP<GO>& elNodeID = wsElNodeID[cell];

    for (std::size_t node = 0; node < this->numNodes; ++node) {
      // Retrieve corresponding 2D node
      const GO base_id = layeredMeshNumbering.getColumnId(elNodeID[node]);
      GO gnode = layeredMeshNumbering.getId(base_id, this->fieldLevel);
      LO lnode = indexer.getLocalElement(gnode);
      (this->field2D)(cell,node) = x_constView[solDOFManager.getLocalDOF(lnode, this->offset)];
    }
  }
}

template<typename Traits>
GatherExtruded2DField<PHAL::AlbanyTraits::Jacobian, Traits>::
GatherExtruded2DField(const Teuchos::ParameterList& p,
                      const Teuchos::RCP<Albany::Layouts>& dl)
 : Gather2DFieldBase<PHAL::AlbanyTraits::Jacobian, Traits>(p,dl)
{
  this->setName("GatherExtruded2DField Jacobian");
}

template<typename Traits>
void GatherExtruded2DField<PHAL::AlbanyTraits::Jacobian, Traits>::
evaluateFields(typename Traits::EvalData workset)
{
  auto nodeID = workset.wsElNodeEqID;
  Teuchos::ArrayRCP<const ST> x_constView = Albany::getLocalData(workset.x);

  TEUCHOS_TEST_FOR_EXCEPTION (workset.disc->getLayeredMeshNumbering().is_null(),
    std::runtime_error, "Error! No layered numbering in the mesh.\n");

  const Albany::LayeredMeshNumbering<GO>& layeredMeshNumbering = *workset.disc->getLayeredMeshNumbering();
  const Albany::NodalDOFManager& solDOFManager = workset.disc->getOverlapDOFManager("ordinary_solution");
  const auto& indexer = *workset.disc->getOverlapGlobalLocalIndexer();

  int numLayers = layeredMeshNumbering.numLayers;
  this->fieldLevel = (this->fieldLevel < 0) ? numLayers : this->fieldLevel;
  const Teuchos::ArrayRCP<Teuchos::ArrayRCP<GO> >& wsElNodeID  = workset.disc->getWsElNodeID()[workset.wsIndex];

  for (std::size_t cell=0; cell < workset.numCells; ++cell ) {
    const Teuchos::ArrayRCP<GO>& elNodeID = wsElNodeID[cell];
    const int neq = nodeID.extent(2);

    for (std::size_t node = 0; node < this->numNodes; ++node) {
      int firstunk = neq * node + this->offset;
      const GO base_id = layeredMeshNumbering.getColumnId(elNodeID[node]);
      GO gnode = layeredMeshNumbering.getId(base_id, this->fieldLevel);
      GO gdof = solDOFManager.getGlobalDOF(gnode, this->offset);
      typename PHAL::Ref<ScalarT>::type val = (this->field2D)(cell,node);

      LO ldof = indexer.getLocalElement(gdof);
      val = FadType(val.size(), x_constView[ldof]);
      val.setUpdateValue(!workset.ignore_residual);
      val.fastAccessDx(firstunk) = workset.j_coeff;
    }
  }
}

template<typename Traits>
GatherExtruded2DField<PHAL::AlbanyTraits::Tangent, Traits>::
GatherExtruded2DField(const Teuchos::ParameterList& p,
                      const Teuchos::RCP<Albany::Layouts>& dl)
 : Gather2DFieldBase<PHAL::AlbanyTraits::Tangent, Traits>(p,dl)
{
  this->setName("GatherExtruded2DField Tangent");
}

template<typename Traits>
GatherExtruded2DField<PHAL::AlbanyTraits::DistParamDeriv, Traits>::
GatherExtruded2DField(const Teuchos::ParameterList& p,
                      const Teuchos::RCP<Albany::Layouts>& dl)
 : Gather2DFieldBase<PHAL::AlbanyTraits::DistParamDeriv, Traits>(p,dl)
{
  this->setName("GatherExtruded2DField DistParamDeriv");
}

} // namespace LandIce
