#include "Albany_FESpace.hpp"

#include "Albany_ThyraUtils.hpp"

// 1D
#include <Intrepid2_HGRAD_LINE_C1_FEM.hpp>
#include <Intrepid2_HGRAD_LINE_Cn_FEM.hpp>
#include <Intrepid2_HVOL_LINE_Cn_FEM.hpp>

// 2D
#include <Intrepid2_HGRAD_TRI_C1_FEM.hpp>
#include <Intrepid2_HGRAD_TRI_C2_FEM.hpp>
#include <Intrepid2_HGRAD_TRI_Cn_FEM.hpp>

#include <Intrepid2_HGRAD_QUAD_C1_FEM.hpp>
#include <Intrepid2_HGRAD_QUAD_C2_FEM.hpp>
#include <Intrepid2_HGRAD_QUAD_Cn_FEM.hpp>

// 3D
#include <Intrepid2_HGRAD_TET_C1_FEM.hpp>
#include <Intrepid2_HGRAD_TET_C2_FEM.hpp>
#include <Intrepid2_HGRAD_TET_Cn_FEM.hpp>

#include <Intrepid2_HGRAD_HEX_C1_FEM.hpp>
#include <Intrepid2_HGRAD_HEX_C2_FEM.hpp>
#include <Intrepid2_HGRAD_HEX_Cn_FEM.hpp>

#include <Intrepid2_HGRAD_WEDGE_C1_FEM.hpp>
#include <Intrepid2_HGRAD_WEDGE_C2_FEM.hpp>

#include <set>

namespace Albany {

FESpace::
FESpace (Teuchos::RCP<const MeshConnectivity> connectivity,
         const EFunctionSpace                 type,
         const int                            degree,
         const EPointType                     pointType)
 : m_connectivity (connectivity)
 , m_type (type)
{
  TEUCHOS_TEST_FOR_EXCEPTION (
      degree<=0 || degree>10, std::runtime_error,
      "Error! Unsupported polynomial degree.\n");
  TEUCHOS_TEST_FOR_EXCEPTION (
      connectivity.is_null(), std::runtime_error,
      "Error! Null connectivity pointer.\n");
  createBasis (degree,pointType);
  initialize ();
}

void FESpace::createBasis (const int degree, const EPointType pointType)
{
  // This makes the code shorter, and it's only local to this method
  using namespace Intrepid2;
  using Device = PHX::Device;

  const auto& ctd = m_connectivity->getCellTopology()->getBaseCellTopologyData();
  const std::string& name = ctd->name;
  const std::string cell_base_name = name.substr(0,name.find('_'));
  const std::string fe_space_type_str = EFunctionSpaceToString(m_type);

  if (cell_base_name=="Line") {
    switch (m_type) {
      case FUNCTION_SPACE_HGRAD:
        m_fe_basis = (degree==1)
                   ? Teuchos::rcp<BasisType>(new Basis_HGRAD_LINE_C1_FEM<Device>())
                   : Teuchos::rcp<BasisType>(new Basis_HGRAD_LINE_Cn_FEM<Device>(degree, pointType));
        break;
      case FUNCTION_SPACE_HVOL:
        m_fe_basis = Teuchos::rcp<BasisType>(new Basis_HVOL_LINE_Cn_FEM<Device>(degree, pointType));
        break;
      default:
        TEUCHOS_TEST_FOR_EXCEPTION (
          true, std::logic_error,
          "Error! Unsupported FE space type '" + fe_space_type_str + "' on cell '" + cell_base_name + "'.\n");
    }
  } else if (cell_base_name=="Triangle") {
    switch (m_type) {
      case FUNCTION_SPACE_HGRAD:
        m_fe_basis = (degree==1)
                   ? Teuchos::rcp<BasisType>(new Basis_HGRAD_TRI_C1_FEM<Device>())
                   : (degree==2
                      ? Teuchos::rcp<BasisType>(new Basis_HGRAD_TRI_C2_FEM<Device>())
                      : Teuchos::rcp<BasisType>(new Basis_HGRAD_TRI_Cn_FEM<Device>(degree, pointType)));
        break;
      case FUNCTION_SPACE_HCURL:
        break;
      case FUNCTION_SPACE_HDIV:
        break;
      case FUNCTION_SPACE_HVOL:
        break;
      default:
        TEUCHOS_TEST_FOR_EXCEPTION (
          true, std::logic_error,
          "Error! Unsupported FE space type '" + fe_space_type_str + "' on cell '" + cell_base_name + "'.\n");
    }
  } else if (cell_base_name=="Quadrilateral") {
    switch (m_type) {
      case FUNCTION_SPACE_HGRAD:
        m_fe_basis = (degree==1)
                   ? Teuchos::rcp<BasisType>(new Basis_HGRAD_QUAD_C1_FEM<Device>())
                   : (degree==2
                      ? Teuchos::rcp<BasisType>(new Basis_HGRAD_QUAD_C2_FEM<Device>())
                      : Teuchos::rcp<BasisType>(new Basis_HGRAD_QUAD_Cn_FEM<Device>(degree, pointType)));
        break;
      case FUNCTION_SPACE_HCURL:
        break;
      case FUNCTION_SPACE_HDIV:
        break;
      case FUNCTION_SPACE_HVOL:
        break;
      default:
        TEUCHOS_TEST_FOR_EXCEPTION (
          true, std::logic_error,
          "Error! Unsupported FE space type '" + fe_space_type_str + "' on cell '" + cell_base_name + "'.\n");
    }
  } else if (cell_base_name=="Tetrahedron") {
    switch (m_type) {
      case FUNCTION_SPACE_HGRAD:
        m_fe_basis = (degree==1)
                   ? Teuchos::rcp<BasisType>(new Basis_HGRAD_TET_C1_FEM<Device>())
                   : (degree==2
                      ? Teuchos::rcp<BasisType>(new Basis_HGRAD_TET_C2_FEM<Device>())
                      : Teuchos::rcp<BasisType>(new Basis_HGRAD_TET_Cn_FEM<Device>(degree, pointType)));
        break;
      case FUNCTION_SPACE_HCURL:
        break;
      case FUNCTION_SPACE_HDIV:
        break;
      case FUNCTION_SPACE_HVOL:
        break;
      default:
        TEUCHOS_TEST_FOR_EXCEPTION (
          true, std::logic_error,
          "Error! Unsupported FE space type '" + fe_space_type_str + "' on cell '" + cell_base_name + "'.\n");
    }
  } else if (cell_base_name=="Hexahedron") {
    switch (m_type) {
      case FUNCTION_SPACE_HGRAD:
        m_fe_basis = (degree==1)
                   ? Teuchos::rcp<BasisType>(new Basis_HGRAD_HEX_C1_FEM<Device>())
                   : (degree==2
                      ? Teuchos::rcp<BasisType>(new Basis_HGRAD_HEX_C2_FEM<Device>())
                      : Teuchos::rcp<BasisType>(new Basis_HGRAD_HEX_Cn_FEM<Device>(degree, pointType)));
        break;
      case FUNCTION_SPACE_HCURL:
        break;
      case FUNCTION_SPACE_HDIV:
        break;
      case FUNCTION_SPACE_HVOL:
        break;
      default:
        TEUCHOS_TEST_FOR_EXCEPTION (
          true, std::logic_error,
          "Error! Unsupported FE space type '" + fe_space_type_str + "' on cell '" + cell_base_name + "'.\n");
    }
  } else if (cell_base_name=="Wedge") {
    switch (m_type) {
      case FUNCTION_SPACE_HGRAD:
        TEUCHOS_TEST_FOR_EXCEPTION(
            degree!=1 && degree!=2, std::runtime_error,
            "Error! Unsupported degree (" << degree << ") for FE space type '" + fe_space_type_str + "' on cell '" + cell_base_name + "'.\n");
        m_fe_basis = (degree==1)
                   ? Teuchos::rcp<BasisType>(new Basis_HGRAD_WEDGE_C1_FEM<Device>())
                   : Teuchos::rcp<BasisType>(new Basis_HGRAD_WEDGE_C2_FEM<Device>());
        break;
      case FUNCTION_SPACE_HCURL:
        break;
      case FUNCTION_SPACE_HDIV:
        break;
      case FUNCTION_SPACE_HVOL:
        break;
      default:
        TEUCHOS_TEST_FOR_EXCEPTION (
          true, std::logic_error,
          "Error! Unsupported FE space type '" + fe_space_type_str + "' on cell '" + cell_base_name + "'.\n");
    }
  } else {
    TEUCHOS_TEST_FOR_EXCEPTION(
        true, std::runtime_error,
        "Error! Unsupported cell '" + cell_base_name + "'.\n");
  }

  TEUCHOS_TEST_FOR_EXCEPTION (
      m_fe_basis.is_null(), std::runtime_error,
      "Error! Something went wrong while creating the FE basis.\n" <<
      "       Here are the details:\n" <<
      "         - FE space type: " << fe_space_type_str << "\n" <<
      "         - degree: " << degree << "\n" << 
      "         - points type: " << EPointTypeToString(pointType) << "\n" <<
      "         - cell base topology: " << cell_base_name << "\n");

}

void FESpace::initialize () {
  // We cannot simply extract all entities with dimension idim (idim=0,...,spaceDim)
  // from the mesh and create dofs based on their gids, since each subcell may
  // have different number of dofs. Probably the only case where this happens
  // is with Wedges: if you have quadratic polynomials, you will have triangles
  // with 0 dofs and quadrilaterals with 1 dof.
  // The only solution is to loop on cells, and ask for the subcells of dimension
  // idim on each cell. Then query the fe basis to see how many dofs we have there.

  std::set<GO> gids_set;
  const auto& topo = *m_connectivity->getCellTopology();
  const int spaceDim = topo.getDimension();
  const int numCells = m_connectivity->getCellCount();

  // First, compute some counters. We place node dofs first, then
  // edge dofs, then face dofs, and finally cell dofs.
  // Since we loop on cells, we need to count how many dofs we have
  // at most on each subcell dimension, as well as store the max entity GID
  // of each dimension across the whole mesh. We will use these info
  // to properly offset GIDs. Notice that our offset is only an
  // upper bound in case the dof count is not constant across subcells
  // of the same dimension (like for p2 WEDGE), otherwise is an
  // exact offset (meaning no GIDs go unused, assuming the unterlying
  // mesh has used all the GIDs).

  Teuchos::Array<GO> subCellMaxGID(spaceDim+1);
  Teuchos::Array<int> numSubCells(spaceDim+1);
  Teuchos::Array<int> subCellMaxDofCount(spaceDim+1,0);
  Teuchos::Array<GO> subCellDofOffset(spaceDim+2,0);
  for (int idim=0; idim<=spaceDim; ++idim) {
    subCellMaxGID[idim] = m_connectivity->getGlobalSubcellMaxGID(idim);
    numSubCells[idim] = topo.getSubcellCount(idim);
    for (int iSubCell=0; iSubCell<numSubCells[idim]; ++iSubCell) {
      subCellMaxDofCount[idim] = std::max(subCellMaxDofCount[idim],m_fe_basis->getDofCount(idim,iSubCell));
    }
    subCellDofOffset[idim+1] = subCellDofOffset[idim] + subCellMaxGID[idim]*subCellMaxDofCount[idim];
  }

  for (int icell=0; icell<numCells; ++icell) {
    for (int idim=0; idim<spaceDim; ++idim) {
      for (int iSubCell=0; iSubCell<numSubCells[idim]; ++iSubCell) {
        // Get how many dofs are on this subcell
        const auto numLocalDofs = m_fe_basis->getDofCount(idim,iSubCell);
        if (numLocalDofs==0) {
          continue;
        }

        const GO subCellGID = m_connectivity->getSubcellGID(icell,idim,iSubCell);
        for (int idof=0; idof<numLocalDofs; ++idof) {
          const GO dof = subCellDofOffset[idim] + subCellGID*subCellMaxDofCount[idim] + idof;
          gids_set.insert(dof);
        }
      }
    }
  }

  // Now that we collected all the dofs, we need to create the vector spaces
  Teuchos::Array<GO> gids;
  gids.reserve(gids_set.size());
  for (auto gid : gids_set) {
    gids.push_back(gid);
  }

  m_overlapped_vs = createVectorSpace(m_connectivity->getComm(),gids());
  m_owned_vs      = createOneToOneVectorSpace(m_overlapped_vs);
}

} // namespace Albany
