#ifndef ALBANY_MESH_CONNECTIVITY_HPP
#define ALBANY_MESH_CONNECTIVITY_HPP

#include "Albany_CommTypes.hpp"
#include "Albany_ScalarOrdinalTypes.hpp"

#include <Shards_CellTopology.hpp>
#include <Teuchos_RCP.hpp>

namespace Albany {

class MeshConnectivity {
public:
  virtual ~MeshConnectivity () = default;

  Teuchos::RCP<const shards::CellTopology> getCellTopology () const;
  Teuchos::RCP<const Teuchos_Comm> getComm () const;

  int getCellCount () const;
  GO getSubcellGID (const int icell, const int subCellDim, const int iSubCell) const;
  GO getGlobalSubcellMaxGID (const int subCellDim) const;
protected:
};

} // namespace Albany

#endif // ALBANY_MESH_CONNECTIVITY_HPP
