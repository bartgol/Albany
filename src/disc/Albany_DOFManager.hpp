#ifndef ALBANY_DOF_MANAGER_HPP
#define ALBANY_DOF_MANAGER_HPP

#include "Albany_MeshConnectivity.hpp"
#include <Albany_CommTypes.hpp>

#include <Teuchos_RCP.hpp>

namespace Albany
{

enum class DOFType {
    
};

class DOFManager {
public:
  DOFManager (Teuchos::RCP<const MeshConnectivity>  connectivity,
              Teuchos::RCP<const Teuchos_Comm>      comm);

  void add_field (const std::string& name);
private:

  Teuchos::RCP<const MeshConnectivity>    m_connectivity;
  Teuchos::RCP<const Teuchos_Comm>        m_comm;
};

} // namespace Albany

#endif // ALBANY_DOF_MANAGER_HPP
