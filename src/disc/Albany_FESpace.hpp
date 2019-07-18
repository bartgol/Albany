#ifndef ALBANY_FE_SPACE_HPP
#define ALBANY_FE_SPACE_HPP

#include "Albany_MeshConnectivity.hpp"
#include "Albany_CommTypes.hpp"
#include "Albany_ThyraTypes.hpp"
#include "Albany_KokkosTypes.hpp"

// Phalanx determines the Kokkos node we use for Tpetra types
#include <Intrepid2_Basis.hpp>
#include <Teuchos_RCP.hpp>

namespace Albany {

class FESpace {
public:
  using EFunctionSpace = Intrepid2::EFunctionSpace;
  using EPointType     = Intrepid2::EPointType;
  using BasisType = Intrepid2::Basis<PHX::Device,RealType,RealType>;

  FESpace (Teuchos::RCP<const MeshConnectivity> connectivity,
           const EFunctionSpace                 type,
           const int                            degree,
           const EPointType                     pointType = Intrepid2::POINTTYPE_EQUISPACED);

  EFunctionSpace  getType   () const { return m_type; }
  int             getdegree () const { return m_fe_basis->getDegree(); }

  Teuchos::RCP<const Thyra_VectorSpace> getVectorSpace ();
  Teuchos::RCP<const Thyra_VectorSpace> getOverlappedVectorSpace ();

protected:

  void createBasis (const int degree, const EPointType pointType);
  void initialize ();
  void addNodalFEs (Teuchos::Array<GO>& gids);

  Teuchos::RCP<const MeshConnectivity>  m_connectivity;
  const EFunctionSpace                  m_type;
  Teuchos::RCP<const Teuchos_Comm>      m_comm;

  Teuchos::RCP<const Thyra_VectorSpace>   m_owned_vs;
  Teuchos::RCP<const Thyra_VectorSpace>   m_overlapped_vs;

  Teuchos::RCP<const BasisType>           m_fe_basis;
};

} // namespace Albany

#endif // ALBANY_FE_SPACE_HPP
