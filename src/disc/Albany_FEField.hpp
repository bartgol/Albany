#ifndef ALBANY_FE_FIELD_HPP
#define ALBANY_FE_FIELD_HPP

#include "Albany_FESpace.hpp"

namespace Albany
{

class FEField {
public:
  FEField (const Teuchos::RCP<const FESpace> fe_space);

protected:

  Teuchos::RCP<const FESpace> m_fe_space;
};

} // namespace Albany

#endif // ALBANY_FE_FIELD_HPP
