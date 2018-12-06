//*****************************************************************//
//    Albany 3.0:  Copyright 2016 Sandia Corporation               //
//    This Software is released under the BSD license detailed     //
//    in the file "license.txt" in the top-level Albany directory  //
//*****************************************************************//

#include "Albany_TikhonovRegularizationResponseFunction.hpp"

#include "PHAL_AlbanyTraits.hpp"

namespace Albany {

TikhonovRegularizationResponseFunction::
TikhonovRegularizationResponseFunction (const Teuchos::ParameterList& params,
                                        const Teuchos::RCP<const Teuchos_Comm>& comm)
 : ScalarResponseFunction(comm)
{
  m_regCoeff   = params.get<double>("Regularization Coefficient");
  m_paramIndex = params.get<int>("Parameter Index");
  m_refValues  = params.get<Teuchos::Array<ST>>("Reference Values");
}

void TikhonovRegularizationResponseFunction::evaluateResponse(
    const double /*current_time*/,
    const Teuchos::RCP<const Thyra_Vector>& /* x */,
    const Teuchos::RCP<const Thyra_Vector>& /*xdot*/,
    const Teuchos::RCP<const Thyra_Vector>& /*xdotdot*/,
		const Teuchos::Array<ParamVec>& params,
		const Teuchos::RCP<Thyra_Vector>& g)
{
  const ParamVec& p = params[m_paramIndex];
  ST reg = 0.0;
  for (unsigned int i=0; i<p.size(); ++i) {
    reg += std::pow((p[i].family->getValue<PHAL::AlbanyTraits::Residual>() - m_refValues[i]),2);
  }
  g->assign(reg*m_regCoeff/2.0);
}

void TikhonovRegularizationResponseFunction::evaluateTangent(
    const double /* alpha */,
		const double /*beta*/,
		const double /*omega*/,
		const double current_time,
		bool /*sum_derivs*/,
    const Teuchos::RCP<const Thyra_Vector>& x,
    const Teuchos::RCP<const Thyra_Vector>& xdot,
    const Teuchos::RCP<const Thyra_Vector>& xdotdot,
    const Teuchos::Array<ParamVec>& params,
		ParamVec* /*deriv_p*/,
    const Teuchos::RCP<const Thyra_MultiVector>& /* Vx */,
    const Teuchos::RCP<const Thyra_MultiVector>& /*Vxdot*/,
    const Teuchos::RCP<const Thyra_MultiVector>& /*Vxdotdot*/,
    const Teuchos::RCP<const Thyra_MultiVector>& /*Vp*/,
    const Teuchos::RCP<Thyra_Vector>& g,
    const Teuchos::RCP<Thyra_MultiVector>& gx,
    const Teuchos::RCP<Thyra_MultiVector>& gp)
{
  // Evaluate response g
  if (!g.is_null()) {
    evaluateResponse(current_time,x,xdot,xdotdot,params,g);
  }

  // Evaluate tangent of g = dg/dx*dx/dp + dg/dxdot*dxdot/dp + dg/dp

  // dg/dx = 0
  if (!gx.is_null()) {
    gx->assign(0.0);
  }

  // dg/dp = m_regCoeff*(params[i]-m_refValues)
  if (!gp.is_null()) {
    const ParamVec& p = params[m_paramIndex];
    for (unsigned int i=0; i<p.size(); ++i) {
      gp->col(i)->assign((p[i].family->getValue<PHAL::AlbanyTraits::Residual>() - m_refValues[i])*m_regCoeff);
    }
  }
}

void TikhonovRegularizationResponseFunction::evaluateGradient(
    const double current_time,
    const Teuchos::RCP<const Thyra_Vector>& x,
    const Teuchos::RCP<const Thyra_Vector>& xdot,
    const Teuchos::RCP<const Thyra_Vector>& xdotdot,
		const Teuchos::Array<ParamVec>& params,
		ParamVec* /*deriv_p*/,
		const Teuchos::RCP<Thyra_Vector>& g,
		const Teuchos::RCP<Thyra_MultiVector>& dg_dx,
		const Teuchos::RCP<Thyra_MultiVector>& dg_dxdot,
		const Teuchos::RCP<Thyra_MultiVector>& dg_dxdotdot,
		const Teuchos::RCP<Thyra_MultiVector>& dg_dp)
{
  // Evaluate response g
  if (!g.is_null()) {
    evaluateResponse(current_time,x,xdot,xdotdot,params,g);
  }
  
  // dg/dx = 0
  if (!dg_dx.is_null()) {
    dg_dx->assign(0.0);
  }

  // dg/dxdot = 0
  if (!dg_dxdot.is_null()) {
    dg_dxdot->assign(0.0);
  }

  // dg/dxdotdot = 0
  if (!dg_dxdotdot.is_null()) {
    dg_dxdotdot->assign(0.0);
  }

  // dg/dp = m_regCoeff*(params[i]-m_refValues)
  if (!dg_dp.is_null()) {
    const ParamVec& p = params[m_paramIndex];
    for (unsigned int i=0; i<p.size(); ++i) {
      dg_dp->col(i)->assign((p[i].family->getValue<PHAL::AlbanyTraits::Residual>() - m_refValues[i])*m_regCoeff);
    }
  }
}

//! Evaluate distributed parameter derivative dg/dp = 0
void TikhonovRegularizationResponseFunction::evaluateDistParamDeriv(
    const double /*current_time*/,
    const Teuchos::RCP<const Thyra_Vector>& /*x*/,
    const Teuchos::RCP<const Thyra_Vector>& /*xdot*/,
    const Teuchos::RCP<const Thyra_Vector>& /*xdotdot*/,
    const Teuchos::Array<ParamVec>& /*param_array*/,
    const std::string& /*dist_param_name*/,
    const Teuchos::RCP<Thyra_MultiVector>& dg_dp)
{
  if (!dg_dp.is_null()) {
    dg_dp->assign(0.0);
  }
}

} // namespace Albany
