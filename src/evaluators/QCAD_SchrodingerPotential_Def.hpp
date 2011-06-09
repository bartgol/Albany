/********************************************************************\
*            Albany, Copyright (2010) Sandia Corporation             *
*                                                                    *
* Notice: This computer software was prepared by Sandia Corporation, *
* hereinafter the Contractor, under Contract DE-AC04-94AL85000 with  *
* the Department of Energy (DOE). All rights in the computer software*
* are reserved by DOE on behalf of the United States Government and  *
* the Contractor as provided in the Contract. You are authorized to  *
* use this computer software for Governmental purposes but it is not *
* to be released or distributed to the public. NEITHER THE GOVERNMENT*
* NOR THE CONTRACTOR MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR      *
* ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE. This notice    *
* including this sentence must appear on any copies of this software.*
*    Questions to Andy Salinger, agsalin@sandia.gov                  *
\********************************************************************/


#include <fstream>
#include "Teuchos_TestForException.hpp"
#include "Phalanx_DataLayout.hpp"
#include "Sacado_ParameterRegistration.hpp"

template<typename EvalT, typename Traits>
QCAD::SchrodingerPotential<EvalT, Traits>::
SchrodingerPotential(Teuchos::ParameterList& p) :
  coordVec(p.get<std::string>("QP Coordinate Vector Name"),
     p.get<Teuchos::RCP<PHX::DataLayout> >("QP Vector Data Layout")),
  psi(p.get<std::string>("QP Variable Name"),
      p.get<Teuchos::RCP<PHX::DataLayout> >("QP Scalar Data Layout")),
  V(p.get<std::string>("QP Potential Name"),
      p.get<Teuchos::RCP<PHX::DataLayout> >("QP Scalar Data Layout"))
{
  Teuchos::ParameterList* psList = p.get<Teuchos::ParameterList*>("Parameter List");

  Teuchos::RCP<const Teuchos::ParameterList> reflist = 
      this->getValidSchrodingerPotentialParameters();
  psList->validateParameters(*reflist,0);

  Teuchos::RCP<PHX::DataLayout> vector_dl =
      p.get< Teuchos::RCP<PHX::DataLayout> >("QP Vector Data Layout");
  std::vector<PHX::DataLayout::size_type> dims;
  vector_dl->dimensions(dims);
  numQPs  = dims[1];
  numDims = dims[2];

  energy_unit_in_eV = p.get<double>("Energy unit in eV");
  length_unit_in_m = p.get<double>("Length unit in m");
  potentialType = psList->get("Type", "defaultType");
  E0 = psList->get("E0", 1.0);
  scalingFactor = psList->get("Scaling Factor", 1.0);

  potentialStateName = p.get<std::string>("QP Potential Name");

  // Add E0 as a Sacado-ized parameter
  Teuchos::RCP<ParamLib> paramLib =
      p.get< Teuchos::RCP<ParamLib> >("Parameter Library", Teuchos::null);
  new Sacado::ParameterRegistration<EvalT, SPL_Traits>(
      "Schrodinger Potential E0", this, paramLib);
  new Sacado::ParameterRegistration<EvalT, SPL_Traits> (
      "Schrodinger Potential Scaling Factor", this, paramLib);

  this->addDependentField(psi);
  this->addDependentField(coordVec);

  this->addEvaluatedField(V);
  this->setName("Schrodinger Potential"+PHX::TypeString<EvalT>::value);
}

// **********************************************************************
template<typename EvalT, typename Traits>
void QCAD::SchrodingerPotential<EvalT, Traits>::
postRegistrationSetup(typename Traits::SetupData d,
                      PHX::FieldManager<Traits>& fm)
{
  this->utils.setFieldData(V,fm);
  this->utils.setFieldData(psi,fm);
  this->utils.setFieldData(coordVec,fm);
}

// **********************************************************************
template<typename EvalT, typename Traits>
void QCAD::SchrodingerPotential<EvalT, Traits>::
evaluateFields(typename Traits::EvalData workset)
{

  //Parabolic potential (test case)
  if (potentialType == "Parabolic") 
  {
    for (std::size_t cell=0; cell < workset.numCells; ++cell) {
      for (std::size_t qp=0; qp < numQPs; ++qp) {
        V(cell, qp) = parabolicPotentialValue(numDims, &coordVec(cell,qp,0));
      }
    }
  }


  //Potential taken from input state "Electric Potential" - now hardcoded, maybe allow variable later?
  else if (potentialType == "FromState") 
  {
    //ANDY - can't access oldState due to const qualifiers, so not just use newstate (since both hold re-init data)
    //const Albany::StateVariables& oldState = *workset.oldState;
    //Intrepid::FieldContainer<RealType>& potentialState = *oldState[potentialStateName];

    Albany::StateVariables& newState = *workset.newState;
    Intrepid::FieldContainer<RealType>& potentialState = *newState[potentialStateName];
    for (std::size_t cell=0; cell < workset.numCells; ++cell) {
      for (std::size_t qp=0; qp < numQPs; ++qp) {
        double d =  potentialState(cell, qp);
        V(cell, qp) = scalingFactor * d;

	//HACK to help anasazi solve
	//if(workset.EBName == "silicon" || scalingFactor < 0) {
	//  V(cell, qp) =  d;
	//}
	//else V(cell, qp) = scalingFactor;
      }
    }
  }    

  /***** otherwise error? ******/
  else 
  {
    TEST_FOR_EXCEPT(true);
  }
}

// **********************************************************************
template<typename EvalT,typename Traits>
typename QCAD::SchrodingerPotential<EvalT,Traits>::ScalarT& 
QCAD::SchrodingerPotential<EvalT,Traits>::getValue(const std::string &n)
{
  if(n == "Schrodinger Potential Scaling Factor") return scalingFactor;
  else if(n == "Schrodinger Potential E0") return E0;
  else TEST_FOR_EXCEPT(true); return E0; //dummy so all control paths return
}

// **********************************************************************
template<typename EvalT,typename Traits>
Teuchos::RCP<const Teuchos::ParameterList>
QCAD::SchrodingerPotential<EvalT,Traits>::getValidSchrodingerPotentialParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> validPL =
       rcp(new Teuchos::ParameterList("Valid Schrodinger Potential Params"));;

  validPL->set<string>("Type", "defaultType", "Switch between different potential types");
  validPL->set<double>("E0", 1.0, "Energy scale - dependent on type");
  validPL->set<double>("Scaling Factor", 1.0, "Constant scaling factor");

  return validPL;
}

// **********************************************************************

//Return potential in energy_unit_in_eV * eV units
template<typename EvalT,typename Traits>
typename QCAD::SchrodingerPotential<EvalT,Traits>::ScalarT
QCAD::SchrodingerPotential<EvalT,Traits>::parabolicPotentialValue(
    const int numDim, const MeshScalarT* coord)
{
  ScalarT val;
  MeshScalarT r2;
  int i;

  //std::cout << "x = " << coord[0] << endl; //in 1D, x-coords run from zero to 1

  /***** define universal constants as double constants *****/
  const double hbar = 1.0546e-34;  // Planck constant [J s]
  const double emass = 9.1094e-31; // Electron mass [kg]
  const double evPerJ = 6.2415e18; // eV per Joule (J/eV)

  const double parabolicFctr = 0.5 * (emass / (evPerJ * hbar*hbar) );

  // prefactor from constant, including scaling due to units
  ScalarT prefactor;  
  
  prefactor = parabolicFctr * E0*E0 * (energy_unit_in_eV * pow(length_unit_in_m,2));  
  for(i=0, r2=0.0; i<numDim; i++)
    r2 += (coord[i]-0.5)*(coord[i]-0.5);
  val = prefactor * r2;  
  
  return scalingFactor * val;
}


// *****************************************************************************
