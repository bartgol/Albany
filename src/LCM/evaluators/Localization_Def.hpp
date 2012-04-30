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


#include "Teuchos_TestForException.hpp"
#include "Phalanx_DataLayout.hpp"

#include "Intrepid_FunctionSpaceTools.hpp"
#include "Tensor.h"

namespace LCM {

//----------------------------------------------------------------------
template<typename EvalT, typename Traits>
Localization<EvalT, Traits>::
Localization(const Teuchos::ParameterList& p) :
  referenceCoords (p.get<std::string>                   ("Reference Coordinates Name"),
                   p.get<Teuchos::RCP<PHX::DataLayout> >("Coordinate Data Layout") ),
  currentCoords   (p.get<std::string>                   ("Current Coordinates Name"),
                   p.get<Teuchos::RCP<PHX::DataLayout> >("Coordinate Data Layout") ),
  cubature        (p.get<Teuchos::RCP<Intrepid::Cubature<RealType> > >("Cubature")),
  intrepidBasis   (p.get<Teuchos::RCP<Intrepid::Basis<RealType, Intrepid::FieldContainer<RealType> > > > ("Intrepid Basis") ),
  cellType        (p.get<Teuchos::RCP<shards::CellTopology> > ("Cell Type")),
  defGrad         (p.get<std::string>                   ("Deformation Gradient Name"),
                   p.get<Teuchos::RCP<PHX::DataLayout> >("QP Tensor Data Layout"))
{
  this->addDependentField(referenceCoords);
  this->addDependentField(currentCoords);
  this->addEvaluatedField(defGrad);

  // Get Dimensions
  Teuchos::RCP<PHX::DataLayout> vert_dl = p.get< Teuchos::RCP<PHX::DataLayout> >("Coordinate Data Layout");
  std::vector<PHX::DataLayout::size_type> dims;
  vert_dl->dimensions(dims);

  int containerSize = dims[0];
  numNodes = dims[1];
  numPlaneNodes = numNodes/2;

  Teuchos::RCP<PHX::DataLayout> defGrad_dl = p.get< Teuchos::RCP<PHX::DataLayout> >("QP Tensor Data Layout");
  defGrad_dl->dimensions(dims);
  numQPs = dims[1];
  numDims = dims[2];

  // Allocate Temporary FieldContainers
  refValues.resize(numNodes, numQPs);
  refGrads.resize(numNodes, numQPs, numDims);
  refPoints.resize(numQPs, numDims);
  refWeights.resize(numQPs);

  // new stuff
  midplaneCoords.resize(containerSize, numPlaneNodes, numDims);
  bases.resize(containerSize, numQPs, numDims, numDims);
  dualBases.resize(containerSize, numQPs, numDims, numDims);
  jacobian.resize(containerSize, numQPs);
  normals.resize(containerSize, numQPs, numDims);
  gap.resize(containerSize, numQPs, numDims);

  // Pre-Calculate reference element quantitites
  std::cout << "Calling Intrepid to get reference quantities" << std::endl;
  cubature->getCubature(refPoints, refWeights);
  intrepidBasis->getValues(refValues, refPoints, Intrepid::OPERATOR_VALUE);
  intrepidBasis->getValues(refGrads, refPoints, Intrepid::OPERATOR_GRAD);
  
  this->setName("Localization"+PHX::TypeString<EvalT>::value);
}

//----------------------------------------------------------------------
template<typename EvalT, typename Traits>
void Localization<EvalT, Traits>::
postRegistrationSetup(typename Traits::SetupData d,
                      PHX::FieldManager<Traits>& fm)
{
  this->utils.setFieldData(referenceCoords,fm);
  this->utils.setFieldData(currentCoords,fm);
  this->utils.setFieldData(defGrad,fm);
}

//----------------------------------------------------------------------
template<typename EvalT, typename Traits>
void Localization<EvalT, Traits>::
evaluateFields(typename Traits::EvalData workset)
{

  for (std::size_t cell(0); cell < workset.numCells; ++cell) 
  {
    // for the reference geometry
    // compute the mid-plane coordinates
    computeMidplaneCoords(referenceCoords, midplaneCoords);

    // compute base vectors
    computeBaseVectors(midplaneCoords, bases);
    
    // compute the dual
    computeDualBaseVectors(midplaneCoords, bases, normals, dualBases);

    // compute the Jacobian
    computeJacobian(bases, dualBases, area, jacobian);
    
    // for the current configuration

    // compute the mid-plane coordinates
    computeMidplaneCoords(currentCoords, midplaneCoords);

    // compute base vectors
    computeBaseVectors(midplaneCoords, bases);
    
    // compute gap
    computeGap(currentCoords, gap);

    // compute deformation gradient
    computeDeformationGradient(bases);

    // call constitutive response

    // compute force
  }
}
//----------------------------------------------------------------------
template<typename EvalT, typename Traits>
void Localization<EvalT, Traits>::
computeMidplaneCoords(PHX::MDField<ScalarT,Cell,Vertex,Dim> coords,
                      FC & midplaneCoords)
{
  std::cout << "In computeMidplaneCoords" << std::endl;

  std::cout << "Dimensions of coords : " 
            << coords.dimension(0) << " " 
            << coords.dimension(1) << " " 
            << coords.dimension(2) << " " 
            << std::endl;
  
  for (int cell(0); cell < midplaneCoords.dimension(0); ++cell) 
  {
    std::cout << "Cell    : " << cell << std::endl;
    // compute the mid-plane coordinates
    for (int node(0); node < numPlaneNodes; ++node)
    {
      int topNode = node + numPlaneNodes;
      std::cout << "node   : " << node << std::endl;
      std::cout << "topNode: " << topNode << std::endl;
      for (int dim(0); dim < numDims; ++dim)
      {
        midplaneCoords(cell, node, dim) = 0.5 * ( coords(cell, node, dim) + coords(cell, topNode, dim) );
      }
    }
  }  
}
//----------------------------------------------------------------------
template<typename EvalT, typename Traits>
void Localization<EvalT, Traits>::
computeBaseVectors(const FC & midplaneCoords, FC & bases)
{
  std::cout << "In computeBaseVectors" << std::endl;
  typedef LCM::Vector<ScalarT> V;

  for (int cell(0); cell < midplaneCoords.dimension(0); ++cell)
  {
    // get the midplane coordinates
    std::vector<LCM::Vector<ScalarT> > midplaneNodes(numPlaneNodes);
    for (std::size_t node(0); node < numPlaneNodes; ++node)
      midplaneNodes.push_back(V(midplaneCoords(cell,node,0),
                                midplaneCoords(cell,node,1),
                                midplaneCoords(cell,node,2)));

    V g0(0.0), g1(0.0), g2(0.0);
    //compute the base vectors
    for (std::size_t pt(0); pt < numQPs; ++pt)
    {
      for (std::size_t node(0); node < numPlaneNodes; ++ node)
      {
        g0 += ScalarT(refGrads(node, pt, 0)) * midplaneNodes[node];
        g1 += ScalarT(refGrads(node, pt, 1)) * midplaneNodes[node];
      }
      g2 = cross(g1,g2)/norm(cross(g1,g2));

      bases(cell,pt,0,0) = g0(0); bases(cell,pt,0,1) = g0(1); bases(cell,pt,0,2) = g0(2);
      bases(cell,pt,1,0) = g1(0); bases(cell,pt,1,1) = g1(1); bases(cell,pt,1,2) = g1(2);
      bases(cell,pt,2,0) = g2(0); bases(cell,pt,2,1) = g2(1); bases(cell,pt,2,2) = g2(2);
    }
  }
}
//----------------------------------------------------------------------
template<typename EvalT, typename Traits>
void Localization<EvalT, Traits>::
computeDualBaseVectors(const FC & midplaneCoords, const FC & bases, FC & normals, FC & dualBases)
{
  std::cout << "In computeDualBaseVectors" << std::endl;
  typedef LCM::Vector<ScalarT> V;
  std::size_t worksetSize = midplaneCoords.dimension(0);

  V g0(0.0), g1(0.0), g2(0.0), G0(0.0), G1(0.0), G2(0.0);
  
  const ScalarT * dataPtr;

  for (std::size_t cell(0); cell < worksetSize; ++cell)
  {
    for (std::size_t pt(0); pt < numQPs; ++pt)
    {
      dataPtr = &bases(cell,pt,0,0);
      g0 = V( dataPtr );
      dataPtr = &bases(cell,pt,1,0);
      g1 = V( dataPtr );
      dataPtr = &bases(cell,pt,1,0);
      g2 = &bases(cell,pt,2,0);
      
      G0 = cross( g1,g2 ) / dot( g0, cross( g1,g2 ) );
      G1 = cross( g0,g2 ) / dot( g1, cross( g0,g2 ) );
      G2 = cross( g0,g1 ) / dot( g2, cross( g0,g1 ) );

      dualBases(cell,pt,0,0) = G0(0); dualBases(cell,pt,0,1) = G0(1); dualBases(cell,pt,0,2) = G0(2);
      dualBases(cell,pt,1,0) = G1(0); dualBases(cell,pt,1,1) = G1(1); dualBases(cell,pt,1,2) = G1(2);
      dualBases(cell,pt,2,0) = G2(0); dualBases(cell,pt,2,1) = G2(1); dualBases(cell,pt,2,2) = G2(2);
    }
  }
}
//----------------------------------------------------------------------
template<typename EvalT, typename Traits>
void Localization<EvalT, Traits>::
computeJacobian(const FC & bases, const FC & dualbases, FC & area, FC & jacobian)
{
  std::cout << "In computeJacobian" << std::endl;
}

//----------------------------------------------------------------------
template<typename EvalT, typename Traits>
void Localization<EvalT, Traits>::
computeGap(const PHX::MDField<ScalarT,Cell,Vertex,Dim> coords, FC & gap)
{
  std::cout << "In computeGap" << std::endl;
}
//----------------------------------------------------------------------
template<typename EvalT, typename Traits>
void Localization<EvalT, Traits>::
computeDeformationGradient(const FC & bases)
{
  std::cout << "In computeDeformationGradient" << std::endl;
  std::size_t worksetSize = bases.dimension(0);
  for (std::size_t cell(0); cell < worksetSize; ++cell)
  {
    for (std::size_t pt(0); pt < numQPs; ++pt)
    {
      for (std::size_t dim1(0); dim1 < numDims; ++dim1)
      {
        for (std::size_t dim2(0); dim2 < numDims; ++dim2)
        {
          defGrad(cell,pt,dim1,dim2) = 0.0;
          if (dim1 == dim2)
            defGrad(cell,pt,dim1,dim2) = 1.0;;
        }
      }
    }
  }
}
//----------------------------------------------------------------------
} //namespace LCM
