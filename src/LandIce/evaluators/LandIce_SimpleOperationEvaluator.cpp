//*****************************************************************//
//    Albany 3.0:  Copyright 2016 Sandia Corporation               //
//    This Software is released under the BSD license detailed     //
//    in the file "license.txt" in the top-level Albany directory  //
//*****************************************************************//

#include "PHAL_AlbanyTraits.hpp"

#include "LandIce_SimpleOperationEvaluator.hpp"
#include "LandIce_SimpleOperationEvaluator_Def.hpp"

PHAL_INSTANTIATE_TEMPLATE_CLASS(LandIce::UnaryScaleOp)
PHAL_INSTANTIATE_TEMPLATE_CLASS(LandIce::UnaryLogOp)
PHAL_INSTANTIATE_TEMPLATE_CLASS(LandIce::UnaryExpOp)
PHAL_INSTANTIATE_TEMPLATE_CLASS(LandIce::UnaryLowPassOp)
PHAL_INSTANTIATE_TEMPLATE_CLASS(LandIce::UnaryHighPassOp)
PHAL_INSTANTIATE_TEMPLATE_CLASS(LandIce::UnaryBandPassOp)

PHAL_INSTANTIATE_TEMPLATE_CLASS(LandIce::BinaryScaleOp)
PHAL_INSTANTIATE_TEMPLATE_CLASS(LandIce::BinarySumOp)
PHAL_INSTANTIATE_TEMPLATE_CLASS(LandIce::BinaryLogOp)
PHAL_INSTANTIATE_TEMPLATE_CLASS(LandIce::BinaryExpOp)
PHAL_INSTANTIATE_TEMPLATE_CLASS(LandIce::BinaryLowPassOp)
PHAL_INSTANTIATE_TEMPLATE_CLASS(LandIce::BinaryHighPassOp)
PHAL_INSTANTIATE_TEMPLATE_CLASS(LandIce::BinaryBandPassFixedLowerOp)
PHAL_INSTANTIATE_TEMPLATE_CLASS(LandIce::BinaryBandPassFixedUpperOp)

PHAL_INSTANTIATE_TEMPLATE_CLASS(LandIce::TernaryBandPassOp)
