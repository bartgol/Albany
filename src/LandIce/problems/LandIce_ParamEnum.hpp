#ifndef LANDICE_PARAM_ENUM_HPP
#define LANDICE_PARAM_ENUM_HPP

#include <string>

namespace LandIce
{

enum class ParamEnum
{
  Homotopy   ,
  Lambda     ,
  MuCoulomb  ,
  MuPowerLaw ,
  Power      ,
  Alpha      ,
  Kappa      ,
  GeoThFlux  ,
  Creep      
};

namespace ParamEnumName
{
  // Homotopy
  static const std::string HomotopyParam = "Homotopy Parameter";
  // Basal friction
  static const std::string Lambda        = "Bed Roughness";
  static const std::string MuCoulomb     = "Coulomb Friction Coefficient";
  static const std::string MuPowerLaw    = "Power Law Coefficient";
  static const std::string Power         = "Power Exponent";
  // Hydrology
  static const std::string Alpha         = "Hydraulic-Over-Hydrostatic Potential Ratio";
  static const std::string Kappa         = "Transmissivity";
  static const std::string GeoThFlux     = "Geothermal Flux";
  static const std::string Creep         = "Creep Closure Coefficient";
} // ParamEnum

} // Namespace LandIce

#endif // LANDICE_PARAM_ENUM_HPP
