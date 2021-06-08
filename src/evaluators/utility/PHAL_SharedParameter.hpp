#ifndef PHAL_SHARED_PARAMETER_HPP
#define PHAL_SHARED_PARAMETER_HPP 1

#include "PHAL_Dimension.hpp"
#include "Albany_SacadoTypes.hpp"
#include "Albany_Utils.hpp"

#include "Phalanx_Evaluator_WithBaseImpl.hpp"
#include "Phalanx_Evaluator_Derived.hpp"
#include "Sacado_ParameterAccessor.hpp"

//IKT 6/3/2020 TODO: implement support for vector parameters, which is not available currently.

namespace PHAL
{

template<typename EvalT, typename Traits, typename ParamNameEnum, ParamNameEnum ParamName>
class SharedParameter : public PHX::EvaluatorWithBaseImpl<Traits>,
                        public PHX::EvaluatorDerived<EvalT, Traits>,
                        public Sacado::ParameterAccessor<EvalT, SPL_Traits>
{
public:

  typedef typename EvalT::ScalarT   ScalarT;
  typedef ParamNameEnum             EnumType;

  SharedParameter (const Teuchos::ParameterList& p, const Teuchos::RCP<Albany::Layouts>& dl)
  {
    param_name   = p.get<std::string>("Parameter Name");
    param_as_field = PHX::MDField<ScalarT,Dim>(param_name,dl->shared_param);

    // Never actually evaluated, but creates the evaluation tag
    this->addEvaluatedField(param_as_field);

    // Sacado-ized parameter
    Teuchos::RCP<ParamLib> paramLib = p.get<Teuchos::RCP<ParamLib>>("Parameter Library");
    this->registerSacadoParameter(param_name, paramLib);
    this->setName("Shared Parameter " + param_name + PHX::print<EvalT>());
  }

  void postRegistrationSetup(typename Traits::SetupData d, PHX::FieldManager<Traits>& fm)
  {
    this->utils.setFieldData(param_as_field,fm);
    d.fill_field_dependencies(this->dependentFields(),this->evaluatedFields(),false);
  }

  static void setNominalValue (const Teuchos::ParameterList& p, double default_value);

  static ScalarT getValue ()
  {
    return value;
  }

  ScalarT& getValue(const std::string &n)
  {
    if (n==param_name)
      return value;

    return dummy;
  }

  void evaluateFields(typename Traits::EvalData /*d*/)
  {
    param_as_field(0) = value;
  }

protected:

  static ScalarT              value;
  static ScalarT              dummy;
  static std::string          param_name;

  PHX::MDField<ScalarT,Dim>   param_as_field;
};

template<typename EvalT, typename Traits, typename ParamNameEnum, ParamNameEnum ParamName>
void SharedParameter<EvalT,Traits,ParamNameEnum,ParamName>::
setNominalValue (const Teuchos::ParameterList& p, double default_value)
{
  // First we scan the Parameter list to see if this parameter is listed in it,
  // in which case we use the nominal value.
  bool found = false;
  if (p.isParameter("Number Of Parameters"))
  {
    int n = p.get<int>("Number Of Parameters");
    for (int i=0; (found==false) && i<n; ++i)
    {
      const Teuchos::ParameterList& pvi = p.sublist(Albany::strint("Parameter",i));
      std::string parameterType = "Scalar";
      if(pvi.isParameter("Type"))
        parameterType = pvi.get<std::string>("Type");
      if (parameterType == "Distributed")
        break; // Pointless to check the remaining parameters as they are all distributed

      if (parameterType == "Scalar") {
        if (!pvi.isParameter("Nominal Value"))
          continue; // Pointless to check the parameter names, since we don't have nominal values
        if (pvi.get<std::string>("Name")==param_name)
        {
          double nom_val = pvi.get<double>("Nominal Value");
          value = nom_val;
          found = true;
          break;
        }
      }
      else {
        int m = pvi.get<int>("Dimension");
        for (int j=0; j<m; ++j)
        {
          const Teuchos::ParameterList& pj = pvi.sublist(Albany::strint("Scalar",j));
          if (!pj.isParameter("Nominal Value"))
            continue; // Pointless to check the parameter names, since we don't have nominal values
          if (pj.get<std::string>("Name")==param_name)
          {
            double nom_val = pj.get<double>("Nominal Value");
            value = nom_val;
            found = true;
            break;
          }
        }
      }
    }
  }
  else if (p.isParameter("Number") && !p.isParameter("Nominal Values"))
  {
    int m = p.get<int>("Number");
    for (int j=0; j<m; ++j)
    {
      if (p.get<std::string>(Albany::strint("Parameter",j))==param_name)
      {
        Teuchos::Array<double> nom_vals = p.get<Teuchos::Array<double>>("Nominal Values");
        value = nom_vals[j];
        found = true;
        break;
      }
    }
  }

  if (!found)
  {
    value = default_value;
  }

  dummy = 0;
}

template<typename EvalT, typename Traits, typename ParamNameEnum, ParamNameEnum ParamName>
typename EvalT::ScalarT SharedParameter<EvalT,Traits,ParamNameEnum,ParamName>::value;

template<typename EvalT, typename Traits, typename ParamNameEnum, ParamNameEnum ParamName>
typename EvalT::ScalarT SharedParameter<EvalT,Traits,ParamNameEnum,ParamName>::dummy;

template<typename EvalT, typename Traits, typename ParamNameEnum, ParamNameEnum ParamName>
std::string SharedParameter<EvalT,Traits,ParamNameEnum,ParamName>::param_name;

} // Namespace PHAL

#endif // PHAL_SHARED_PARAMETER_HPP
