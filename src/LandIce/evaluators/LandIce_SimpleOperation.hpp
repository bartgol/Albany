//*****************************************************************//
//    Albany 3.0:  Copyright 2016 Sandia Corporation               //
//    This Software is released under the BSD license detailed     //
//    in the file "license.txt" in the top-level Albany directory  //
//*****************************************************************//

#ifndef LANDICE_SIMPLE_OPERATION_HPP
#define LANDICE_SIMPLE_OPERATION_HPP 1

#include <Teuchos_ParameterList.hpp>

#include <cmath>
#include <algorithm>

namespace LandIce
{

namespace UnaryOps
{

struct Scale
{
  void setup (const Teuchos::ParameterList& p) {
    factor = p.get<double>("Scaling Factor");
  }

  template<typename ScalarT>
  ScalarT operator() (const ScalarT& x) const {
    return factor*x;
  }

private:
  double factor;
};

struct Log
{
  void setup (const Teuchos::ParameterList& p) {
    a = p.isParameter("Factor") ? p.get<double>("Factor") : 0.0;
  }

  template<typename ScalarT>
  ScalarT operator() (const ScalarT& x) const {
    return std::log(a*x);
  }
private:
  double a;
};

struct Exp
{
  void setup (const Teuchos::ParameterList& p) {
    tau = p.isParameter("Tau") ? p.get<double>("Tau") : 0.0;
  }

  template<typename ScalarT>
  ScalarT operator() (const ScalarT& x) const {
    return std::exp(tau*x);
  }

private:
  double tau;
};

struct LowPass
{
  void setup (const Teuchos::ParameterList& p) {
    threshold_up = p.get<double>("Upper Threshold");
  }

  template<typename ScalarT>
  ScalarT operator() (const ScalarT& x) const {
    return std::min(x,threshold_up);
  }

private:
  double threshold_up;
};

struct HighPass
{
  void setup (const Teuchos::ParameterList& p) {
    threshold_lo = p.get<double>("Lower Threshold");
  }

  template<typename ScalarT>
  ScalarT operator() (const ScalarT& x) const {
    return std::max(x,threshold_lo);
  }

private:
  double threshold_lo;
};

struct BandPass
{
  void setup (const Teuchos::ParameterList& p) {
    threshold_lo = p.get<double>("Lower Threshold");
    threshold_up = p.get<double>("Upper Threshold");
  }

  template<typename ScalarT>
  ScalarT operator() (const ScalarT& x) const {
    return std::max(std::min(x,threshold_up),threshold_lo);
  }

private:
  double threshold_lo;
  double threshold_up;
};

} // namespace UnaryOps

namespace BinaryOps
{

struct Scale
{
  void setup (const Teuchos::ParameterList& /*p*/) {}
  
  template<typename Scalar1, typename Scalar2>
  auto operator() (const Scalar1& x, const Scalar2& factor) const 
  -> decltype(factor*x) {
    return factor*x;
  }
};

struct Sum
{
  void setup (const Teuchos::ParameterList& /*p*/) {}
  template<typename Scalar1, typename Scalar2>
  auto operator() (const Scalar1& x, const Scalar2& y) const
  -> decltype(x+y) {
    return x+y;
  }
};

struct Log
{
  void setup (const Teuchos::ParameterList& /*p*/) {}
  template<typename ScalarT>
  auto operator() (const ScalarT& x, const ScalarT& a) const
  -> decltype(std::log(a*x)) {
    return std::log(a*x);
  }
};

struct Exp
{
  void setup (const Teuchos::ParameterList& /*p*/) {}
  template<typename Scalar1, typename Scalar2>
  auto operator() (const Scalar1& x, const Scalar2& tau) const
  -> decltype(std::exp(tau*x)) {
    return std::exp(tau*x);
  }
};

struct LowPass
{
  void setup (const Teuchos::ParameterList& /*p*/) {}
  template<typename Scalar1, typename Scalar2>
  auto operator() (const Scalar1& x, const Scalar1& threshold_up) const
  -> decltype(std::min(threshold_up,x)) {
    return std::min(x,threshold_up);
  }
};

struct HighPass
{
  void setup (const Teuchos::ParameterList& /*p*/) {}
  template<typename Scalar1, typename Scalar2>
  auto operator() (const Scalar1& x, const Scalar2& threshold_lo) const
  -> decltype(std::max(threshold_lo,x)) {
    return std::max(x,threshold_lo);
  }
};

struct BandPassFixedUpper
{
  void setup (const Teuchos::ParameterList& p) {
    threshold_up = p.get<double>("Upper Threshold");
  }
  template<typename Scalar1, typename Scalar2>
  auto operator() (const Scalar1& x, const Scalar2& threshold_lo) const
  -> decltype(std::max(std::min(x,1.0),threshold_lo)) {
    return std::max(std::min(x,threshold_up),threshold_lo);
  }
private:
  double threshold_up;
};

struct BandPassFixedLower
{
  void setup (const Teuchos::ParameterList& p) {
    threshold_lo = p.get<double>("Lower Threshold"); 
  }
  template<typename Scalar1, typename Scalar2>
  auto operator() (const Scalar1& x, const Scalar2& threshold_up) const
  -> decltype(std::max(std::min(x,threshold_up),1.0)) {
    return std::max(std::min(x,threshold_up),threshold_lo);
  }
private:
  double threshold_lo;
};

} // namespace BinaryOps

namespace TernaryOps
{

struct BandPass
{
  template<typename Scalar1, typename Scalar2, typename Scalar3>
  auto operator() (const Scalar1& x, const Scalar2& threshold_lo, const Scalar3& threshold_up) const
  -> decltype(std::max(std::min(x,threshold_up),threshold_lo)) {
    return std::max(std::min(x,threshold_up),threshold_lo);
  }
};

} // namespace TernaryOps

} // namespace LandIce

#endif // LANDICE_SIMPLE_OPERATION_HPP
