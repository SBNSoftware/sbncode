/*

ThreeBodyIntegator.cc

Class for conducting three-body integration

Joshua Berger 01/22/2023

*/

#include "ThreeBodyIntegrator.h"
#include <gsl/gsl_integration.h>

namespace IntegratorWrapper {
  extern "C" {
    struct inner_params {
      double m122;
      ThreeBodyIntegrator *integrator;
    };

    double SquaredAmp_wrapper(double m232, void * params) {
      inner_params& this_params = *(inner_params *)params;
      double m122 = this_params.m122;
      ThreeBodyIntegrator *this_integrator = this_params.integrator;
      double result = this_integrator->GetSquaredAmp(m122, m232);
      return result;
    }

    double FirstIntegral_wrapper(double m122, void * params) {
      ThreeBodyIntegrator *this_integrator = (ThreeBodyIntegrator *)params;
      double result = this_integrator->GetFirstIntegral(m122);
      return result;
    }
  }
}


ThreeBodyIntegrator::ThreeBodyIntegrator(std::array<double, 4> masses, std::function<double (double, double)> amp):
  fAmp(amp)
{
  fMasses = masses;
  fM122Min = pow(masses[1] + masses[2], 2.0);
  fM122Max = pow(masses[0] - masses[3], 2.0);
  fWidth = -1.0;
}

ThreeBodyIntegrator::~ThreeBodyIntegrator() {
}

double ThreeBodyIntegrator::GetE2st(double m122) {
  double m1 = fMasses[1];
  double m2 = fMasses[2];
  return 0.5 * (m122 - pow(m1, 2.0) + pow(m2, 2.0)) / sqrt(m122);
}

double ThreeBodyIntegrator::GetE3st(double m122) {
  double M = fMasses[0];
  double m3 = fMasses[3];
  return 0.5 * (pow(M, 2.0) - m122 - pow(m3, 2.0)) / sqrt(m122);
}

double ThreeBodyIntegrator::M232Min(double m122) {
  double m2 = fMasses[2];
  double m3 = fMasses[3];
  double E2st = GetE2st(m122);
  double E3st = GetE3st(m122);
  return pow(E2st + E3st, 2.0) - pow(sqrt(E2st*E2st - m2*m2) + sqrt(E3st*E3st - m3*m3), 2.0);
}

double ThreeBodyIntegrator::M232Max(double m122) {
  double m2 = fMasses[2];
  double m3 = fMasses[3];
  double E2st = GetE2st(m122);
  double E3st = GetE3st(m122);
  return pow(E2st + E3st, 2.0) - pow(sqrt(E2st*E2st - m2*m2) - sqrt(E3st*E3st - m3*m3), 2.0);
}

double ThreeBodyIntegrator::GetSquaredAmp(double m122, double m232) {
  double amp = fAmp(m122, m232);
  return amp*amp;
}

double ThreeBodyIntegrator::GetFirstIntegral(double m122) {
  double m232min = M232Min(m122);
  double m232max = M232Max(m122);
  gsl_integration_workspace * w = gsl_integration_workspace_alloc(1000);
  gsl_function F;
  F.function = &IntegratorWrapper::SquaredAmp_wrapper;
  IntegratorWrapper::inner_params params;
  params.m122 = m122;
  params.integrator = this;
  F.params = &params;
  double result, error;
  gsl_integration_qags(&F, m232min, m232max, 0, 1.0e-7, 1000, w, &result, &error);
  gsl_integration_workspace_free(w);
  return result;
}

double ThreeBodyIntegrator::Integrate() {
  if (fWidth >= 0.0)
    return fWidth;
  else {
    double M = fMasses[0];
    gsl_integration_workspace * w = gsl_integration_workspace_alloc(1000);
    gsl_function F;
    F.function = &IntegratorWrapper::FirstIntegral_wrapper;
    F.params = (void *)this;
    double result, error;
    gsl_integration_qags(&F, fM122Min, fM122Max, 0, 1.0e-7, 1000, w, &result, &error);
    gsl_integration_workspace_free(w);
    fWidth = result / (128.0 * pow(M_PI, 3.0) * pow(M, 2.0));
    return fWidth;
  }
}
