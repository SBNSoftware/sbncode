/*

ThreeBodyIntegator.h

Class for conducting three-body integration

Joshua Berger 01/22/2023

*/

#include <array>
#include <functional>

// Integrate a 1->3 decay width over the relevant phase space,
// defined in terms of the lorentz invariants m_{12}^2, m_{23}^2
// where "0" is the parent particle and {1,2,3} are the daughters
class ThreeBodyIntegrator {
public:
  // Inputs: 
  //   array of 4 masses m0, m1, m2, m3
  //   function taking m_{12}^2 and m_{23}^2 are returning the amplitude
  ThreeBodyIntegrator(std::array<double, 4> masses, std::function<double (double, double)> amp);
  ~ThreeBodyIntegrator();

  // Do the integral
  double Integrate();
  // Integrand
  double GetSquaredAmp(double m122, double m232);
  double GetFirstIntegral(double m122);

private:
  // Internal helper methods
  double GetE2st(double m122);
  double GetE3st(double m122);
  double M232Min(double m122);
  double M232Max(double m122);

  std::array<double, 4> fMasses;
  std::function<double (double, double)> fAmp;
  double fM122Min;
  double fM122Max;
  double fWidth;
};
