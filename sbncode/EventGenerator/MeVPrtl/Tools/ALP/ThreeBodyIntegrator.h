/*

ThreeBodyIntegator.h

Class for conducting three-body integration

Joshua Berger 01/22/2023

*/

#include <array>
#include <functional>

class ThreeBodyIntegrator {
public:
  ThreeBodyIntegrator(std::array<double, 4> masses, std::function<double (double, double)> amp);
  ~ThreeBodyIntegrator();

  double Integrate();
  double GetSquaredAmp(double m122, double m232);
  double GetFirstIntegral(double m122);

private:
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
