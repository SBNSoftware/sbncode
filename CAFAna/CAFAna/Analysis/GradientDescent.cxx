#include "CAFAna/Analysis/GradientDescent.h"

#include "Minuit2/FCNGradientBase.h"
#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnUserParameters.h"

#include <iostream>

// Globals set to the default settings
ROOT::Minuit2::MnUserParameterState gState;
ROOT::Minuit2::MnStrategy gStrat;

namespace ana
{
  //----------------------------------------------------------------------
  GradientDescent::
  GradientDescent(const ROOT::Minuit2::FCNGradientBase& func,
                  const ROOT::Minuit2::MnUserParameters& pars)
    : ROOT::Minuit2::MnApplication(func, gState, gStrat),
      fFunc(func), fPars(pars)
  {
  }

  //----------------------------------------------------------------------
  ROOT::Minuit2::FunctionMinimum GradientDescent::
  operator()(unsigned int /*maxfcn*/, double /*tolerance*/)
  {
    // Initialize at the input parameters
    std::vector<double> pt = fPars.Params();
    const unsigned int N = pt.size();

    // Use the user errors as a crude way to set the initial step size
    double step = Magnitude(fPars.Errors());

    // Usually we reuse the chisq from the last iteration, so need to
    // initialize it once out here.
    double chi = fFunc(pt);
    int ncalls = 1;

    while(true){
      //      std::cout << "New direction, chisq = " << chi << std::endl;

      // Evaluate gradient to figure out a direction to step in
      std::vector<double> grad = fFunc.Gradient(pt);
      ++ncalls;
      // Make gradient into a unit vector but keep the magnitude around
      const double gradmag = Magnitude(grad);
      MakeUnit(grad);

      // Take step in that direction
      std::vector<double> trialpt = pt;
      for(unsigned int i = 0; i < N; ++i) trialpt[i] -= grad[i]*step;

      // And evaluate the chisq there
      const double trialchi = fFunc(trialpt);
      ++ncalls;

      //      std::cout << "step = " << step << " -> chisq " << trialchi << std::endl;

      // Estimate the second derivative from the two points and one
      // gradient
      const double d2 = (trialchi+gradmag*step-chi)/(step*step);
      // We expect the function to be convex
      if(d2 > 0){
        //        std::cout << " d2 = " << d2;
        // If so, we can compute a better step size, one that would place us
        // right at the minimum if the function only had quadratic terms.
        step = gradmag/(2*d2);
        //        std::cout << " new step = " << step << std::endl;
      }

      while(true){
        // Keep trying steps until we find one that reduces the chisq
        std::vector<double> newpt = pt;
        for(unsigned int i = 0; i < N; ++i) newpt[i] -= grad[i]*step;

        const double newchi = fFunc(newpt);
        ++ncalls;
        //        std::cout << "  step = " << step << " -> chisq = " << newchi << std::endl;

        if(newchi > chi){
          // If the chisq went up instead, try again with smaller step
          step /= 2;
          if(step < 1e-8){
            //            std::cout << "Step too small!" << std::endl;
            // Maybe it's just because we're already in the minimum. If we
            // can't improve even with very tiny steps, call it a day.
            return Package(pt, chi, ncalls);
          }
          continue;
        }
        else{
          if(chi-newchi < 1e-5){
            //            std::cout << "Small chisq step" << std::endl;
            // If the improvement we did get is really tiny we're also likely
            // already at the minimum
            return Package(newpt, newchi, ncalls);
          }

          // In all other cases (ie we took a reasonable step and found a
          // reasonably better chisq) we want to update our state, take another
          // look at the gradient vector to figure out which direction to go
          // next, and preserve our step size, which is some kind of good
          // estimate of the scale of the function.
          pt = newpt;
          chi = newchi;
          break;
        }
      } // end line search
    } // end while (searching for better pt)
  }

  //----------------------------------------------------------------------
  ROOT::Minuit2::FunctionMinimum GradientDescent::
  Package(const std::vector<double>& pt, double chi, int ncalls) const
  {
    // In practice we only use the chisq and the parameter values, so don't be
    // too careful about setting this all up.
    const unsigned int N = pt.size();
    ROOT::Minuit2::MnAlgebraicVector vec(N);
    for(unsigned int i = 0; i < N; ++i) vec(i) = pt[i];
    ROOT::Minuit2::MinimumParameters params(vec, chi);    
    ROOT::Minuit2::MinimumState state(params, 0, ncalls);
    ROOT::Minuit2::MnUserTransformation trans(pt, std::vector<double>(N));
    ROOT::Minuit2::MinimumSeed seed(state, trans);
    return ROOT::Minuit2::FunctionMinimum(seed, {state}, chi);
  }

  //----------------------------------------------------------------------
  double GradientDescent::Magnitude(const std::vector<double>& xs) const
  {
    double ret = 0;
    for(double x: xs) ret += x*x;
    return sqrt(ret);
  }

  //----------------------------------------------------------------------
  void GradientDescent::MakeUnit(std::vector<double>& xs) const
  {
    const double mag = Magnitude(xs);
    for(double& x: xs) x /= mag;
  }
}
