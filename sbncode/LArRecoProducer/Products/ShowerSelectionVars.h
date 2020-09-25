#ifndef sbncode_Showerselectionvars_H
#define sbncode_Showerselectionvars_H

namespace sbn {
  class ShowerDensityFit {
    public:

      ShowerDensityFit(double densityGrad=-999, double densityPow=-999):
        mDensityGrad(densityGrad),
        mDensityPow(densityPow)
    {};

      double mDensityGrad;
      double mDensityPow;
  }; // end sbnd::ShowerDensityFit

  class ShowerTrackFit {
    public:

      ShowerTrackFit(double trackLength=-999, double trackWidth=-999, unsigned int numHits=0):
        mTrackLength(trackLength),
        mTrackWidth(trackWidth),
        mNumHits(numHits)
    {}
      double mTrackLength;
      double mTrackWidth;
      unsigned int mNumHits;
  }; // end sbnd::ShowerTrackFit
}

#endif
