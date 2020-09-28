////////////////////////////////////////////////////////////////////////
// \file    SRTrackCalo.h
////////////////////////////////////////////////////////////////////////
#ifndef SRTRACKCALO_H
#define SRTRACKCALO_H

namespace caf
{
  /// Track PID from dE/dx v. residual range Chi2
  class SRTrackCalo
    {
    public:

      SRTrackCalo();
      virtual ~SRTrackCalo();

      int            nhit;        //!< Number of hits on this plane counted in the calorimetry
      float          ke;          //!< Kinetic energy deposited on this plane [GeV]
      float          charge;      //!< Deposited charge as seen by wireplane (pre recombination and electric lifetime corrections) [#electrons]
      void setDefault();
    };

} // end namespace

#endif // SRTRACKCALO_H
//////////////////////////////////////////////////////////////////////////////
