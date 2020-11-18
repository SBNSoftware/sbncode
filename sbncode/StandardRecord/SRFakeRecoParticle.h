//  SRFakeRecoParticle.h
// \author  grayputnam@uchicago.edu
////////////////////////////////////////////////////////////////////////
#ifndef SRFAKERECOPARTICLE_H
#define SRFAKERECOPARTICLE_H

#include "sbncode/StandardRecord/SRTruthMatch.h"
#include "sbncode/StandardRecord/SRTrueParticle.h"

namespace caf
{
  /// The SRFakeRecoParticle is a faked reconstruction using estimates from the SBN proposal 
  class SRFakeRecoParticle
  {
  public:
    SRFakeRecoParticle();
    ~SRFakeRecoParticle() {  }

    float ke; ///! Fake-reco kinetic energy [GeV]
    float costh; ///! Fake-reco cosine of angle w.r.t. beam direction
    float len; ///! Fake-reco particle length [cm]
    int pid;    ///! Fake-reco particle ID
    bool contained; ///! Whether contained
  };

} // end namespace

#endif // SRFAKERECOPARTICLE_H
//////////////////////////////////////////////////////////////////////////////
