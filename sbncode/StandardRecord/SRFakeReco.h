//  SRFakeReco.h
// \author  grayputnam@uchicago.edu
////////////////////////////////////////////////////////////////////////
#ifndef SRFAKERECO_H
#define SRFAKERECO_H

#include "sbncode/StandardRecord/SRVector3D.h"
#include "sbncode/StandardRecord/SRFakeRecoParticle.h"

namespace caf
{
  /// The SRFakeReco is a faked reconstruction using estimates from the SBN proposal 
  class SRFakeReco
  {
  public:
    SRFakeReco();
    ~SRFakeReco() {  }

    float nuE; ///! Fake-reco neutrino Energy [GeV]
    SRVector3D vtx; ///! Interaction vertex in detector coordinates [cm] 
    SRFakeRecoParticle lepton; ///! Fake-reco lepton information 
    std::vector<SRFakeRecoParticle> hadrons; ///! Fake-reco information on hadronic state
    int nhad; ///! Number of hadrons
    float wgt; ///! Weight for this interaction
  };

} // end namespace

#endif // SRFAKERECO_H
//////////////////////////////////////////////////////////////////////////////
