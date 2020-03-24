#include <iostream>
#include "SBNAna/Cuts/NumuCuts.h"
#include "StandardRecord/Proxy/SRProxy.h"
#include <cassert>

namespace ana{

  const Cut kNumuBasicQual([](const caf::SRProxy* sr)
			   { // 
			     bool hastrk = ( sr->reco.ntrk > 0 );
			     if (!hastrk) return hastrk;
			 
			     unsigned int muIdx = (unsigned int)kPrimMuonIdx(sr);
 			     double len = sr->reco.trk[muIdx].len;
			     return (len > 0 );
			   });

  const Cut kHasFlashMatch([](const caf::SRProxy* sr)
			{ 
			  return ( sr->slc.fmatch.present );
			});

  const Cut kFlashMatchScore([](const caf::SRProxy* sr)
			{ 
			  return ( sr->slc.fmatch.time > 0 && 
				   sr->slc.fmatch.score > 6   );
			});

  const Cut kFlashMatchNumuCut = kHasFlashMatch  && kFlashMatchScore;

  
  const Cut kNumuTrkLen([](const caf::SRProxy* sr)
			{ // TO DO: Prim pandora trk len primary_track.length > fConfig.TrackLength
			  bool hastrk = ( sr->reco.ntrk > 0 );
			  if (!hastrk) return hastrk;
			  
			  unsigned int muIdx = (unsigned int)kPrimMuonIdx(sr);
			  double len = sr->reco.trk[muIdx].len;
			  return (len > 50 );
			});
  
  // const Cut kNumuContND([](const caf::SRProxy* sr)
  // 			{ // TO DO: Trk end-point containment
  // 			  bool hastrk = ( sr->reco.ntrk > 0 );
  // 			  if (!hastrk) return hastrk;
  // 			  // CosmicContain from sbnanalysis OscRecoPostprocess
  // 			  //  ytop: 40	 ybottom: 15  zfront: 15   zback: 15
  // 			  unsigned int muIdx = (unsigned int)kPrimMuonIdx(sr);
  // 			  return ( sr->reco.trk[muIdx].start.x > ( kNDTopX -15 ) &&
  // 				   sr->reco.trk[muIdx].start.x < 15 &&
  // 				   sr->reco.trk[muIdx].start.y > ( kNDTopY -40 ) &&
  // 				   sr->reco.trk[muIdx].start.y < 15 &&
  // 				   sr->reco.trk[muIdx].start.z > ( kNDBack -15 ) &&
  // 				   sr->reco.trk[muIdx].start.z < 15 &&
  // 				   sr->reco.trk[muIdx].end.x > ( kNDTopX -15 ) &&
  // 				   sr->reco.trk[muIdx].end.x < 15 &&
  // 				   sr->reco.trk[muIdx].end.y > ( kNDTopY -40 ) &&
  // 				   sr->reco.trk[muIdx].end.y < 15 &&
  // 				   sr->reco.trk[muIdx].end.z > ( kNDBack -15 ) &&
  // 				   sr->reco.trk[muIdx].end.z < 15 && )
  // 			});
  
  // const Cut kNumuContFD([](const caf::SRProxy* sr)
  // 			{ // TO DO: Trk end-point containment
  // 			  bool hastrk = ( sr->reco.ntrk > 0 );
  // 			  return ( hastrk );
  // 			});
  
  // const Cut kNumuVtxFidND([](const caf::SRProxy* sr)
  // 			  { // TO DO: Vtx Containment
  // 			  bool hastrk = ( sr->reco.ntrk > 0 );
  // 			    return ( hastrk );

  // 			  });
  
  // const Cut kNumuVtxFidFD([](const caf::SRProxy* sr)
  // 			  { // TO DO: Vtx Containment
  // 			  bool hastrk = ( sr->reco.ntrk > 0 );
  // 			    return ( hastrk );
  // 			  });  
    
  // const Cut kNumuMCSLen([](const caf::SRProxy* sr)
  // 			{ // TO DO: Prim trk MCS length -  primary_track.mcs_momentum < 7. /*garbage value*/&& 			  (fConfig.MCSTrackLength <0. || 			   primary_track.length > fConfig.MCSTrackLength
  // 			  bool hastrk = ( sr->reco.ntrk > 0 );
  // 			  return ( hastrk );
  // 			});
  
  
  // const Cut kNumuCRTHit([](const caf::SRProxy* sr)
  // 			{ // TO DO
  // 			  bool hastrk = ( sr->reco.ntrk > 0 );
  // 			  return ( hastrk );
  // 			});
  
  // const Cut kNumuCRTTrk([](const caf::SRProxy* sr)
  // 			{ // TO DO
  // 			  bool hastrk = ( sr->reco.ntrk > 0 );
  // 			  return ( hastrk );
  // 			});
  
  // const Cut kNumuCRTAct([](const caf::SRProxy* sr)
  // 			{ // TO DO
  // 			  bool hastrk = ( sr->reco.ntrk > 0 );
  // 			  return ( hastrk );
  // 			});  

}
