#include "SBNAna/Vars/NumuVars.h"
#include "CAFAna/Core/Utilities.h"
#include "StandardRecord/Proxy/SRProxy.h"
#include <cassert>

namespace ana
{


  const Var kPrimMuonIdx([](const caf::SRProxy *sr)
                        {       //Find the most muon-like track			  
			  if( (int)sr->reco.ntrk == 0 ) return -5.0;

                          double best_idx   = 0;
                          //double best_score = -5.0;
                          double best_len   = -5.0;			  
			  for( unsigned int trkIdx = 0; trkIdx < sr->reco.ntrk; trkIdx++ ){
			    auto &trk = sr->reco.trk[trkIdx];
			    if( (int)trk.chi2pid.size() < 1 ) return -5.0;
			    
			    //Find longest trk w/Chi2 for muon < Chi2 for pion
			    bool isMuonLike = trk.nchi2pid > 0 && trk.chi2pid[2].chi2_pion > trk.chi2pid[2].chi2_muon;
			    if( isMuonLike && trk.len > best_len ){ //Chi2 muon < Chi2 pion
			      // best_score = score;
			      best_len = trk.len;
			      best_idx = trkIdx;
			    }
			  } // end loop over trks
			  
                          return best_idx;
                        });


  const Var kPrimTrkLen(
			[](const caf::SRProxy *sr)
			{
			  double len = -5.0;
			  if ( sr->reco.ntrk > 0 ){
			    int muIdx = (int)kPrimMuonIdx(sr);
                            if (muIdx >= 0) { 
  			      len = sr->reco.trk[muIdx].len;
                            }
			  }
			  return len;
			});

}
