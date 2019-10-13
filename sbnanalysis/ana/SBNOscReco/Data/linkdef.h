#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link C++ nestedclasses;

#pragma link C++ class numu::RecoEvent+;
#pragma link C++ class numu::RecoInteraction+;
#pragma link C++ class numu::RecoParticle+;
#pragma link C++ class numu::RecoTrack+;
#pragma link C++ class numu::TruthMatch+;
#pragma link C++ class numu::TrackTruthMatch+;
#pragma link C++ class std::vector<numu::RecoInteraction>+;
#pragma link C++ class std::map<size_t, numu::RecoParticle>+;
#pragma link C++ class std::map<size_t, numu::RecoTrack>+;
#pragma link C++ class numu::RecoSlice+;
#pragma link C++ class numu::CRTMatch+;
#pragma link C++ class numu::CRTMatch::Track+;
#pragma link C++ class numu::CRTMatch::Hit+;
#pragma link C++ class numu::FlashMatch+;
#pragma link C++ class std::optional<numu::CRTMatch>+;
#pragma link C++ class std::optional<numu::FlashMatch>+;
#endif
