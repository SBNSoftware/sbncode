#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link C++ nestedclasses;

#pragma link C++ class numu::RecoEvent+;
#pragma link C++ class numu::RecoInteraction+;
#pragma link C++ class numu::RecoParticle+;
#pragma link C++ class numu::TrueParticle+;
#pragma link C++ class numu::RecoTrack+;
#pragma link C++ class numu::MCSFitResult+;
#pragma link C++ class numu::SliceTruth+;
#pragma link C++ class numu::TrackTruth+;
#pragma link C++ class std::vector<numu::RecoInteraction>+;
#pragma link C++ class std::vector<numu::CRTHit>+;
#pragma link C++ class std::map<size_t, numu::RecoParticle>+;
#pragma link C++ class std::map<size_t, numu::RecoTrack>+;
#pragma link C++ class numu::RecoSlice+;
#pragma link C++ class numu::CRTMatch+;
#pragma link C++ class numu::CRTMatch::Track+;
#pragma link C++ class numu::CRTMatch::HitMatch+;
#pragma link C++ class numu::FlashMatch+;
#pragma link C++ class std::optional<numu::CRTMatch>+;
#pragma link C++ class std::optional<numu::FlashMatch>+;
#pragma link C++ class numu::CRTHit;

#pragma link C++ class numu::flat::FlatInteraction;
#pragma link C++ class numu::flat::PrimaryTrack;
#pragma link C++ class numu::flat::Slice;
#pragma link C++ class numu::flat::TrueNeutrino;
#pragma link C++ class numu::flat::TrackTruth;
#pragma link C++ class numu::flat::EventInfo;
#pragma link C++ class numu::flat::EventMeta;
#pragma link C++ class std::vector<float>+;
#pragma link C++ class std::vector<std::vector<float>>+;
#endif
