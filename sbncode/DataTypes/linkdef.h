#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link C++ nestedclasses;

#pragma link C++ class event::Event+;
#pragma link C++ class event::RecoEvent+;
#pragma link C++ class SubRun+;
#pragma link C++ class FileMeta+;
#pragma link C++ class event::Metadata+;
#pragma link C++ class event::Interaction+;
#pragma link C++ class event::RecoInteraction+;
#pragma link C++ class event::Neutrino+;
#pragma link C++ class event::FinalStateParticle+;
#pragma link C++ class std::map<std::string, std::vector<float> >+;
#pragma link C++ class std::vector<event::Interaction>+;
#pragma link C++ class std::vector<event::RecoInteraction>+;
#pragma link C++ class std::vector<event::FinalStateParticle>+;

#pragma link C++ class std::vector<TVector3>+;

#pragma link C++ class map<string,vector<float> >+;
#pragma link C++ class vector<map<string,vector<float> > >+;

#endif

