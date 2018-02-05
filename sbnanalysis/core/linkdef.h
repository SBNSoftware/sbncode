#ifdef __CINT__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link C++ nestedclasses;

#pragma link C++ class Event+;
#pragma link C++ class Event::Metadata+;
#pragma link C++ class Event::Interaction+;
#pragma link C++ class Event::Neutrino+;
#pragma link C++ class Event::FinalStateParticle+;
#pragma link C++ class std::vector<Event::Interaction>+;
#pragma link C++ class std::vector<Event::FinalStateParticle>+;
#pragma link C++ class std::map<std::string, std::vector<double> >+;

#pragma link C++ class TVector3+;
#pragma link C++ class std::vector<TVector3>+;
#pragma link C++ class std::vector<TLorentzVector>+;

#endif

