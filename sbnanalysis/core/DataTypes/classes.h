#include "sbncode/sbnanalysis/core/DataTypes/Event.hh"
#include "sbncode/sbnanalysis/core/DataTypes/Experiment.hh"
#include "sbncode/sbnanalysis/core/DataTypes/SubRun.hh"

template class std::vector<sbnanalysis::Event>;
template class std::vector<sbnanalysis::SubRun>;
template class std::vector<sbnanalysis::Event::Metadata>;
template class std::vector<sbnanalysis::Event::Interaction>;
template class std::vector<sbnanalysis::Event::Neutrino>;
template class std::vector<sbnanalysis::Event::FinalStateParticle>;
template class std::map<std::string, std::vector<double> >;
template class std::vector<std::map<std::string,std::vector<double> > >;

