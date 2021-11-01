#include <map>
#include <string>
#include "canvas/Persistency/Common/Wrapper.h"

#include "dk2nu/tree/dk2nu.h"
#include "dk2nu/tree/NuChoice.h"
#include "nusimdata/SimulationBase/MCTruth.h"

template class art::Wrapper<art::Ptr<bsim::Dk2Nu> >;
template class art::Wrapper<art::Ptr<bsim::NuChoice> >;
template class art::Wrapper<art::Ptr<simb::MCTruth> >;

