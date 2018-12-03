#include "nusimdata/SimulationBase/MCTruth.h"
#include "sbncode/NuMuSelection/DataTypes/AssocTruthInfo.h"
#include "canvas/Persistency/Common/Wrapper.h"
#include "canvas/Persistency/Common/Assns.h"
#include <vector>

namespace {
  struct dictionary {
    numuselection::AssocTruthInfo ato;
    std::vector<numuselection::AssocTruthInfo> v_ato;
    art::Wrapper<numuselection::AssocTruthInfo> w_ato;
    art::Wrapper<std::vector<numuselection::AssocTruthInfo>> w_v_ato;
    art::Assns<numuselection::AssocTruthInfo,simb::MCTruth, void> a_ato_mct;
    art::Wrapper<art::Assns<numuselection::AssocTruthInfo,simb::MCTruth, void>> w_a_ato_mct;
    art::Assns<simb::MCTruth,numuselection::AssocTruthInfo, void> a_mct_ato;
    art::Wrapper<art::Assns<simb::MCTruth,numuselection::AssocTruthInfo, void>> w_a_mct_ato;
  };
}
