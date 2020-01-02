#include "Flatten.h"
#include "SetEvent.h"

#include "TClassTable.h"
#include "TClass.h"
#include "TDictionary.h"
#include "TDictAttributeMap.h"
#include "TProtoClass.h"
#include "TSystem.h"
#include "TParameter.h"

#include "../NumuReco/TrackAlgo.h"

std::vector<std::string> ProcessClass(TClass *tclass, const std::string &prefix) {
  std::vector<std::string> ret;
  TList *members = tclass->GetListOfDataMembers();
  TIterator *m_iterator = members->MakeIterator();
  TObject *obj;
  while ((obj = m_iterator->Next()) != NULL) {
    TDataMember *member = (TDataMember *)obj;
    TDataType *basic_type = member->GetDataType();
    if (basic_type != NULL) {
      // only accept floats
      assert(basic_type->GetType() == kFloat_t);
      // check if array
      if (member->GetArrayDim() == 0) {
        // not array -- just add
        ret.push_back(prefix + "." + std::string(obj->GetName()));
      }
      else {
        // sort through array
        unsigned max_ind = 1;
        for (unsigned i = 0; i < member->GetArrayDim(); i++) {
          max_ind = max_ind * member->GetMaxIndex(i);
        }
        for (unsigned i = 0; i < max_ind; i++) {
          std::vector<unsigned> each_ind(member->GetArrayDim(), 0);
          unsigned ind = i;
          for (int j = member->GetArrayDim()-1; j >= 0; j--) {
            unsigned this_i = ind % member->GetMaxIndex(j);
            each_ind[j] = this_i;
            ind = ind / member->GetMaxIndex(j);
          }
          std::string base = prefix + "." + std::string(obj->GetName()); 
          for (unsigned ind: each_ind) base += ("." + std::to_string(ind));
          ret.push_back(base);
        }
      } 

    }
    // no enums in TNtuple
    assert(!member->IsEnum());
    // no pointers in TNtuple
    assert(!member->IsaPointer());
    // no containers in TNtuple
    assert(!member->IsSTLContainer());
    // if not basic -- descend
    if (!member->IsBasic()) {
      DictFuncPtr_t next = gClassTable->GetDict(member->GetTypeName());
      if (next == NULL) {
        std::cout << "Missing dictionary for class type: " << member->GetTypeName() << std::endl;
        assert(next != NULL);
      }
      std::vector<std::string> this_ret = ProcessClass(next(), prefix + "." + std::string(obj->GetName()));
      ret.insert(ret.end(), this_ret.begin(), this_ret.end());
    }
  }
  return ret;
}

std::string ClasstoNTupleName(const std::string &classname, const std::string &prefix) {
  DictFuncPtr_t classdict = gClassTable->GetDict(classname.c_str());
  std::vector<std::string> var_names = ProcessClass(classdict(), prefix);
  std::string ret;
  for (unsigned i = 0; i < var_names.size(); i++) {
    ret += var_names[i];
    if (i < var_names.size() -1) ret += ":";
  }
  return ret;
}

void ana::SBNOsc::Flatten::Initialize(fhicl::ParameterSet *config) {
  fCuts.Initialize(config->get<fhicl::ParameterSet>("Cuts"), fProviderManager->GetGeometryProvider());
  std::string output_name = config->get<std::string>("OutputFile", "output.root");

  fOutputFiles = std::vector<TFile *>(NWorkers());
  fNtuples = std::vector<TNtuple *>(NWorkers());
  fInteractions = std::vector<numu::flat::FlatInteraction>(NWorkers());
  fRecoEvents = std::vector<numu::RecoEvent *>(NWorkers(), 0);
  fMCTypes = std::vector<numu::MCType>(NWorkers());
  std::string ntuple = ClasstoNTupleName("numu::flat::FlatInteraction", "interaction");

  for (unsigned i = 0; i < NWorkers(); i++) {
    std::string name = output_name;
    if (i > 0) name += std::to_string(i);
    fOutputFiles[i] = new TFile(name.c_str(), "RECREATE");
    fOutputFiles[i]->cd();
    fNtuples[i] = new TNtuple("interaction", "interaction", ntuple.c_str());
  }
}

void ana::SBNOsc::Flatten::InitializeThread() {
  unsigned i = WorkerID();
  fOutputFiles[i]->cd();
}

void ana::SBNOsc::Flatten::FileSetup(TFile *f, TTree *eventTree) {
  unsigned ind = WorkerID();

  eventTree->SetBranchAddress("reco_event", &fRecoEvents[ind]);

  TParameter<int> *file_type = (TParameter<int> *)f->Get("MCType");
  fMCTypes[ind] = (numu::MCType) file_type->GetVal();  

  fInteractions[ind].meta.mc_type =  fMCTypes[ind];
  fInteractions[ind].meta.detector = fConfigExperimentID; 
}

void ana::SBNOsc::Flatten::ProcessEvent(const event::Event *ev) {
  unsigned index = WorkerID();
  ana::SBNOsc::SetEvent(*fRecoEvents[index], *ev, fCuts, fMCTypes[index], true);
  for (const numu::RecoInteraction &interaction: fRecoEvents[index]->reco) {
    // set stuff
    // 
    // Primary Track
    numu::flat::PrimaryTrack track;
    track.length = interaction.primary_track.length; 
    track.costh = interaction.primary_track.costh;
    track.range_momentum = numu::RangeMomentum(interaction.primary_track);
    track.mcs_momentum = numu::MCSMomentum(interaction.primary_track);
    track.crt_hit_distance = (interaction.primary_track.crt_match.hit_match.present) ? interaction.primary_track.crt_match.hit_match.distance : -1;
    track.crt_hit_time = interaction.primary_track.crt_match.hit_match.time;
    track.crt_track_angle = (interaction.primary_track.crt_match.track.present) ? interaction.primary_track.crt_match.track.angle: -1;
    track.crt_track_time = interaction.primary_track.crt_match.track.time;
    interaction.primary_track.start.GetXYZ(track.start);
    interaction.primary_track.end.GetXYZ(track.end);
    if (interaction.primary_track.match.has_match) {
      const numu::TrueParticle &truth = fRecoEvents[index]->particles.at(interaction.primary_track.match.mcparticle_id);
      truth.start_momentum.GetXYZ(track.truth.momentum);
      track.truth.wall_enter = truth.wall_enter;
      track.truth.wall_exit = truth.wall_exit;
      track.truth.time = truth.start_time;
      track.truth.is_cosmic = truth.is_cosmic;
      track.truth.pdgid = truth.pdgid;
      track.truth.is_contained = truth.is_contained;
    }
    else track.truth.pdgid = -1; 
   
    // Slice
    numu::flat::Slice slice;
    slice.flash_score = interaction.slice.flash_match.score;
    slice.flash_time = interaction.slice.flash_match.time;
  
    // True Neutrino
    numu::flat::TrueNeutrino neutrino;
    if (interaction.match.has_match) {
      neutrino.E = ev->truth[interaction.match.mctruth_track_id].neutrino.energy;
      neutrino.Q2 = ev->truth[interaction.match.mctruth_track_id].neutrino.Q2;
      ev->truth[interaction.match.mctruth_track_id].neutrino.position.GetXYZ(neutrino.vertex);
    }
    else neutrino.E = -1;

    // Event Info
    numu::flat::EventInfo event_info;
    for (unsigned i = 0; i < 10; i++) event_info.crt_hit_pes[i] = -1;
    // Event Info
    unsigned i = 0;
    for (const numu::CRTHit &hit: fRecoEvents[index]->in_time_crt_hits) {
      assert(i < 10);
      event_info.crt_hit_pes[i] = hit.pes; 
      event_info.crt_hit_times[i] = hit.time;
      i += 1;
    }
    event_info.pass_trig = fCuts.PassFlashTrigger(*fRecoEvents[index]);

    // Event Meta data -- set previously

    fNtuples[index]->Fill((float *)&fInteractions[index]);

    // overwrite normalization stuff so it isn't duplicated
    fInteractions[index].meta.pot = 0;
    fInteractions[index].meta.n_gen_events = 0;
  }

}

void ana::SBNOsc::Flatten::ProcessSubRun(const SubRun *subrun) {
  unsigned index = WorkerID();
  fInteractions[index].meta.pot += subrun->totgoodpot; 
}

void ana::SBNOsc::Flatten::ProcessFileMeta(const FileMeta *meta) {
  unsigned index = WorkerID();
  fInteractions[index].meta.n_gen_events += meta->n_gen_events;
}

void ana::SBNOsc::Flatten::Finalize() {
  TList *merge = new TList;
  fOutputFiles[0]->cd();
  for (unsigned i = 1; i < NWorkers(); i++) {
    merge->Add(fNtuples[i]);
  }
  fNtuples[0]->Merge(merge);
  fNtuples[0]->Write();
  delete merge;

  for (unsigned i = 1; i < NWorkers(); i++) {
    const char *name = fOutputFiles[i]->GetName();
    fOutputFiles[i]->Close();
    gSystem->Unlink(name);
  }
 
}



DECLARE_SBN_POSTPROCESSOR(ana::SBNOsc::Flatten);
