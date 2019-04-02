#ifndef __sbnanalysis_core_ProcessorBase__
#define __sbnanalysis_core_ProcessorBase__

/**
 * \file ProcessorBase.hh
 *
 * A generic processor that writes an sbnanalysis tree.
 *
 * Author: A. Mastbaum <mastbaum@uchicago.edu>, 2018/01/25
 */

#include <string>
#include <vector>
#include "art/Framework/Principal/Event.h"
#include "sbncode/sbnanalysis/core/DataTypes/Event.hh"

class TBranch;
class TFile;
class TTree;
class Event;
class SubRun;

namespace fhicl {
  class ParameterSet;
}

/** Core framework functionality. */
namespace sbnanalysis {

/**
 * \class core::ProcessorBase
 * \brief A generic tree-writing event-by-event processor.
 */
class ProcessorBase {
public:
  /** Constructor */
  ProcessorBase();

  /** Destructor */
  virtual ~ProcessorBase();

  /**
   * Fill the tree and increment the event index.
   */
  virtual void FillTree();

  /**
   * Cleanup any objects that were filled per event
   */
  virtual void EventCleanup();

  /**
   * Add a branch to the output tree.
   *
   * Called in user subclasses to augment the default event tree.
   * This mirrors the TTree::Branch API.
   *
   * \param name The branch name
   * \param obj A pointer to the object
   * \returns A pointer to the created TBranch (we retain ownership)
   */
  template<class T>
  TBranch* AddBranch(std::string name, T* obj) {
    return fTree->Branch(name.c_str(), obj);
  }

  /**
   * Process one event.
   *
   * This also serves as a filter: if the function results false, it acts as a
   * filter and the event is not written out.
   *
   * \param ev The event, as a art::Event
   * \param reco Reco interactions, to be populated by the user
   * \returns True if event passes filter
   */
  virtual bool ProcessEvent(const art::Event& ev,
                            const std::vector<Event::Interaction> &truth,
                            std::vector<Event::RecoInteraction>& reco) = 0;

  /** Pointer to reco event information */
  std::vector<Event::RecoInteraction>* fReco;  //!< Reco interaction list

  /**
   * Perform user-level initialization.
   *
   * \param config A configuration, as a JSON object.
   */
  virtual void Initialize(fhicl::ParameterSet* config=NULL) = 0;

  /** Perform user-level finalization. */
  virtual void Finalize() = 0;

  /**
   * Perform framework-level initialization.
   *
   * \param config A configuration as a JSON object
   */
  virtual void Setup(fhicl::ParameterSet* config=NULL);

  /** Perform framework-level finalization. */
  virtual void Teardown();

  /**
   * Populate the default event tree variables.
   *
   * \param ev The current art event
  */
  void BuildEventTree(const art::Event& ev);

  /**
   * Update subrun list to include subruns for this event's file.
   *
   * \param ev The current art event
   */
  void UpdateSubRuns(const art::Event& ev);

  /** Get the event pointer. */
  Event* GetEvent() { return fEvent; }

protected:
  unsigned long fEventIndex;  //!< An incrementing index
  Experiment fExperimentID;  //!< Experiment identifier
  TFile* fOutputFile;  //!< The output ROOT file
  TTree* fTree;  //!< The output ROOT tree
  Event* fEvent;  //!< The standard output event data structure
  TTree* fSubRunTree;  //!< Subrun output tree
  SubRun* fSubRun;  //!< Standard output subrun structure
  std::set<std::pair<int, int> > fSubRunCache;  //!< Cache stored subruns
  art::InputTag fTruthTag;  //!< art tag for MCTruth information
  art::InputTag fFluxTag;  //!< art tag for MCFlux information
  art::InputTag fMCTrackTag; //!< art tag for MCTrack
  art::InputTag fMCShowerTag; //!< art tag for MCShower
  art::InputTag fMCParticleTag; //!< art tag for MCParticle
  std::vector<art::InputTag> fWeightTags;  //!< art tag(s) for MCEventWeight information
};

}  // namespace sbnanalysis

#endif  // __sbnanalysis_core_ProcessorBase__

