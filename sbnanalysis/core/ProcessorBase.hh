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
#include <map>
#include "lardataobj/RecoBase/MCSFitResult.h"
#include "larreco/RecoAlg/TrajectoryMCSFitter.h"
#include "larreco/RecoAlg/TrackMomentumCalculator.h"
#include "gallery/Event.h"
#include "Loader.hh"
#include "Event.hh"

class TBranch;
class TFile;
class TTree;
class Event;
class SubRun;

namespace fhicl {
  class ParameterSet;
}

/** Core framework functionality. */
namespace core {

class ProviderManager;

/**
 * \class core::ProcessorBase
 * \brief A generic tree-writing event-by-event processor.
 */
class ProcessorBase {
friend class ProcessorBlock;
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
   * \param ev The event, as a gallery::Event
   * \param reco Reco interactions, to be populated by the user
   * \returns True if event passes filter
   */
  virtual bool ProcessEvent(const gallery::Event& ev,
                            const std::vector<Event::Interaction> &truth,
                            std::vector<Event::RecoInteraction>& reco) = 0;

  /** Pointer to reco event information */
  std::vector<Event::RecoInteraction>* fReco;  //!< Reco interaction list

protected:
  /**
   * Perform user-level initialization.
   *
   * \param config A configuration, as a JSON filename.
   */
  virtual void Initialize(char* config=NULL);

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
   * \param config A configuration as a JSON filename.
   */
  virtual void Setup(char* config=NULL);

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
   * \param ev The current gallery event
  */
  void BuildEventTree(gallery::Event& ev);


  /**
   * Setup LArSoft services
   *
   * \param ev The current gallery event
  */
  void SetupServices(gallery::Event& ev);

  /**
   * Update subrun list to include subruns for this event's file.
   *
   * \param ev The current gallery event
   */
  void UpdateSubRuns(gallery::Event& ev);

  unsigned long fEventIndex;                          //!< An incrementing index
  Experiment fExperimentID;                           //!< Experiment identifier
  ProviderManager* fProviderManager;                  //!< Interface for provider access
  std::string fOutputFilename;                        //!< The output filename
  std::string fProviderConfig;                        //!< A custom provider config fcl file
  TFile* fOutputFile;                                 //!< The output ROOT file
  TTree* fTree;                                       //!< The output ROOT tree
  Event* fEvent;                                      //!< The standard output event data structure
  TTree* fSubRunTree;                                 //!< Subrun output tree
  SubRun* fSubRun;                                    //!< Standard output subrun structure
  std::set<std::pair<int, int> > fSubRunCache;        //!< Cache stored subruns
  art::InputTag fFluxTag;                             //!< art tag for MCTruth information
  art::InputTag fTruthTag;                            //!< art tag for MCTruth information
  std::vector<art::InputTag> fWeightTags;             //!< art tag(s) for MCEventWeight information
  art::InputTag fMCTrackTag;                          //!< art tag for MCTrack
  art::InputTag fMCShowerTag;                         //!< art tag for MCShower
  art::InputTag fMCParticleTag;                       //!< art tag for MCParticle
  art::InputTag fHitTag;                              //!< art tag for Hits
  art::InputTag fPFParticleTag;                       //!< art tag for PFParticles
  art::InputTag fRecoTrackTag;                        //!< art tag for Reconstructed tracks  
  art::InputTag fRecoShowerTag;                       //!< art tag for Reconstructed showers  
  art::InputTag fVertexTag;                           //!< art tag for Vertices
  art::InputTag fRecoTrackCalorimetryTag;             //!< art tag for Reconstructed track calorimetry
  art::InputTag fRecoTrackParticleIDTag;              //!< art tag for Reconstructed track pid
  trkf::TrajectoryMCSFitter::Config fMCSFitterConfig; //!< MCS momentum fitter config
  trkf::TrajectoryMCSFitter *fMCSFitter;              //!< MCS momentum fitter
  trkf::TrackMomentumCalculator fRangeFitter;         //!< Range momentum fitter
  int fMCSFitterPID;
  int fMCSFitterNumSegments;
  double fMCSFitterLengthSegment;
  int fMCSFitterHitsSegment;
  int fMCSFitterElossSteps;   
  int fMCSFitterElossMode;    
  double fMCSFitterPMin;        
  double fMCSFitterPMax;      
  double fMCSFitterPStep;          
  double fMCSFitterAngResol;
};

}  // namespace core


/** Macro to create plugin library for user-defined Processors. */
#define DECLARE_SBN_PROCESSOR(classname) extern "C" { \
core::ProcessorBase* CreateProcessorObject() { return new classname; } \
void DestroyProcessorObject(core::ProcessorBase* o) { delete o; } \
struct core::export_table exports = { CreateProcessorObject, DestroyProcessorObject };}

#endif  // __sbnanalysis_core_ProcessorBase__

