#ifndef __sbnanalysis_ana_AutoVetoAnalysis_AutoVetoSelection__
#define __sbnanalysis_ana_AutoVetoAnalysis_AutoVetoSelection__

/**
 * \file AutoVetoSelection.h
 *
 * An example event selection processor.
 *
 * This is an implementation of the core::SelectionBase class. We define
 * the methods called for initialization, finalization, and event-by-event
 * processing.
 *
 * Author: A. Mastbaum <mastbaum@uchicago.edu>
 */

// Includes
#include <iostream>
#include "canvas/Utilities/InputTag.h"
#include "core/SelectionBase.hh"
#include "core/Event.hh"

// Forward declarations
class TH2D;
class TH1D;
class TH1F;

/** All analysis code is defined in namespace "ana" */
namespace ana {

  /** Code specific to the AutoVetoAnalysis. */
  namespace AutoVetoAnalysis {

/**
 * \class AutoVetoSelection
 * \brief An example selection analysis
 *
 * This selection analysis doesn't actually select events, it just
 * demonstrates the framework!
 */
class AutoVetoSelection : public core::SelectionBase {
public:
  /** Constructor. */
  AutoVetoSelection();

  /**
   * Initialization.
   *
   * Here we load configuration parameters, set up histograms for output, and
   * add our own branches to the output tree.
   *
   * \param config A configuration, as FHiCL ParameterSet object
   */
  void Initialize(fhicl::ParameterSet* config=NULL);

  /** Finalize and write objects to the output file. */
  void Finalize();

  /**
   * Process one event.
   *
   * \param ev A single event, as a gallery::Event
   * \param reco Reconstructed interactions
   * \return True to keep event
   **/
  bool ProcessEvent(const gallery::Event& ev,
                    const std::vector<event::Interaction>& truth,
                    std::vector<event::RecoInteraction>& reco);

protected:
  unsigned fEventCounter;  //!< Count processed events

  /** Configuration parameters */
  art::InputTag fTruthTag;  //!< art tag for MCTruth information
  art::InputTag fDetTag;    //!< art tag for DetSim information
  double fFidXOut;          //!< fiducial cut in drift direction from outside TPC boundary[cm]
  double fFidXIn;           //!< fiducial cut in drift direction from cathode [cm]
  double fFidYTop;          //!< fiducial cut in vertical direction from top boundary [cm]
  double fFidYBot;          //!< fiducial cut in vertical direction from bottom boundary [cm]
  double fFidZUp;           //!< fiducial cut in beam direction from upstream boundary [cm]
  double fFidZDown;         //!< fiducial cut in beam direction from downstream boundary [cm]

  /** Custom data branches */
  int fNuCount;             //!< Number of neutrino interactions in the event, vertex anywhere in cryostat
  int fCRTTag;              //!< Number of neutrino events tagged by CRT, vertex anywhere in cryostat
  int fCCNC;                //!< Whether nuetrino interaction was CC(0) or NC(1)
  int fNuReg;               //!< Location of neutrino vertex in LAr: LAr inactive volume (0), LAr AV(1), LAr FV(2)

  /** Histograms */
  TH2D* fNuVertexXZHist;    //!< Neutrino vertex XZ projection, anywhere in cryostat
  TH2D* fNuVertexXYHist;    //!< Neutrino vertex XY projection, anywhere in cryostat
  TH2D* fNuVertexYZHist;    //!< Neutrino vertex YZ projection, anywhere in cryostat
  TH2D* fNuVertexXZAVHist;  //!< Neutrino vertex XZ projection, anywhere in LAr AV
  TH2D* fNuVertexXYAVHist;  //!< Neutrino vertex XY projection, anywhere in LAr AV
  TH2D* fNuVertexYZAVHist;  //!< Neutrino vertex YZ projection, anywhere in LAr AV
  TH2D* fNuVertexXZFVHist;  //!< Neutrino vertex XZ projection, anywhere in LAr FV
  TH2D* fNuVertexXYFVHist;  //!< Neutrino vertex XY projection, anywhere in LAr FV
  TH2D* fNuVertexYZFVHist;  //!< Neutrino vertex YZ projection, anywhere in LAr FV

  TH1D* fNuE;               //!< Neutrino energy, vertex anywhere in cryostat
  TH1D* fNuEAV;             //!< Neutrino energy, vertex anywhere in LAr AV
  TH1D* fNuEAVNotFV;        //!< Neutrino energy, vertex anywhere in LAr AV except FV
  TH1D* fNuEFV;             //!< Neutrino energy, vertex anywhere in LAr FV
  TH1D* fNuEVeto;           //!< Auto-vetoed neutrino energy, vertex anywhere in cryostat
  TH1D* fNuEVetoAV;         //!< Auto-vetoed neutrino energy, vertex anywhere in LAr AV
  TH1D* fNuEVetoAVNotFV;    //!< Auto-vetoed neutrino energy, vertex anywhere in LAr AV except FV
  TH1D* fNuEVetoFV;         //!< Auto-vetoed neutrino energy, vertex anywhere in LAr FV

  TH1F* fCRTReg;            //!< First CRT region hit, vertex anywhere in cryostat
  TH1F* fCRTRegAV;          //!< First CRT region hit, vertex anywhere in LAr AV
  TH1F* fCRTRegAVNotFV;     //!< First CRT region hit, vertex anywhere in LAr AV except FV
  TH1F* fCRTRegFV;          //!< First CRT region hit, vertex anywhere in LAr FV
  TH1F* fNCRTReg;           //!< Number of different CRT regions hit, vertex anywhere in cryostat
  TH1F* fNCRTRegAV;         //!< Number of different CRT regions hit, vertex anywhere in LAr AV
  TH1F* fNCRTRegAVNotFV;    //!< Number of different CRT regions hit, vertex anywhere in LAr AV, except FV
  TH1F* fNCRTRegFV;         //!< Number of different CRT regions hit, vertex anywhere in LAr FV
};

  }  // namespace AutoVetoAnalysis
}  // namespace ana

#endif  // __sbnanalysis_ana_ExampleAnalysis_ExampleSelection__

