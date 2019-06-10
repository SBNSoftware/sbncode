/**
 * \file GetFV.cc
 *
 *
 * Author:
 */

#include <iostream>
#include <array>

#include "canvas/Utilities/InputTag.h"
#include "core/SelectionBase.hh"
#include "core/Event.hh"
#include "core/Experiment.hh"

#include "TH1D.h"
#include "TDatabasePDG.h"
#include "TGraph.h"

#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/MCBase/MCTrack.h"

#include "larcorealg/Geometry/TPCGeo.h"
#include "larcorealg/Geometry/CryostatGeo.h"
#include "larcorealg/Geometry/GeometryCore.h"

#include "core/ProviderManager.hh"

namespace ana {
  namespace SBNOsc {

/**
 * \class GetFV
 * \brief Electron neutrino event selection
 */
class GetFV : public core::SelectionBase {
public:
  /** Constructor. */
  GetFV() {}

  /**
   * Initialization.
   *
   * \param config A configuration, as a FHiCL ParameterSet object
   */
  void Initialize(fhicl::ParameterSet* config=NULL) {
    _managerSBND = new core::ProviderManager(kExpSBND);
    _managerUBOONE = new core::ProviderManager(kExpMicroBooNE);
    _managerICARUS = new core::ProviderManager(kExpICARUS);

    {
    // print out SBND
    int cryo_i = 0;
    int tpc_i = 0;
    std::cout << std::endl << "SBND" << std::endl;
    for (auto const &cryo: _managerSBND->GetGeometryProvider()->IterateCryostats()) { 
      cryo_i ++;
      geo::GeometryCore::TPC_iterator iTPC = _managerSBND->GetGeometryProvider()->begin_TPC(cryo.ID()),
                                      tend = _managerSBND->GetGeometryProvider()->end_TPC(cryo.ID());
      while (iTPC != tend) {
        geo::TPCGeo const& TPC = *iTPC;
        auto this_volume = TPC.ActiveBoundingBox();
        std::cout << "    { " << std::endl;
        std::cout << "        xmin: " << this_volume.MinX() + 8.25 << std::endl;
        std::cout << "        xmax: " << this_volume.MaxX() - 8.25 << std::endl;
        std::cout << "        ymin: " << this_volume.MinY() + 15 << std::endl;
        std::cout << "        ymax: " << this_volume.MaxY() - 15 << std::endl;
        std::cout << "        zmin: " << this_volume.MinZ() + 15 << std::endl;
        std::cout << "        zmax: " << this_volume.MaxZ() - 80 << std::endl;
        std::cout << "    }, " << std::endl;
        tpc_i ++;
        iTPC ++;
      }
      tpc_i = 0;
    }
    }
   
    {
    // print out uBooNE
    int cryo_i = 0;
    int tpc_i = 0;
    std::cout << std::endl << "uBooNE" << std::endl;
    for (auto const &cryo: _managerUBOONE->GetGeometryProvider()->IterateCryostats()) { 
      cryo_i ++;
      geo::GeometryCore::TPC_iterator iTPC = _managerUBOONE->GetGeometryProvider()->begin_TPC(cryo.ID()),
                                      tend = _managerUBOONE->GetGeometryProvider()->end_TPC(cryo.ID());
      while (iTPC != tend) {
        geo::TPCGeo const& TPC = *iTPC;
        auto this_volume = TPC.ActiveBoundingBox();
        std::cout << "    { " << std::endl;
        std::cout << "        xmin: " << this_volume.MinX() + 15 << std::endl;
        std::cout << "        xmax: " << this_volume.MaxX() - 15 << std::endl;
        std::cout << "        ymin: " << this_volume.MinY() + 15 << std::endl;
        std::cout << "        ymax: " << this_volume.MaxY() - 15 << std::endl;
        std::cout << "        zmin: " << this_volume.MinZ() + 15 << std::endl;
        std::cout << "        zmax: " << this_volume.MaxZ() - 80 << std::endl;
        std::cout << "    }, " << std::endl;
        tpc_i ++;
        iTPC ++;
      }
      tpc_i = 0;
    }
    }

    {
    // print out ICARUS
    int cryo_i = 0;
    int tpc_i = 0;
    std::cout << std::endl << "ICARUS" << std::endl;
    std::cout << "NCHANNELS: " << _managerICARUS->GetGeometryProvider()->Nchannels() << std::endl;
    for (auto const &cryo: _managerICARUS->GetGeometryProvider()->IterateCryostats()) { 
      cryo_i ++;
      geo::GeometryCore::TPC_iterator iTPC = _managerICARUS->GetGeometryProvider()->begin_TPC(cryo.ID()),
                                      tend = _managerICARUS->GetGeometryProvider()->end_TPC(cryo.ID());
      while (iTPC != tend) {
        geo::TPCGeo const& TPC = *iTPC;
        auto this_volume = TPC.ActiveBoundingBox();
        std::cout << "    { " << std::endl;
        std::cout << "        xmin: " << this_volume.MinX() + 8.25 << std::endl;
        std::cout << "        xmax: " << this_volume.MaxX() - 8.25 << std::endl;
        std::cout << "        ymin: " << this_volume.MinY() + 15 << std::endl;
        std::cout << "        ymax: " << this_volume.MaxY() - 15 << std::endl;
        std::cout << "        zmin: " << this_volume.MinZ() + 15 << std::endl;
        std::cout << "        zmax: " << this_volume.MaxZ() - 80 << std::endl;
        std::cout << "    }, " << std::endl;
        tpc_i ++;
        iTPC ++;
      }
      tpc_i = 0;
    }
  }
  }

  /** Finalize and write objects to the output file. */
  void Finalize() {}

  /**
   * Process one event.
   *
   * \param ev A single event, as a gallery::Event
   * \param Reconstructed interactions
   * \return True to keep event
   */
  bool ProcessEvent(const gallery::Event& ev, const std::vector<Event::Interaction> &truth, std::vector<Event::RecoInteraction>& reco) {return false; }

protected:
  core::ProviderManager *_managerICARUS;
  core::ProviderManager *_managerSBND;
  core::ProviderManager *_managerUBOONE;
};

  }  // namespace SBNOsc
}  // namespace ana
DECLARE_SBN_PROCESSOR(ana::SBNOsc::GetFV)

