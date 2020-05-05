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
#include "lardataalg/DetectorInfo/DetectorPropertiesStandard.h"
#include "lardataalg/DetectorInfo/DetectorClocksStandard.h"

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
    _managerUBOONE = NULL;//new core::ProviderManager(kExpMicroBooNE);
    _managerICARUS = NULL;//new core::ProviderManager(kExpICARUS);

    {
    // print out SBND
    int cryo_i = 0;
    int tpc_i = 0;
    std::cout << std::endl << "SBND" << std::endl;
    std::cout << "DRIFT V: " << _managerSBND->GetDetectorPropertiesProvider()->DriftVelocity() << std::endl;
    std::cout << "E lifetime: " << _managerSBND->GetDetectorPropertiesProvider()->ElectronLifetime() << std::endl;
    std::cout << "Density: " << _managerSBND->GetDetectorPropertiesProvider()->Density() << std::endl;
    std::cout << "EField: " <<  _managerSBND->GetDetectorPropertiesProvider()->Efield() << std::endl;
    std::cout << "Sampling Rate: " <<  _managerSBND->GetDetectorPropertiesProvider()->SamplingRate() << std::endl;
    std::cout << "TPC tick period: " << _managerSBND->GetDetectorClocksProvider()->TPCClock().TickPeriod() << std::endl;
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

        for (auto const &plane: TPC.IteratePlanes()) {
          std::cout << "TPC: " << iTPC << " plane: " << plane.ID() << " view: " << plane.View() << " pitch: " << plane.WirePitch() << " sinphi: "<< plane.SinPhiZ() <<  " cosphi: " << plane.CosPhiZ() << std::endl;

          double pi = 3.14159265358979323846;
          std::cout << "Angle to vert: " <<  _managerSBND->GetGeometryProvider()->WireAngleToVertical( _managerSBND->GetGeometryProvider()->View(plane.ID()), plane.ID()) - 0.5*pi << std::endl;

          for (auto const &wire: plane.IterateWires()) {
            TVector3 dir = wire.GetEnd() - wire.GetStart();
            dir = dir.Unit();
            std::cout << "Wire direction x: " << dir.X() << " y: " << dir.Y() << " z: " << dir.Z() << std::endl;
            std::cout << "Wire Length: " << wire.Length() << std::endl;
            break;
          }
        }
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
    for (auto const &planeID: _managerICARUS->GetGeometryProvider()->IteratePlaneIDs()) {
      std::cout << "Plane ID: " << planeID.deepestIndex() << " type: " << _managerICARUS->GetGeometryProvider()->SignalType(planeID) << std::endl;
      for (auto const &wireID: _managerICARUS->GetGeometryProvider()->IterateWireIDs(planeID)) {
        std::cout << _managerICARUS->GetGeometryProvider()->PlaneWireToChannel(wireID) << " ";
      }
      std::cout << std::endl;
    }

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
  bool ProcessEvent(const gallery::Event& ev, const std::vector<event::Interaction> &truth, std::vector<event::RecoInteraction>& reco) {return false; }

protected:
  core::ProviderManager *_managerICARUS;
  core::ProviderManager *_managerSBND;
  core::ProviderManager *_managerUBOONE;
};

  }  // namespace SBNOsc
}  // namespace ana
DECLARE_SBN_PROCESSOR(ana::SBNOsc::GetFV)

