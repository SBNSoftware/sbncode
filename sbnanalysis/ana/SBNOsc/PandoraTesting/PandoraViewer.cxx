#include <iostream>
#include <cmath>
#include <vector>

#include "fhiclcpp/ParameterSet.h"

#include <TH2D.h>
#include <TH1D.h>
#include <TCanvas.h>

#include "gallery/ValidHandle.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/FindMany.h"

#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"

#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"

#include "core/Event.hh"
#include "PandoraViewer.h"

#include "../Utilities.h"

#include "sbndcode/CRT/CRTUtils/CRTAnaUtils.h"
#include "sbndcode/CRT/CRTProducts/CRTHit.hh"

namespace ana {
  namespace SBNOsc {


void PandoraViewer::Initialize(fhicl::ParameterSet* config) {
  if (config) {
    // setup active volume bounding boxes
    std::vector<fhicl::ParameterSet> AVs = \
      config->get<std::vector<fhicl::ParameterSet> >("active_volumes");
    for (auto const& AV : AVs) {
      double xmin = AV.get<double>("xmin");
      double ymin = AV.get<double>("ymin");
      double zmin = AV.get<double>("zmin");
      double xmax = AV.get<double>("xmax");
      double ymax = AV.get<double>("ymax");
      double zmax = AV.get<double>("zmax");
      _config.active_volumes.emplace_back(xmin, ymin, zmin, xmax, ymax, zmax);
    }

    // get tag names
    _config.HitTag = config->get<std::string>("HitTag", "gaushit");
    _config.RecoTrackTag = config->get<std::string>("RecoTrackTag", "pandoraTrack");
    _config.RecoVertexTag = config->get<std::string>("RecoVertexTag", "pandora");

    // setup services manager
    // _services_manager = new art::ServicesManager(config->get<fhicl::ParameterSet>("services"), _registry);
  }

  _event_index = 0;
}

void PandoraViewer::Finalize() {
  fOutputFile->cd();

  std::cout << "x: " << aaBoxesMin(_config.active_volumes, 0) << " " << aaBoxesMax(_config.active_volumes, 0) << std::endl;

  // save the multigraphs
  // save everything as a canvas
  // alost now set the range
  for (auto const &plots: _event_plots) {
    plots.xy.Draw(_config, 0, 1);
    plots.xz.Draw(_config, 0, 2);
    plots.yz.Draw(_config, 1, 2);
  }

}

bool PandoraViewer::ProcessEvent(const gallery::Event& ev, const std::vector<Event::Interaction> &truth, std::vector<Event::RecoInteraction>& reco) {
  // Get truth
  auto const& mctruths = \
    *ev.getValidHandle<std::vector<simb::MCTruth> >(fTruthTag);
  // get tracks and showers
  auto const& mctracks = \
    *ev.getValidHandle<std::vector<sim::MCTrack> >(fMCTrackTag);
  auto const& mcshowers = \
    *ev.getValidHandle<std::vector<sim::MCShower> >(fMCShowerTag);

  // make the plots
  TrackPlot plots(_event_index);

  // get the PFParticles
  auto const &pfp_handle = \
    ev.getValidHandle<std::vector<recob::PFParticle>>("pandora");

  // use these to get all the vertices
  art::FindMany<recob::Vertex> pfp_vertices(pfp_handle, ev, "pandora");
  // and the tracks
  art::FindMany<recob::Track> pfp_tracks(pfp_handle, ev, "pandoraTrack");
  // and the metadata
  art::FindMany<larpandoraobj::PFParticleMetadata> pfp_metadatas(pfp_handle, ev, "pandora");

  // get the CRT hits
  gallery::ValidHandle<std::vector<sbnd::crt::CRTHit>> crt_hits = \
    ev.getValidHandle<std::vector<sbnd::crt::CRTHit>>("crthit");

  // build the CRT Tracks
  // std::vector<std::vector<sbnd::crt::CRTHit>> crt_zeros = sbnd::CRTAnaUtils::CreateCRTTzeros(*crt_hits, 0.2);
  // std::vector<sbnd::crt::CRTTrack> crt_tracks = sbnd::CRTAnaUtils::CreateCRTTracks(std::move(crt_zeros), 
  //   30., false, 25.);

  // iterate over all of the pfparticles
  for (size_t i = 0; i < pfp_handle->size(); i++) {
    // get the PFParticle
    const recob::PFParticle& this_pfp = pfp_handle->at(i);
    // get the metadata
    const larpandoraobj::PFParticleMetadata* this_metadata = pfp_metadatas.at(i).at(0);
    // and the properties dict
    auto const &properties = this_metadata->GetPropertiesMap();
    // get the reco vertices
    const std::vector<const recob::Vertex*> &this_vertices = pfp_vertices.at(i);
    // and the tracks 
    const std::vector<const recob::Track*> &this_tracks = pfp_tracks.at(i);

    // figure out if neutrino
    bool is_neutrino = false;
    if (properties.count("IsNeutrino")) {
      is_neutrino = properties.at("IsNeutrino");
    }

    // and if cosmic
    bool is_cosmic = false;
    if (properties.count("IsClearCosmic")) {
      is_cosmic = properties.at("IsClearCosmic");
    }

    // iterate through the reco tracks
    for (auto const &track: this_tracks) {
      // make a tgraph for all the reco tracks
      std::vector<TVector3> points;
      for (unsigned i = 0; i < track->NumberTrajectoryPoints(); i++) {
        auto pos = track->TrajectoryPoint(i).position;

        // TODO: why are some poitns at -999?
        if (pos.X() < -900 && pos.Y() < -900 &&pos.Z() < -900) continue;

        TVector3 vect(pos.X(), pos.Y(), pos.Z());
        points.push_back(vect);
      }
      std::array<TGraph *, 3> graphs = plots.addTrack(points);

      // format the graphs
      for (int i = 0; i < 3; i++) {
        if (is_neutrino) {
          graphs[i]->SetLineColor(kGreen);
          graphs[i]->SetLineWidth(4.);
        }
        else if (is_cosmic) {
          graphs[i]->SetLineColor(kBlue);
          graphs[i]->SetLineWidth(2.);
        }
        else {
          graphs[i]->SetLineColor(kRed);
          graphs[i]->SetLineWidth(2.);
        }

        graphs[i]->SetLineStyle(2);
      }
    }

    // draw the reco vertices
    for (auto const &vertex: this_vertices) {
      auto pos = vertex->position();
      TVector3 rv(pos.X(), pos.Y(), pos.Z());
      std::array<TGraph *, 3> graphs = plots.addVertex(rv);
      for (int i = 0; i < 3; i++) {
        if (is_neutrino) {
          graphs[i]->SetMarkerColor(kGreen);
        }
        else if (is_cosmic) {
          graphs[i]->SetMarkerColor(kBlue);
        }
        else {
          graphs[i]->SetMarkerColor(kRed);
        }
        graphs[i]->SetMarkerSize(3.);
        graphs[i]->SetMarkerStyle(5);
      }
    }
  }

  // Iterate through the CRT Tracks
  /*for (auto const &track: crt_tracks) {
    std::cout << std::endl << "New track" << std::endl;
    // make a tgraph of trajectory points for each track
    std::cout << "x1: " << track.x1_pos<< std::endl;
    std::cout << "y1: " << track.y1_pos<< std::endl;
    std::cout << "z1: " << track.z1_pos<< std::endl;

    std::cout << "x2: " << track.x2_pos<< std::endl;
    std::cout << "y2: " << track.y2_pos<< std::endl;
    std::cout << "z2: " << track.z2_pos<< std::endl;

    std::cout << "length: " << track.length << std::endl;

    std::cout << "time 0: " << track.ts0_s << " [s] " << track.ts0_ns << " [ns]" << std::endl;
    std::cout << "time 1: " << track.ts0_s << " [s] " << track.ts1_ns << " [ns]" << std::endl;

  }*/

  // Iterate through the truth tracks
  for (auto const &track: mctracks) {
    // make a tgraph of trajectory points for each track

    // ignore zero-sized tracks
    if (track.size() == 0) continue;
    // ignore neutral tracks
    if (PDGCharge(track.PdgCode()) < 1e-4) continue;

    // get points
    std::vector<TVector3> points;
    for (auto const &trajPoint: track) {
      TVector3 vect = trajPoint.Position().Vect();
      points.push_back(vect);
    }
    // save them in the multigraph
    std::array<TGraph *, 3> graphs = plots.addTrack(points);

    for (int i = 0; i < 3; i++) {
      graphs[i]->SetLineColor(kBlack);
      graphs[i]->SetLineWidth(2.);
      graphs[i]->SetLineStyle(3);
    }
  }

  // now do the vertices
  for (auto const &mctruth: mctruths) {
    TVector3 vertex = mctruth.GetNeutrino().Nu().Trajectory().Position(0).Vect();
    std::array<TGraph *, 3> graphs = plots.addVertex(vertex);
    for (int i = 0; i < 3; i++) {
      graphs[i]->SetMarkerColor(kBlack);
      graphs[i]->SetMarkerSize(3.);
      graphs[i]->SetMarkerStyle(4);
    }
  }

  // keep track of the plots if there were any tracks
  if (plots.xy.graph->GetListOfGraphs() != NULL && plots.xy.graph->GetListOfGraphs()->GetSize() > 0) {
    _event_plots.push_back(std::move(plots));
  }
  
  _event_index ++;

  return true; // always select
}

  }  // namespace SBNOsc
}  // namespace ana


DECLARE_SBN_PROCESSOR(ana::SBNOsc::PandoraViewer)

