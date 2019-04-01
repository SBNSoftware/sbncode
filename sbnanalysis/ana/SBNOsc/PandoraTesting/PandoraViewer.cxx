#include <iostream>
#include <cmath>
#include <vector>

#include "fhiclcpp/ParameterSet.h"

#include <TH2D.h>
#include <TH1D.h>
#include <TCanvas.h>

#include "gallery/ValidHandle.h"
#include "gallery/Handle.h"
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

bool PandoraViewer::containedInAV(TVector3& point) {
  for (auto const &av: _config.active_volumes) {
    if (av.Contain(point)) return true;
  }
  return false;
}

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
    _config.CRTTrackTag = config->get<std::string>("CRTTrackTag", "crttrack");

  }

  _event_index = 0;
}

void PandoraViewer::Finalize() {
  fOutputFile->cd();

  // save the multigraphs
  // save everything as a canvas
  // alost now set the range
  for (auto const &plots: _event_plots) {
    plots.xy.Draw(_config, 0, 1);
    plots.xz.Draw(_config, 0, 2);
    plots.yz.Draw(_config, 1, 2);
    for (auto const &graph: plots.tpc_planes) {
      graph.DrawNoLimit(_config);
    }
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
  TrackPlot plots(_event_index, fProviderManager);

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

  // get the CRT tracks
  gallery::Handle<std::vector<sbnd::crt::CRTTrack>> crt_tracks;
  ev.getByLabel<std::vector<sbnd::crt::CRTTrack>>(_config.CRTTrackTag, crt_tracks);

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
      if (points.size() == 0) continue;

      std::array<TGraph *, 3> graphs = plots.addTrack(points);
      std::vector<TGraph *> tpc_graphs = plots.addTrackTPC(points);
      std::vector<TGraph *> all_graphs(std::begin(graphs), std::end(graphs)); 
      all_graphs.insert(all_graphs.end(), tpc_graphs.begin(), tpc_graphs.end());

      // save points
      std::array<TText *, 3> text = plots.addTrackText(points, track->ID());
      for (int i = 0; i < 3; i++) {
        // text[i]->SetTextSize(1);
        text[i]->SetTextColor(kRed);
      }

      // format the graphs
      for (int i = 0; i < all_graphs.size(); i++) {
        if (is_neutrino) {
          all_graphs[i]->SetLineColor(kGreen);
          all_graphs[i]->SetLineWidth(4.);
        }
        else if (is_cosmic) {
          all_graphs[i]->SetLineColor(kBlue);
          all_graphs[i]->SetLineWidth(2.);
        }
        else {
          all_graphs[i]->SetLineColor(kRed);
          all_graphs[i]->SetLineWidth(2.);
        }

        all_graphs[i]->SetLineStyle(2);
      }
    }

    // draw the reco vertices
    for (auto const &vertex: this_vertices) {
      auto pos = vertex->position();
      TVector3 rv(pos.X(), pos.Y(), pos.Z());
      std::array<TGraph *, 3> graphs = plots.addVertex(rv);
      std::vector<TGraph *> tpc_graphs = plots.addVertexTPC(rv);

      std::vector<TGraph *> all_graphs(std::begin(graphs), std::end(graphs)); 
      all_graphs.insert(all_graphs.end(), tpc_graphs.begin(), tpc_graphs.end());
      for (int i = 0; i < all_graphs.size(); i++) {
        if (is_neutrino) {
          all_graphs[i]->SetMarkerColor(kGreen);
        }
        else if (is_cosmic) {
          all_graphs[i]->SetMarkerColor(kBlue);
        }
        else {
          all_graphs[i]->SetMarkerColor(kRed);
        }
        all_graphs[i]->SetMarkerSize(3.);
        all_graphs[i]->SetMarkerStyle(5);
      }


    }
  }

  // Iterate through the CRT Tracks
  if (crt_tracks.isValid()) {
    for (auto const &track: *crt_tracks) {
      // get the start and end 
      TVector3 p1(track.x1_pos, track.y1_pos, track.z1_pos);
      TVector3 p2(track.x2_pos, track.y2_pos, track.z2_pos);
      
      std::vector<TVector3> points {p2, p1};
      
      // get the errors
      TVector3 e1(track.x1_err, track.y1_err, track.z1_err);
      TVector3 e2(track.x2_err, track.y2_err, track.z2_err);
      std::vector<TVector3> errors {e2, e1};
      
      std::array<TGraphErrors *, 3> graphs = plots.addTrackwErrors(points, errors);
      for (int i = 0; i < 3; i++) {
	graphs[i]->SetLineColor(kGray);
	graphs[i]->SetLineWidth(2.);
	graphs[i]->SetFillColor(kGray);
	graphs[i]->SetFillStyle(3002);
      }
    }
  }

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
      if (containedInAV(vect)) {
        points.push_back(vect);
      }
    }
    if (points.size() == 0) continue;

    // save them in the multigraph
    std::array<TGraph *, 3> graphs = plots.addTrack(points);
    std::vector<TGraph *> tpc_graphs = plots.addTrackTPC(points);
    std::vector<TGraph *> all_graphs(std::begin(graphs), std::end(graphs)); 
    all_graphs.insert(all_graphs.end(), tpc_graphs.begin(), tpc_graphs.end());

    // save text for the truth tracks
    std::array<TText *, 3> text = plots.addTrackText(points, track.TrackID());
    for (int i = 0; i < 3; i++) {
      // text[i]->SetTextSize(1);
    }


    // format the graphs
    for (int i = 0; i < all_graphs.size(); i++) {
      all_graphs[i]->SetLineColor(kBlack);
      all_graphs[i]->SetLineWidth(2.);
      all_graphs[i]->SetLineStyle(1);
    }
  }

  // now do the vertices
  for (auto const &mctruth: mctruths) {
    if (!mctruth.NeutrinoSet()) continue;

    TVector3 vertex = mctruth.GetNeutrino().Nu().Trajectory().Position(0).Vect();
    std::array<TGraph *, 3> graphs = plots.addVertex(vertex);
    std::vector<TGraph *> tpc_graphs = plots.addVertexTPC(vertex);

    std::vector<TGraph *> all_graphs(std::begin(graphs), std::end(graphs)); 
    all_graphs.insert(all_graphs.end(), tpc_graphs.begin(), tpc_graphs.end());
    for (int i = 0; i < all_graphs.size(); i++) {
      all_graphs[i]->SetMarkerColor(kBlack);
      all_graphs[i]->SetMarkerSize(3.);
      all_graphs[i]->SetMarkerStyle(4);
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

