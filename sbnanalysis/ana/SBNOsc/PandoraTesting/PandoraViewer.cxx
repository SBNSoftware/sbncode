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

#include "core/Event.hh"
#include "PandoraViewer.h"
#include "Utilities.h"

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
    // drawing options
    _config.draw_str = config->get<std::string>("draw_str", "A");

    // get tag names
    _config.HitTag = config->get<std::string>("HitTag", "gaushit");
    _config.RecoTrackTag = config->get<std::string>("RecoTrackTag", "pandoraTrack");
  }

  _event_index = 0;
}

void PandoraViewer::Finalize() {
  fOutputFile->cd();

  const char *draw_str = _config.draw_str.c_str();

  std::cout << "x: " << aaBoxesMin(_config.active_volumes, 0) << " " << aaBoxesMax(_config.active_volumes, 0) << std::endl;

  // save the multigraphs
  // save everything as a canvas
  // alost now set the range
  for (auto const &plots: _event_plots) {
    plots.xy.Draw(_config, draw_str, 0, 1);
    plots.xz.Draw(_config, draw_str, 0, 2);
    plots.yz.Draw(_config, draw_str, 1, 2);
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

  // get the hits
  auto const& gaus_hits = \
    *ev.getValidHandle<std::vector<recob::Hit>>(_config.HitTag);

  // get the reco tracks
  auto const &reco_tracks = \
    *ev.getValidHandle<std::vector<recob::Track>>(_config.RecoTrackTag);

  // make the plots
  TrackPlot plots(_event_index);

  // iterate through the reco tracks
  for (auto const &track: reco_tracks) {
    // make a tgraph for all the reco tracks
    std::vector<TVector3> points;
    for (unsigned i = 0; i < track.NumberTrajectoryPoints(); i++) {
      auto pos = track.TrajectoryPoint(i).position;

      // TODO: why are some poitns at -999?
      if (pos.X() < -900 && pos.Y() < -900 &&pos.Z() < -900) continue;

      TVector3 vect(pos.X(), pos.Y(), pos.Z());
      points.push_back(vect);
    }
    ProjectionGraphs graphs(points, &plots);

    // format the graphs
    for (int i = 0; i < 3; i++) {
      graphs.graphs[i]->SetLineColor(kRed);
      graphs.graphs[i]->SetLineWidth(2.);
      graphs.graphs[i]->SetLineStyle(2);
    }
  }

  // Iterate through the truth tracks
  for (auto const &track: mctracks) {
    // make a tgraph of trajectory points for each track

    // ignore zero-sized tracks
    if (track.size() == 0) continue;

    // get points
    std::vector<TVector3> points;
    for (auto const &trajPoint: track) {
      TVector3 vect = trajPoint.Position().Vect();
      points.push_back(vect);
    }
    // save them in the multigraph
    ProjectionGraphs graphs(points, &plots);

    for (int i = 0; i < 3; i++) {
      graphs.graphs[i]->SetLineColor(kBlack);
      graphs.graphs[i]->SetLineWidth(2.);
      graphs.graphs[i]->SetLineStyle(3);
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

