#ifndef __sbnanalysis_ana_SBNOsc_NumuSelection__
#define __sbnanalysis_ana_SBNOsc_NumuSelection__

/**
 * \file NumuSelection.h
 *
 * SBN nue selection.
 *
 * Author:
 */

#include <iostream>
#include <sstream>

#include "canvas/Utilities/InputTag.h"
#include "core/SelectionBase.hh"
#include "core/Event.hh"

#include "../Utilities.h"

#include "TH1D.h"
#include "TDatabasePDG.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"

#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/MCBase/MCTrack.h"
#include "sbndcode/CRT/CRTProducts/CRTTrack.hh"

// take the geobox stuff from uboonecode
#include "ubcore/LLBasicTool/GeoAlgo/GeoAABox.h"
#include "ubcore/LLBasicTool/GeoAlgo/GeoAlgo.h"

class TH2D;

namespace ana {
  namespace SBNOsc {

double aaBoxesMin(const std::vector<geoalgo::AABox> &boxes, unsigned dim) {
  return std::min_element(boxes.begin(), boxes.end(), [dim](auto &lhs, auto &rhs) { return lhs.Min()[dim] < rhs.Min()[dim]; })->Min()[dim];
}

double aaBoxesMax(const std::vector<geoalgo::AABox> &boxes, unsigned dim) {
  return std::max_element(boxes.begin(), boxes.end(), [dim](auto &lhs, auto &rhs) { return lhs.Max()[dim] < rhs.Max()[dim]; })->Max()[dim];
}


/**
 * \class PandoraViewer
 * \brief little module for making plots of pandora tracking
 */
class PandoraViewer : public core::SelectionBase {
public:
  // framework functions
  void Initialize(fhicl::ParameterSet* config);
  void Finalize();
  bool ProcessEvent(const gallery::Event& ev, const std::vector<Event::Interaction> &truth, std::vector<Event::RecoInteraction>& reco);

  class Config {
  public:
    std::vector<geoalgo::AABox> active_volumes;
    std::string draw_str;
    std::string HitTag;
    std::string RecoTrackTag;
    std::string RecoVertexTag;
    std::string CRTTrackTag;
  }; 

  class TrackProjection {
    public:
    std::vector<TGraph *> tracks; //!< TGraph of each individual graph
    std::unique_ptr<TMultiGraph> graph; //!< graph of all the tracks
    std::vector<TGraph *> vertices; //!< collection of vertices
    std::vector<TGraphErrors *> error_tracks; //!< tracks with errors
    std::string name;

    TrackProjection(std::string name):
      graph(new TMultiGraph(name.c_str(), name.c_str())),
      name(name)
      {}

    void Draw(Config &config, unsigned axis0, unsigned axis1) const {
      TCanvas c(name.c_str(), name.c_str(), 1200, 800);
      c.cd();
      // draw the multigraph thingy
      graph->Draw("ALP");
      c.Update();
      graph->GetHistogram()->GetXaxis()->SetLimits(aaBoxesMin(config.active_volumes, axis0), aaBoxesMax(config.active_volumes, axis0));
      graph->GetHistogram()->GetYaxis()->SetRangeUser(aaBoxesMin(config.active_volumes, axis1), aaBoxesMax(config.active_volumes, axis1));
      // draw all the vertex plots
      for (size_t i = 0; i < vertices.size(); i++) {
        vertices[i]->Draw("P");
      }
      // draw all the error graphs
      for (size_t i = 0; i < error_tracks.size(); i++) {
        error_tracks[i]->Draw("SAME 3");
      }
      c.Update();
      c.Write();
    }
  };

  class TrackPlot {
    public:
    TrackProjection xy;
    TrackProjection xz;
    TrackProjection yz;

    TrackPlot(unsigned event_index):
      xy("xy proj event: " + std::to_string(event_index)),
      xz("xz proj event: " + std::to_string(event_index)),
      yz("yz proj event: " + std::to_string(event_index)) {}

    std::array<TGraph *, 3> addVertex(TVector3 &point) {
      double x = point.X();
      double y = point.Y();
      double z = point.Z();

      // a graph with one point...
      TGraph *g_xy = new TGraph(1, &x, &y);
      TGraph *g_xz = new TGraph(1, &x, &z);
      TGraph *g_yz = new TGraph(1, &y, &z);

      xy.vertices.push_back(g_xy);
      xz.vertices.push_back(g_xz);
      yz.vertices.push_back(g_yz);

      return {g_xy, g_xz, g_yz};
    }

    std::array<TGraphErrors *, 3> addTrackwErrors(std::vector<TVector3> &points, std::vector<TVector3> &errors) {
      // first get each position in x, y, z
      std::vector<double> xs;
      std::vector<double> ys;
      std::vector<double> zs;

      for (TVector3 const &point: points) {
        xs.push_back(point.X());
        ys.push_back(point.Y());
        zs.push_back(point.Z());
      }

      // and the errors
      std::vector<double> ex;
      std::vector<double> ey;
      std::vector<double> ez;

      for (TVector3 const &err: errors) {
        ex.push_back(err.X());
        ey.push_back(err.Y());
        ez.push_back(err.Z());
      }

      // set up the tgraphs
      TGraphErrors *g_xy = new TGraphErrors(xs.size(), &xs[0], &ys[0], &ex[0], &ey[0]);
      TGraphErrors *g_xz = new TGraphErrors(xs.size(), &xs[0], &zs[0], &ex[0], &ez[0]);
      TGraphErrors *g_yz = new TGraphErrors(xs.size(), &ys[0], &zs[0], &ey[0], &ez[0]);

      // add to the multigraph -- which now takes ownership of the TGraph
      xy.error_tracks.push_back(g_xy);
      xz.error_tracks.push_back(g_xz);
      yz.error_tracks.push_back(g_yz);
      

      // save them for the user if necessary
      return {g_xy, g_xz, g_yz};
    }

    std::array<TGraph *, 3> addTrack(std::vector<TVector3> &points) {
      // first get each position in x, y, z
      std::vector<double> xs;
      std::vector<double> ys;
      std::vector<double> zs;

      for (TVector3 const &point: points) {
        xs.push_back(point.X());
        ys.push_back(point.Y());
        zs.push_back(point.Z());
      }

      // set up the tgraphs
      TGraph *g_xy = new TGraph(xs.size(), &xs[0], &ys[0]);
      TGraph *g_xz = new TGraph(xs.size(), &xs[0], &zs[0]);
      TGraph *g_yz = new TGraph(xs.size(), &ys[0], &zs[0]);

      // add to the multigraph -- which now takes ownership of the TGraph
      xy.graph->Add(g_xy, "");
      xz.graph->Add(g_xz, "");
      yz.graph->Add(g_yz, "");
      xy.tracks.push_back(g_xy);
      xz.tracks.push_back(g_xz);
      yz.tracks.push_back(g_yz);
      

      // save them for the user if necessary
      return {g_xy, g_xz, g_yz};
     }

  };

  std::vector<TrackPlot> _event_plots;
  unsigned _event_index;
  Config _config;

  //art::ServicesManager *_services_manager;
  //art::ActivityRegistry _registry;

};

  }  // namespace SBNOsc
}  // namespace ana

#endif  // __sbnanalysis_ana_SBNOsc_NumuSelection__

