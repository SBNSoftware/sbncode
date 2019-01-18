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
#include "TMultiGraph.h"

#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/MCBase/MCTrack.h"

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
  }; 

  class TrackProjection {
    public:
    std::vector<TGraph *> tracks; //!< TGraph of each individual graph
    std::unique_ptr<TMultiGraph> graph; //!< graph of all the tracks
    std::string name;

    TrackProjection(std::string name):
      graph(new TMultiGraph(name.c_str(), name.c_str())),
      name(name)
      {}

    void Draw(Config &config, const char *draw_str, unsigned axis0, unsigned axis1) const {
      TCanvas c(name.c_str(), name.c_str(), 600, 400);
      c.cd();
      graph->Draw(draw_str);
      c.Update();
      graph->GetHistogram()->GetXaxis()->SetLimits(aaBoxesMin(config.active_volumes, axis0), aaBoxesMax(config.active_volumes, axis0));
      graph->GetHistogram()->GetYaxis()->SetRangeUser(aaBoxesMin(config.active_volumes, axis1), aaBoxesMax(config.active_volumes, axis1));
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
  };


  class ProjectionGraphs {
    public:
    std::array<TGraph *, 3> graphs; 
    ProjectionGraphs(std::vector<TVector3> &points, TrackPlot *plots) {
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
      TGraph *xy = new TGraph(xs.size(), &xs[0], &ys[0]);
      TGraph *xz = new TGraph(xs.size(), &xs[0], &zs[0]);
      TGraph *yz = new TGraph(xs.size(), &ys[0], &zs[0]);

      // add to the multigraph -- which now takes ownership of the TGraph
      if (plots != NULL) {
        plots->xy.graph->Add(xy, "");
        plots->xz.graph->Add(xz, "");
        plots->yz.graph->Add(yz, "");
        plots->xy.tracks.push_back(xy);
        plots->xz.tracks.push_back(xz);
        plots->yz.tracks.push_back(yz);
      }

      // save them for the user if necessary
      graphs = {xy, xz, yz};
    }
  };
/*
  class WireProjectionGraphs {
    public:
      std::vector<std::vector<TGraph *>> graphs;
      WireProjectionGraphs(std::vector<recob::Hit> &hits, unsigned n_planes, TrackPlot *plots):
        graphs({}, n_planes)
      {
        // data for each plane
        std::vector<std::vector<double>> times({}, n_planes);
        std::vector<std::vector<double>> channels({}, n_planes);
        for (auto const &hit: hits) {
          // get plane index
          unsigned plane = hit.WireID().planeID();
          // get channel and time
          double channel = hit.Channel();
          double time = hit.PeakTime();
          // store
          channels[plane].push_back(channel);
          times[plane].push_back(time);
        }

        // save all the plots
        for (unsigned i = 0; i < times.size(); i++) {
          TGraph *this_graph = new TGraph(channels[i].size(), &channels[0], &times[0]);
          if (plots != NULL) {
            plots->planes[i].Add(this_graph, "");
          }
          // for later if necessary
          graphs[i].push_back(this_graph);
        }
      }
  };
*/

  std::vector<TrackPlot> _event_plots;
  unsigned _event_index;
  Config _config;

};

  }  // namespace SBNOsc
}  // namespace ana

#endif  // __sbnanalysis_ana_SBNOsc_NumuSelection__

