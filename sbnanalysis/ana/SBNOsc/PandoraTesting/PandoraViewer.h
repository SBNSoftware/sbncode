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

#include "lardataalg/DetectorInfo/DetectorClocks.h"
#include "lardataalg/DetectorInfo/DetectorClocksStandard.h"

#include "canvas/Utilities/InputTag.h"
#include "core/SelectionBase.hh"
#include "core/Event.hh"
#include "core/ProviderManager.hh"

#include "../Utilities.h"

#include "TH1D.h"
#include "TDatabasePDG.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TText.h"

#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/MCBase/MCTrack.h"
#include "sbndcode/CRT/CRTProducts/CRTTrack.hh"

// take the geobox stuff from uboonecode
#include "ubcore/LLBasicTool/GeoAlgo/GeoAABox.h"
#include "ubcore/LLBasicTool/GeoAlgo/GeoAlgo.h"

#include "larcorealg/Geometry/TPCGeo.h"
#include "larcorealg/Geometry/CryostatGeo.h"
#include "larcorealg/Geometry/GeometryCore.h"

#include "lardataobj/RecoBase/SpacePoint.h"

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
    bool drawText;
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
    std::vector<TGraph *> points; //!< collection of points
    std::vector<TGraphErrors *> error_tracks; //!< tracks with errors
    std::vector<TText *> text;
    std::string name;

    TrackProjection(std::string name):
      graph(new TMultiGraph(name.c_str(), name.c_str())),
      name(name)
      {}

    void DrawHits(Config &config, const geo::PlaneID &plane_id, core::ProviderManager *manager) const {
      // check if the graph has stuff
      if (graph->GetListOfGraphs() == NULL || graph->GetListOfGraphs()->GetEntries() == 0) return;
      TCanvas c(name.c_str(), name.c_str(), 1200, 800);
      c.cd();
      // draw the multigraph thingy
      graph->Draw("AP");
      c.Update();

      // set the limits
      // xlimit is the number of wires
      graph->GetHistogram()->GetXaxis()->SetLimits(0, manager->GetGeometryProvider()->Nwires(plane_id));
      // time limit is about 0->2000us
      graph->GetHistogram()->GetYaxis()->SetRangeUser(0., 2000.);

      /*
      // draw all the vertex plots
      for (size_t i = 0; i < vertices.size(); i++) {
        vertices[i]->Draw("P");
      }
      // draw all the error graphs
      for (size_t i = 0; i < error_tracks.size(); i++) {
        error_tracks[i]->Draw("SAME 3");
      }*/
      c.Update();
      c.Write();
    }

    void Draw(Config &config, unsigned axis0, unsigned axis1) const {
      TCanvas c(name.c_str(), name.c_str(), 1200, 800);
      c.cd();
      bool axes = false;
      TH1F *draw_histo = NULL;
      // check if the graph has stuff
      if (!(graph->GetListOfGraphs() == NULL || graph->GetListOfGraphs()->GetEntries() == 0)) {
        axes = true;
        graph->Draw("ALP");
        c.Update();
        draw_histo = graph->GetHistogram();
      }
      // draw all the vertex plots
      for (size_t i = 0; i < points.size(); i++) {
        if (axes) {
          points[i]->Draw("P");
        }
        else {
          points[i]->Draw("AP");
          axes = true;
          draw_histo = points[i]->GetHistogram();
        }
      }
      // draw all the error graphs
      for (size_t i = 0; i < error_tracks.size(); i++) {
        if (axes) {
          error_tracks[i]->Draw("SAME 3");
        }
        else {
          error_tracks[i]->Draw("A3");
          axes = true;
          draw_histo = error_tracks[i]->GetHistogram();
        }
      }

      if (draw_histo != NULL) {
        draw_histo->GetXaxis()->SetLimits(aaBoxesMin(config.active_volumes, axis0), aaBoxesMax(config.active_volumes, axis0));
        draw_histo->GetYaxis()->SetRangeUser(aaBoxesMin(config.active_volumes, axis1), aaBoxesMax(config.active_volumes, axis1));
      }

      // draw all the text
      for (size_t i = 0; i < text.size(); i++) {
        text[i]->Draw();
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
    std::vector<TrackProjection> plane_graphs;
    core::ProviderManager *_manager;
    std::vector<geo::PlaneID> _plane_ids;

    TrackPlot(unsigned event_index, core::ProviderManager *manager):
      xy("xy proj event: " + std::to_string(event_index)),
      xz("xz proj event: " + std::to_string(event_index)),
      yz("yz proj event: " + std::to_string(event_index)),
      _manager(manager)
    {
      for (const geo::PlaneID &plane_id: _manager->GetGeometryProvider()->IteratePlaneIDs()) {
        // get strring for induction or collection
        const char *sig_type = (_manager->GetGeometryProvider()->SignalType(plane_id) == geo::kCollection) ? "C " : "I ";
        plane_graphs.emplace_back(sig_type + plane_id.toString() + " event: " + std::to_string(event_index));
        _plane_ids.push_back(plane_id);

      }
    }

    std::array<TGraph *, 3> addVertex(TVector3 &point) {
      double x = point.X();
      double y = point.Y();
      double z = point.Z();

      // a graph with one point...
      TGraph *g_xy = new TGraph(1, &x, &y);
      TGraph *g_xz = new TGraph(1, &x, &z);
      TGraph *g_yz = new TGraph(1, &y, &z);

      xy.points.push_back(g_xy);
      xz.points.push_back(g_xz);
      yz.points.push_back(g_yz);

      return {g_xy, g_xz, g_yz};
    }

    /*
    std::vector<TGraph *> addVertexTPC(TVector3 &point) {
      std::vector<TGraph *> ret;
      // Add point for each TPC 
      unsigned graph_id = 0;
      for (auto const &cryo: _manager->GetGeometryProvider()->IterateCryostats()) {
        geo::GeometryCore::TPC_iterator iTPC = _manager->GetGeometryProvider()->begin_TPC(cryo.ID()),
                                        tend = _manager->GetGeometryProvider()->end_TPC(cryo.ID());
        while (iTPC != tend) {
          geo::TPCGeo const& TPC = *iTPC;
          if (TPC.ContainsPosition(point)) {
            // get the distance to the plane
            double height = TPC.DistanceFromReferencePlane(point);
            // and the wire in the plane
            int plane_no = 0;
            for (auto const &plane_id: _manager->GetGeometryProvider()->IteratePlaneIDs(TPC.ID())) {
              int wire_id;
              try {
                wire_id = _manager->GetGeometryProvider()->NearestWire(point, plane_id); //.deepestIndex();
              }
              catch (geo::InvalidWireError const& e) {
                continue;
              }
              //std::cout << " Height: " << height << " ID: " << (double) wire_id << std::endl;
              double xval = (double) wire_id;
              int this_graph_id = graph_id + plane_no;
              TGraph *graph = new TGraph(1, &xval, &height);
              tpc_plane[this_graph_id].vertices.push_back(graph);
              ret.push_back(graph);
              plane_no ++;
            }
          }
          graph_id += TPC.Nplanes();
          iTPC++;
        }
      }

      return ret;
    }*/

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

    std::array<TText *, 3> addTrackText(std::vector<TVector3> &points, int id) {
      if (points.size() == 0) return {NULL, NULL, NULL};
      // get the middle position
      int middle = points.size() / 2;
      double x = points[middle].X();
      double y = points[middle].Y();
      double z = points[middle].Z();

      std::string text = std::to_string(id);
      // make all the text
      TText *t_xy = new TText(x, y, text.c_str());
      TText *t_xz = new TText(x, z, text.c_str());
      TText *t_yz = new TText(y, z, text.c_str());

      xy.text.push_back(t_xy);
      xz.text.push_back(t_xz);
      yz.text.push_back(t_yz);
      return {t_xy, t_xz, t_yz};
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

    std::array<TGraph *, 3> addSpacePoints(const std::vector<art::Ptr<recob::SpacePoint>> &points) {
      // first get each position in x, y, z
      std::vector<double> xs;
      std::vector<double> ys;
      std::vector<double> zs;

      for (const art::Ptr<recob::SpacePoint> &point: points) {
        // get the xyz
        const double *xyz = point->XYZ();
        xs.push_back(xyz[0]);
        ys.push_back(xyz[1]);
        zs.push_back(xyz[2]);
      }

      // set up the tgraphs
      TGraph *g_xy = new TGraph(xs.size(), &xs[0], &ys[0]);
      TGraph *g_xz = new TGraph(xs.size(), &xs[0], &zs[0]);
      TGraph *g_yz = new TGraph(xs.size(), &ys[0], &zs[0]);

      // add to the multigraph -- which now takes ownership of the TGraph
      xy.graph->Add(g_xy, "");
      xz.graph->Add(g_xz, "");
      yz.graph->Add(g_yz, "");
      xy.points.push_back(g_xy);
      xz.points.push_back(g_xz);
      yz.points.push_back(g_yz);

      // save them for the user if necessary
      return {g_xy, g_xz, g_yz};

    }

    std::vector<TGraph *> addTrackHits(const std::vector<art::Ptr<recob::Hit>> &hits, core::ProviderManager *manager) {
      // setup to keep track of the points in each Plane
      std::vector<std::array<std::vector<double>, 2>> plane_graph_points;
      for (int i = 0; i < _plane_ids.size(); i++) {
        std::array<std::vector<double>, 2> this_graph { std::vector<double>(), std::vector<double>() };
        plane_graph_points.emplace_back(this_graph);
      }

      // add hits into each TPC graph
      for (const art::Ptr<recob::Hit> &hit: hits) {
        // check if valid
        if (!hit->WireID()) continue;

        // get the plane ID
        const geo::PlaneID this_plane_id = hit->WireID().asPlaneID();
        int index = 0;
        bool set = false;
        for (const geo::PlaneID &plane_id: _plane_ids) {
          if (this_plane_id.cmp(plane_id) == 0) { 
            set = true;
            break;
          }
          index ++;
        }
        assert(set);

        // add the point
        // this gets just the wire # part of the ID
        double wire_no = (double) hit->WireID().deepestIndex();
        // get the average time point
        double start_time = manager->GetDetectorClocksProvider()->TPCTick2Time(hit->StartTick());
        double end_time = manager->GetDetectorClocksProvider()->TPCTick2Time(hit->EndTick());
        double time = (start_time + end_time) / 2.;

        // make sure range is ok
        assert(time > 0 && time < 2000.);
        assert(wire_no >= 0 && wire_no < manager->GetGeometryProvider()->Nwires(this_plane_id));

        plane_graph_points[index][0].push_back(wire_no);
        plane_graph_points[index][1].push_back(time);
      }

      /*
      // also add points for each TPC
      unsigned graph_id = 0;
      for (auto const &cryo: _manager->GetGeometryProvider()->IterateCryostats()) {
        geo::GeometryCore::TPC_iterator iTPC = _manager->GetGeometryProvider()->begin_TPC(cryo.ID()),
                                        tend = _manager->GetGeometryProvider()->end_TPC(cryo.ID());
        while (iTPC != tend) {
          geo::TPCGeo const& TPC = *iTPC;
          for (TVector3 const &point: points) {
            if (TPC.ContainsPosition(point)) {
              // get the distance to the plane
              double height = TPC.DistanceFromReferencePlane(point);
              // and the wire in the plane
              int plane_no = 0;
              for (auto const &plane_id: _manager->GetGeometryProvider()->IteratePlaneIDs(TPC.ID())) {
                int wire_id;
		try {
                  wire_id = _manager->GetGeometryProvider()->NearestWire(point, plane_id); //.deepestIndex();
		}
                catch (geo::InvalidWireError const& e) {
                  continue;
                }
                int this_graph_id = graph_id + plane_no;
                plane_graph_points[this_graph_id][0].push_back((double) wire_id);
                plane_graph_points[this_graph_id][1].push_back(height);
                plane_no ++;
              }
            }
          }
          graph_id += TPC.Nplanes();
          iTPC++;
        }
      }*/

      std::vector<TGraph *> ret;
      unsigned ind = 0;
      for (auto const &plane_graph_xy: plane_graph_points) {
        TGraph *graph = new TGraph(plane_graph_xy[0].size(), &plane_graph_xy[0][0], &plane_graph_xy[0][1]);
        plane_graphs[ind].graph->Add(graph, "");
        plane_graphs[ind].tracks.push_back(graph);
        ret.push_back(graph);
        ind += 1;
      }

      return ret;
    }
  };

  bool containedInAV(TVector3& point);

  std::vector<TrackPlot> _event_plots;
  unsigned _event_index;
  Config _config;

};

  }  // namespace SBNOsc
}  // namespace ana

#endif  // __sbnanalysis_ana_SBNOsc_NumuSelection__

