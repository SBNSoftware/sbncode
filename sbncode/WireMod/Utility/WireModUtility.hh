#ifndef sbncode_WireMod_Utility_WireModUtility_hh
#define sbncode_WireMod_Utility_WireModUtility_hh

//std includes
#include <vector>

//ROOT includes
#include "TFile.h"
#include "TSpline.h"
#include "TGraph2DErrors.h"
#include "TNtuple.h"

//Framework includes
#include "larcorealg/Geometry/GeometryCore.h"
#include "lardataalg/DetectorInfo/DetectorClocksData.h"
#include "lardataalg/DetectorInfo/DetectorPropertiesData.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/TrackHitMeta.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/Simulation/SimEnergyDeposit.h"
#include "larcoreobj/SimpleTypesAndConstants/PhysicalConstants.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"

#include <limits>
#include <memory>

namespace sys {
  class WireModUtility{
    // for now make everything public, though it's probably a good idea to think about what doesn't need to be
    public:
      // detector constants, should be set by geometry service
      // the notes at the end refer to their old names in the MircoBooNE code which preceded this
      // TODO: how best to initialize the splines/graphs?
      const geo::GeometryCore* geometry;                  // save the TPC geometry
      const detinfo::DetectorClocksData& detClockData;    // save the detector clock data
      const detinfo::DetectorPropertiesData& detPropData; // save the detector property data
      bool   applyXScale;                                 // do we scale with X? (fApplyXScale)
      bool   applyYZScale;                                // do we scale with YZ? (fApplyYZScale)
      bool   applyXZAngleScale;                           // do we scale with XZ angle? (fApplyXZAngleScale)
      bool   applyYZAngleScale;                           // do we scale with YZ angle? (fApplyYZAngleScale)
      bool   applydEdXScale;                              // do we scale with dEdx? (fApplydEdXScale)
      double readoutWindowTicks;                          // how many ticks are in the readout window?
      double tickOffset;                                  // do we want an offset in the ticks? (fTickOffset)

      std::vector<TSpline3*> splines_Charge_X;            // the splines for the charge correction in X (fTSplines_Charge_X)
      std::vector<TSpline3*> splines_Sigma_X;             // the splines for the width correction in X (fTSplines_Sigma_X)
      std::vector<TSpline3*> splines_Charge_XZAngle;      // the splines for the charge correction in XZ angle (fTSplines_Charge_XZAngle)
      std::vector<TSpline3*> splines_Sigma_XZAngle;       // the splines for the width correction in XZ angle (fTSplines_Sigma_XZAngle)
      std::vector<TSpline3*> splines_Charge_YZAngle;      // the splines for the charge correction in YZ angle (fTSplines_Charge_YZAngle)
      std::vector<TSpline3*> splines_Sigma_YZAngle;       // the splines for the width correction in YZ angle (fTSplines_Sigma_YZAngle)
      std::vector<TSpline3*> splines_Charge_dEdX;         // the splines for the charge correction in dEdX (fTSplines_Charge_dEdX)
      std::vector<TSpline3*> splines_Sigma_dEdX;          // the splines for the width correction in dEdX (fTSplines_Sigma_dEdX)
      std::vector<TGraph2DErrors*> graph2Ds_Charge_YZ;    // the graphs for the charge correction in YZ (fTGraph2Ds_Charge_YZ)
      std::vector<TGraph2DErrors*> graph2Ds_Sigma_YZ;     // the graphs for the width correction in YZ (fTGraph2Ds_Sigma_YZ)

      // lets try making a constructor here
      // assume we can get a geometry service, a detector clcok, and a detector properties
      // pass the CryoStat and TPC IDs because it's IDs all the way down
      // set some optional args fpr the booleans, the readout window, and the offset
      WireModUtility(const geo::GeometryCore* geom, const detinfo::DetectorClocksData& detClock, const detinfo::DetectorPropertiesData& detProp,
                     const bool& arg_ApplyXScale = true, const bool& arg_ApplyYZScale = true, const bool& arg_ApplyXZAngleScale = true, const bool& arg_ApplyYZAngleScale = true, const bool& arg_ApplydEdXScale = true,
                     const double& arg_TickOffset = 0)
      : geometry(geom),
        detClockData(detClock),
        detPropData(detProp),
        applyXScale(arg_ApplyXScale),
        applyYZScale(arg_ApplyYZScale),
        applyXZAngleScale(arg_ApplyXZAngleScale),
        applyYZAngleScale(arg_ApplyYZAngleScale),
        applydEdXScale(arg_ApplydEdXScale),
        readoutWindowTicks(detProp.ReadOutWindowSize()),                                               // the default A2795 (ICARUS TPC readout board) readout window is 4096 samples
        tickOffset(arg_TickOffset)                                                                     // tick offset is for MC truth, default to zero and set only as necessary
      {
      }

      // typedefs
      typedef std::pair<unsigned int,unsigned int>  ROI_Key_t;
      typedef std::pair<ROI_Key_t, unsigned int> SubROI_Key_t;

      typedef struct ROIProperties
      {
        ROI_Key_t key;
        unsigned int plane;
        float begin;
        float end;
        float total_q;
        float center;   //charge weighted center of ROI
        float sigma;    //charge weighted RMS of ROI
      } ROIProperties_t;

      typedef struct SubROIProperties
      {
        SubROI_Key_t key;
        unsigned int plane;
        float total_q;
        float center;
        float sigma;
      } SubROIProperties_t;

      typedef struct ScaleValues
      {
        double r_Q;
        double r_sigma;
      } ScaleValues_t;

      typedef struct TruthProperties
      {
        float x;
        float x_rms;
        float x_rms_noWeight;
        float tick;
        float tick_rms;
        float tick_rms_noWeight;
        float total_energy;
        float x_min;
        float x_max;
        float tick_min;
        float tick_max;
        float y;
        float z;
        double dxdr;
        double dydr;
        double dzdr;
        double dedr;
        ScaleValues_t scales_avg[3];
      } TruthProperties_t;

      // internal containers
      std::map< ROI_Key_t,std::vector<size_t> > ROIMatchedEdepMap;
      std::map< ROI_Key_t,std::vector<size_t> > ROIMatchedHitMap;

      // some useful functions
      // geometries
      double planeXToTick(double xPos, int planeNumber, const geo::TPCGeo& tpcGeom, double offset = 0) { return detPropData.ConvertXToTicks(xPos, tpcGeom.Plane(planeNumber).ID()) + offset; }

      bool planeXInWindow(double xPos, int planeNumber, const geo::TPCGeo& tpcGeom, double offset = 0)
      {
        double tick = planeXToTick(xPos, planeNumber, tpcGeom, offset);
        return (tick > 0 && tick <= detPropData.ReadOutWindowSize());
      }
      
      // for this function: in the future if we want to use non-gaussian functions make this take a vector of parameters
      // the another wiremod utility could overwrite the ``fitFunc'' with some non-standard function
      // would require a fair bit of remodling (ie q and sigma would need to be replace with, eg, funcVar[0] and funcVar[1] and probs a bunch of loops)
      // so lets worry about that later
      double gausFunc(double t, double mean, double sigma, double a = 1.0)
      {
        return (a / (sigma * std::sqrt(2 * util::pi()))) * std::exp(-0.5 * std::pow((t - mean)/sigma, 2));
      }

      double FoldAngle(double theta)
      {
        return (std::abs(theta) > 0.5 * util::pi()) ? util::pi() - std::abs(theta) : std::abs(theta);
      }

      double ThetaXZ_PlaneRel(double dxdr, double dydr, double dzdr, double planeAngle)
      {
        double planeAngleRad = planeAngle * (util::pi() / 180.0);
        double sinPlaneAngle = std::sin(planeAngleRad);
        double cosPlaneAngle = std::cos(planeAngleRad);

        //double dydrPlaneRel = dydr * cosPlaneAngle - dzdr * sinPlaneAngle; // don't need to rotate Y for this angle
        double dzdrPlaneRel = dzdr * cosPlaneAngle + dydr * sinPlaneAngle;

        double theta = std::atan2(dxdr, dzdrPlaneRel);
        return FoldAngle(theta);
      }

      double ThetaYZ_PlaneRel(double dxdr, double dydr, double dzdr, double planeAngle)
      {
        double planeAngleRad = planeAngle * (util::pi() / 180.0);
        double sinPlaneAngle = std::sin(planeAngleRad);
        double cosPlaneAngle = std::cos(planeAngleRad);

        double dydrPlaneRel = dydr * cosPlaneAngle - dzdr * sinPlaneAngle;
        double dzdrPlaneRel = dzdr * cosPlaneAngle + dydr * sinPlaneAngle;

        double theta = std::atan2(dydrPlaneRel, dzdrPlaneRel);
        return FoldAngle(theta);
      }
      
      // theste are set in the .cc file
      ROIProperties_t CalcROIProperties(recob::Wire::RegionsOfInterest_t::datarange_t const&);

      std::vector<std::pair<unsigned int, unsigned int>> GetTargetROIs(sim::SimEnergyDeposit const&, double offset);
      std::vector<std::pair<unsigned int, unsigned int>> GetHitTargetROIs(recob::Hit const&);

      void FillROIMatchedEdepMap(std::vector<sim::SimEnergyDeposit> const&, std::vector<recob::Wire> const&, double offset);
      void FillROIMatchedHitMap(std::vector<recob::Hit> const&, std::vector<recob::Wire> const&);

      std::vector<SubROIProperties_t> CalcSubROIProperties(ROIProperties_t const&, std::vector<const recob::Hit*> const&);

      std::map<SubROI_Key_t, std::vector<const sim::SimEnergyDeposit*>> MatchEdepsToSubROIs(std::vector<SubROIProperties_t> const&, std::vector<const sim::SimEnergyDeposit*> const&, double offset);

      TruthProperties_t CalcPropertiesFromEdeps(std::vector<const sim::SimEnergyDeposit*> const&, double offset);

      ScaleValues_t GetScaleValues(TruthProperties_t const&,ROIProperties_t const&);
      ScaleValues_t GetScaleValues(TruthProperties_t const&,unsigned int const&);

      void ModifyROI(std::vector<float> &,
                     ROIProperties_t const &,
                     std::vector<SubROIProperties_t> const&,
                     std::map<SubROI_Key_t, ScaleValues_t> const&);
  }; // end class
} // end namespace

#endif // sbncode_WireMod_Utility_WireModUtility_hh
