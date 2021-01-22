////////////////////////////////////////////////////////////////////////
// Class:       PCAngleKinkFinderFinder
// Plugin Type: producer (art v3_02_06)
// File:        PCAngleKinkFinder_module.cc
//
// Generated at Wed Feb 19 17:38:21 2020 by Gray Putnam using cetskelgen
// from cetlib version v3_07_02.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "lardata/Utilities/AssociationUtil.h"

#include "lardataobj/RecoBase/PFParticle.h"

#include <memory>

#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "lardataalg/DetectorInfo/DetectorPropertiesStandard.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataobj/RecoBase/PFParticle.h"

#include "sbncode/TPCReco/Products/PCAnglePlane.h"
#include "sbncode/TPCReco/Products/PCAngleKink.h"

#include "TFitter.h"
#include "TMinuit.h"

namespace sbn {
  class PCAngleKinkFinder;
}


class sbn::PCAngleKinkFinder : public art::EDProducer {
public:
  explicit PCAngleKinkFinder(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  PCAngleKinkFinder(PCAngleKinkFinder const&) = delete;
  PCAngleKinkFinder(PCAngleKinkFinder&&) = delete;
  PCAngleKinkFinder& operator=(PCAngleKinkFinder const&) = delete;
  PCAngleKinkFinder& operator=(PCAngleKinkFinder&&) = delete;

  // Required functions.
  void produce(art::Event& e) override;

  struct FitResult {
    float chi2;
    float angle;
    float pitch;
  };

private:
  // config
  std::string fPFParticleLabel;
  std::string fAngleLabel;
  float fHitThreshold;
  float fMinDist;
  bool fAllowIncomplete;
  float fPCARadius;
  float fFitResolution;


  // private helpers
  FitResult Fit(const std::vector<sbn::PCAngle> &angles, unsigned lo, unsigned hi, unsigned center, float est_angle, float est_pitch);

};

// TODO: FIX
// TMinuit sux, so we have to define the data as a static global variable
// private data
static std::vector<float> fFitDataX;
static std::vector<float> fFitDataY;
static unsigned fFitCenter;

// static helper functions
sbn::PCAngleKink BuildKink(const sbn::PCAngle &max, const sbn::PCAngle &lo, const sbn::PCAngle &hi, const sbn::PCAngleKinkFinder::FitResult &fit, float est_angle) {
  sbn::PCAngleKink kink;
  kink.maxWire = max.wire;
  
  kink.fwhm_distance = sqrt((hi.hit_time_cm - lo.hit_time_cm) * (hi.hit_time_cm - lo.hit_time_cm) 
      + (hi.hit_wire_cm - lo.hit_wire_cm) * (hi.hit_wire_cm - lo.hit_wire_cm)); 
  kink.max_angle = max.angle;
  kink.est_angle = est_angle;

  kink.position_max = {max.hit_time_cm, max.hit_wire_cm};
  kink.position_lo = {lo.hit_time_cm, lo.hit_wire_cm};
  kink.position_hi = {hi.hit_time_cm, hi.hit_wire_cm};
  
  kink.vec_lo_at_max = max.lo_vector;
  kink.vec_hi_at_max = max.hi_vector;
  
  kink.vec_lo_at_halfmax_lo = lo.lo_vector;
  kink.vec_hi_at_halfmax_lo = lo.hi_vector;
  
  kink.vec_lo_at_halfmax_hi = hi.lo_vector;
  kink.vec_hi_at_halfmax_hi = hi.hi_vector;

  kink.fit_chi2 = fit.chi2;
  kink.fit_angle = fit.angle;
  kink.fit_pitch = fit.pitch;

  return kink;
}

// flatten the structure of each branch
std::tuple<std::vector<sbn::PCAngle>, unsigned> FlattenBranch(const sbn::PCAnglePlane &plane, unsigned i_branch, bool allow_incomplete, float min_dist) {
  std::vector<sbn::PCAngle> ret;

  std::cout << "Flattening angles for branch: " << i_branch << " plane: " << plane.plane.Plane << " hierarchy: ";
  for (int id: plane.branchHierarchy[i_branch]) std::cout << id << " ";
  std::cout << std::endl;

  // start at the end of the branch and count backwards
  for (unsigned i_branch_hierarchy = 0; i_branch_hierarchy < plane.branchHierarchy[i_branch].size(); i_branch_hierarchy++) {
    // look up the index of the current branch
    unsigned this_branch = std::distance(plane.branchIDs.begin(), std::find(plane.branchIDs.begin(), plane.branchIDs.end(), plane.branchHierarchy[i_branch][i_branch_hierarchy])); 
    if (this_branch == plane.branchIDs.size()) {
      std::cout << "BAD: UNFOUND BRANCH ID\n";
      continue;
    }
    std::cout << "BRANCH CHECK: " << plane.branchHierarchy[i_branch][i_branch_hierarchy] << " " << this_branch << " " << i_branch << std::endl;
    for (int i = plane.angles[this_branch].size()-1; i >= 0; i--) {
      // invalid if incomplete
      bool is_complete = (plane.angles[this_branch][i].complete || allow_incomplete) &&
           plane.angles[this_branch][i].dist_to_lo > min_dist && plane.angles[this_branch][i].dist_to_hi > min_dist &&
           plane.angles[this_branch][i].angle > -99.;
       
      if (!is_complete) continue;

      // Check if this hit is replicated in the following branch
      bool found_hit = false;
      if (i_branch_hierarchy > 0) {
        int i_check_branch_hierarchy = i_branch_hierarchy - 1;
        unsigned chk_branch = plane.branchIDs.size();
        int check_hitID = plane.angles[this_branch][i].hitID;
        while (i_check_branch_hierarchy >= 0 && chk_branch == plane.branchIDs.size()) {
          chk_branch = std::distance(plane.branchIDs.begin(), std::find(plane.branchIDs.begin(), plane.branchIDs.end(), plane.branchHierarchy[i_branch][i_check_branch_hierarchy]));
          i_check_branch_hierarchy --;
        }
        if (chk_branch < plane.branchIDs.size()) {
          for (unsigned j = 0; j < plane.angles[chk_branch].size(); j++) {
            if (plane.angles[chk_branch][j].hitID == check_hitID) {
              found_hit = true;
              break;
            }
          }
        }
      }
      // invalid if replicated
      if (found_hit) continue;

      // Valid!
      ret.push_back(plane.angles[this_branch][i]);
    }
  }

  // reverse the angles
  std::reverse(ret.begin(), ret.end());

  // get the number of valid hits in the main branch
  unsigned n_branch_hits = 0;
  for (unsigned i = 0; i < plane.angles[i_branch].size(); i++) {
    n_branch_hits += (plane.angles[i_branch][i].complete || allow_incomplete) &&
           plane.angles[i_branch][i].dist_to_lo > min_dist && plane.angles[i_branch][i].dist_to_hi > min_dist &&
           plane.angles[i_branch][i].angle > -99.;
  }

  // get the index of the start of this branch
  unsigned branch_start = ret.size() - n_branch_hits;
  return {ret, branch_start};
}

// Estimate peak height -- try the mean of the three largest points
float EstimateHeight(const std::vector<sbn::PCAngle> &angles, unsigned lo, unsigned hi) {
  static const unsigned NAVG = 3;
  static const unsigned IAVG = (NAVG - 1) / 2;

  // if the "peak is too small -- ignore
  if (hi - lo + 1 < NAVG) return -1.;

  std::vector<float> avgd;
  // average out the vector
  for (unsigned i = IAVG + lo; i <= hi - IAVG; i++) {
    float sum = 0.;
    for (unsigned j = i - IAVG; j < i - IAVG + NAVG; j++) {
      sum += angles[j].angle;
    }

    avgd.push_back(sum / NAVG);
  }

  return *std::max_element(avgd.begin(), avgd.end());
}


int FindMaxIndex(const std::vector<sbn::PCAngle> &angles, unsigned lo, unsigned hi) {
  static const unsigned NAVG = 3;
  static const unsigned IAVG = (NAVG - 1) / 2;

  // if the "peak is too small -- ignore
  if (hi - lo + 1 < (int)NAVG) return -1.;

  std::vector<float> avgd;
  // average out the vector
  for (unsigned i = lo + IAVG; i <= hi - IAVG; i++) {
    float sum = 0.;
    for (unsigned j = i - IAVG; j < i - IAVG + NAVG; j++) {
      sum += angles[j].angle;
    }

    avgd.push_back(sum / NAVG);
  }

  // index of the averaged
  unsigned avgd_index = std::distance(avgd.begin(), std::max_element(avgd.begin(), avgd.end()));

  // convert this to index in angles list
  unsigned angles_index = avgd_index + IAVG + lo;

  return angles_index;
} 

int FindFWHMHiIndex(const std::vector<sbn::PCAngle> &angles, unsigned start, unsigned end, float amp, float threshold) {
  // handle bad amplitude
  if (amp < 0.) return -1; 

  // go below half-the amplitude
  int halfmax_ind = -1;

  // forward track to the end of this branch
  //
  // We always identify hits on the last part of any branch, so we don't need 
  // to worry about forward tracking to any other part of the branch
  for (unsigned i = end; i < angles.size(); i++) {
    if (angles[i].angle < amp/2.) {
      halfmax_ind = i;
      break;
    } 
  }
  return halfmax_ind;
} 

int FindFWHMLoIndex(const std::vector<sbn::PCAngle> &angles, unsigned start, unsigned end, float amp, float threshold) {
  // handle bad amplitude
  if (amp < 0.) return -1; 

  // backtrack
  int halfmax_ind = -1;

  for (int i = start; i >= 0; i--) {
    if (angles[i].angle < amp / 2.) {
      halfmax_ind = i;      
      break;
    }
  }

  return halfmax_ind;
}

sbn::PCAngleKinkFinder::PCAngleKinkFinder(fhicl::ParameterSet const& p)
  : EDProducer{p},
    fPFParticleLabel(p.get<std::string>("PFParticleLabel")),
    fAngleLabel(p.get<std::string>("AngleLabel")),
    fHitThreshold(p.get<float>("HitThreshold")),
    fMinDist(p.get<float>("MinDist")),
    fAllowIncomplete(p.get<bool>("AllowIncomplete")),
    fPCARadius(p.get<float>("PCARadius")),
    fFitResolution(p.get<float>("FitResolution", 0.05))
{
  produces<std::vector<sbn::PCAngleKink>>();
  produces<art::Assns<recob::PFParticle, sbn::PCAngleKink>>();

}

double ExpectedCosAngle(double dist_to_center, double radius, double mean_distance, double angle) {
  if (dist_to_center >= radius) return 1.;

  double theta = angle;
  double A = dist_to_center;
  double B = radius - dist_to_center;
  double R = radius;
  double l = mean_distance;

  // double M01 = (B*(B + lB)*(3*A*(A + lA)*lB + (B*B*lA + B*(4*A - lA)*lB + 2*A*lB*lB)*cos(theta))*sin(theta))/(12.*lB*(B*lA + A*lB));
  //double M00 = (A*(A + lA)*lB*(2*B*lA*(2*A + lA) + A*(A - lA)*lB) + B*(B + lB)*cos(theta)*(6*A*lA*(A + lA)*lB + lA*(B*B*lA + B*(4*A - lA)*lB + 2*A*lB*lB)*cos(theta)))/(12.*lA*lB*(B*lA + A*lB));
  //double M11 = (B*(B + lB)*(B*B*lA + B*(4*A - lA)*lB + 2*A*lB*lB)*sin(theta)*sin(theta))/(12.*lB*(B*lA + A*lB));

  double M01 = (B*(B + l)*(-6*A*(A + l)*(A + B - 2*R)*sin(theta) + (3*B*(A + B)*(B + l) - 6*B*(B + l)*R + 2*(2*B + l)*R*R)*sin(2*theta)))/(24.*l*R*R);
  double M00 = (A*(A + l)*(3*A*(A + B)*(A + l) - 6*A*(A + l)*R + 2*(2*A + l)*R*R) + B*(B + l)*cos(theta)*(-6*A*(A + l)*(A + B - 2*R) + (3*B*(A + B)*(B + l) - 6*B*(B + l)*R + 2*(2*B + l)*R*R)*cos(theta)))/(12.*l*R*R);
  double M11 = (B*(B + l)*(3*B*(A + B)*(B + l) - 6*B*(B + l)*R + 2*(2*B + l)*R*R)*sin(theta)*sin(theta))/(12.*l*R*R);

  double trace = M00 + M11;
  double det = M00 * M11 - M01 * M01;
  double eigen = 0.5 * (trace + sqrt(trace * trace - 4 *det));
  double vecA = M01;
  double vecB = eigen - M00;

  return vecA / sqrt(vecA * vecA + vecB * vecB);
}

void FitLikelihood(int &npar, double *g, double &result, double *par, int flag) {
  double angle = par[0];
  double res = par[1];
  double radius = par[2];
  double mean_distance = par[3];

  double chi2 = 0.;
  for (unsigned i = 0; i < fFitDataX.size(); i++) {
    double expectation = acos(ExpectedCosAngle(fFitDataX[i], radius, mean_distance, angle));
    chi2 += (expectation - fFitDataY[i]) * (expectation - fFitDataY[i]) / (fFitDataY[i] * fFitDataY[i] * res * res);
  }

  result = chi2 / fFitDataX.size();
}

sbn::PCAngleKinkFinder::FitResult sbn::PCAngleKinkFinder::Fit(const std::vector<sbn::PCAngle> &angles, unsigned lo, unsigned hi, unsigned center, float est_angle, float est_pitch) {
  // setup the fitdata
  fFitDataX.clear();
  fFitDataY.clear();
  for (unsigned i = lo; i <= hi; i++) {
    fFitDataX.push_back( sqrt((angles[i].hit_time_cm - angles[center].hit_time_cm) * (angles[i].hit_time_cm - angles[center].hit_time_cm) +
      (angles[i].hit_wire_cm - angles[center].hit_wire_cm) * (angles[i].hit_wire_cm - angles[center].hit_wire_cm)) );
    fFitDataY.push_back(angles[i].angle);
  }
  fFitCenter = center;

  // setup the fitter
  TFitter fitter(5); // 5 param
  fitter.SetFCN(&FitLikelihood);
  fitter.SetParameter(0, "angle", est_angle, 0.1, 0, M_PI);
  fitter.SetParameter(1, "resolution", fFitResolution, 0., fFitResolution, fFitResolution);
  fitter.SetParameter(2, "radius", fPCARadius, 0., fPCARadius, fPCARadius);
  fitter.SetParameter(3, "pitch", est_pitch, 0., est_pitch, est_pitch);
  fitter.GetMinuit()->FixParameter(1);
  fitter.GetMinuit()->FixParameter(2);
  fitter.GetMinuit()->FixParameter(3);

  fitter.ExecuteCommand("MIGRAD", 0, 0);

  sbn::PCAngleKinkFinder::FitResult ret;

  ret.angle = fitter.GetParameter(0);
  ret.pitch = fitter.GetParameter(3);

  double chi2;
  double edm;
  double errdef;
  int nvpar;
  int nparx;

  fitter.GetStats(chi2, edm, errdef, nvpar, nparx);
  std::cout << "Fit chi2: " << chi2 << std::endl;
  std::cout << "Fit angle: " << ret.angle << std::endl;

  ret.chi2 = chi2;

  return ret;
}

void sbn::PCAngleKinkFinder::produce(art::Event& evt)
{
  // output stuff
  std::unique_ptr<std::vector<sbn::PCAngleKink>> outKinks(new std::vector<sbn::PCAngleKink>);
  std::unique_ptr<art::Assns<recob::PFParticle, sbn::PCAngleKink>> assn(new art::Assns<recob::PFParticle, sbn::PCAngleKink>);

  art::PtrMaker<sbn::PCAngleKink> kinkPtrMaker {evt};

  // get the PFParticle's and the associated data
  art::Handle<std::vector<recob::PFParticle>> pfparticle_handle;
  evt.getByLabel(fPFParticleLabel, pfparticle_handle);

  std::vector<art::Ptr<recob::PFParticle>> pfparticles;
  art::fill_ptr_vector(pfparticles, pfparticle_handle);

  // organize the PFPlist into a map
  std::map<unsigned, art::Ptr<recob::PFParticle>> id_to_pfp;
  for (unsigned i = 0; i < pfparticles.size(); i++) {
    id_to_pfp[pfparticles[i]->Self()] = pfparticles[i];
  }

  art::FindManyP<sbn::PCAnglePlane> pfparticleAngles(pfparticles, evt, fAngleLabel);

  // iterate over each particle
  for (unsigned i_part = 0; i_part < pfparticles.size(); i_part++) {
    art::Ptr<recob::PFParticle> thisPFP = pfparticles[i_part];
    std::cout << "NEW Particle: " << thisPFP->Self() << std::endl;
    const std::vector<art::Ptr<sbn::PCAnglePlane>> &thisAngles = pfparticleAngles.at(i_part);

    if (!thisAngles.size()) continue;

    std::vector<sbn::PCAngleKink> thisKinks;
    for (const art::Ptr<sbn::PCAnglePlane> &angle: thisAngles) {
      for (unsigned i_branch = 0; i_branch < angle->nBranches; i_branch++) {
        // flatten the structure of this branch
        auto [branchAngles, branch_start] = FlattenBranch(*angle, i_branch, fAllowIncomplete, fMinDist);
        bool in_hit = false;
        int hit_start_ind = -1;
        for (unsigned i_hit = branch_start; i_hit < branchAngles.size(); i_hit++) {
          bool this_in_hit = branchAngles[i_hit].angle > fHitThreshold;
          // New hit! set it up
          if (!in_hit && this_in_hit) {
            hit_start_ind = i_hit;
          }
          // end of hit! finish it up
          else if (in_hit && (!this_in_hit || i_hit+1 == branchAngles.size())) {
            float height = EstimateHeight(branchAngles, hit_start_ind, i_hit);

            int hi_hit_ind = FindFWHMHiIndex(branchAngles, hit_start_ind, i_hit, height, fHitThreshold);
            int lo_hit_ind = FindFWHMLoIndex(branchAngles, hit_start_ind, i_hit, height, fHitThreshold);
            // invalid hit
            // check if hit is valid
            if (hi_hit_ind >= 0 && lo_hit_ind >= 0) {
              std::cout << "NEW KINK: Particle: " << thisPFP->Self() << " Branch Ind: " << i_branch << " Plane: " << angle->plane.Plane << std::endl;
              std::cout << "EST HEIGHT: " << height << " IND LO: " << lo_hit_ind << " IND HI: " << hi_hit_ind << std::endl;

              // build out the kink information
              int max_ind = FindMaxIndex(branchAngles, lo_hit_ind, hi_hit_ind);

              // estimate the track pitch -- distance between successive hits
              float dist_lo_sum = 0.;
              for (int i = lo_hit_ind; i < max_ind; i++) {
                dist_lo_sum += sqrt((branchAngles[i].hit_time_cm - branchAngles[i+1].hit_time_cm) * (branchAngles[i].hit_time_cm - branchAngles[i+1].hit_time_cm) +
                    (branchAngles[i].hit_wire_cm - branchAngles[i+1].hit_wire_cm) * (branchAngles[i].hit_wire_cm - branchAngles[i+1].hit_wire_cm)); 
              }
              // float pitch_lo = dist_lo_sum / (max_ind - lo_hit_ind);

              float dist_hi_sum = 0.;
              for (int i = max_ind; i < hi_hit_ind; i++) {
                dist_hi_sum += sqrt((branchAngles[i].hit_time_cm - branchAngles[i+1].hit_time_cm) * (branchAngles[i].hit_time_cm - branchAngles[i+1].hit_time_cm) +
                    (branchAngles[i].hit_wire_cm - branchAngles[i+1].hit_wire_cm) * (branchAngles[i].hit_wire_cm - branchAngles[i+1].hit_wire_cm)); 
              }
              // float pitch_hi = dist_hi_sum / (hi_hit_ind - max_ind);
              float pitch = (dist_lo_sum + dist_hi_sum) / (hi_hit_ind - lo_hit_ind);


              sbn::PCAngleKinkFinder::FitResult fit = Fit(branchAngles, lo_hit_ind, hi_hit_ind, max_ind, height, pitch);

              const sbn::PCAngle &max = branchAngles[max_ind];
              const sbn::PCAngle &lo = branchAngles[lo_hit_ind];
              const sbn::PCAngle &hi = branchAngles[hi_hit_ind];

              sbn::PCAngleKink kink = BuildKink(max, lo, hi, fit, height);
              thisKinks.push_back(kink);
              // update the hit index to the end of this hit
              i_hit = hi_hit_ind;
            }
          }

          in_hit = this_in_hit;
        } 
      }
    }

    // save the kinks
    for (const sbn::PCAngleKink &k: thisKinks) {
      outKinks->push_back(k);
      art::Ptr<sbn::PCAngleKink> thisKinkPtr = kinkPtrMaker(outKinks->size()-1);
      assn->addSingle(thisPFP, thisKinkPtr);
    }
  }

  evt.put(std::move(outKinks));
  evt.put(std::move(assn));

}

DEFINE_ART_MODULE(sbn::PCAngleKinkFinder)
