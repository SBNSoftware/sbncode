/**
 * Dark Photon Decay Tool
 */

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Persistency/Common/PtrMaker.h"
#include "art/Utilities/ToolMacros.h"
#include "cetlib/cpu_timer.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"


#include "sbnobj/Common/EventGen/MeVPrtl/MeVPrtlFlux.h"
#include "sbncode/EventGenerator/MeVPrtl/Tools/IMeVPrtlDecay.h"
#include "sbncode/EventGenerator/MeVPrtl/Tools/Constants.h"

#include "CLHEP/Random/RandFlat.h"
#include "TDatabasePDG.h"

#include <fstream>
#include <cmath>
#include <vector>
#include <string>

namespace evgen {
  namespace ldm {

    class DPMakeDecay : public IMeVPrtlDecay {

    public:
      DPMakeDecay(fhicl::ParameterSet const &pset);
      ~DPMakeDecay() {}

      void configure(fhicl::ParameterSet const &pset) override;

      bool Decay(const MeVPrtlFlux &flux,
		 const TVector3 &in,
		 const TVector3 &out,
		 MeVPrtlDecay &decay,
		 double &weight) override;

      double MaxWeight() override { return 1.0; }

    private:

      double fEpsilon;
      bool fVerbose;

      std::vector<double> fMassBR;
      std::vector<double> fBRElec;
      std::vector<double> fBRHad;
      std::vector<double> fBRMu;

      double AlphaEM() const { return 1.0/137.0; }

      double GammaLL(double m, double mlep) const {
	if (m <= 2*mlep) return 0.0;

	double coupsq = fEpsilon*fEpsilon*AlphaEM();
	return (coupsq/3.0)*m*
	  (1 + 2*mlep*mlep/(m*m))*
	  std::sqrt(1 - 4*mlep*mlep/(m*m));
      }

      double InterpBR(const std::vector<double>& mass,
		      const std::vector<double>& br,
		      double m) const {

	if (mass.empty()) return 0.0;

	size_t best = 0;
	double mindiff = std::abs(mass[0] - m);

	for (size_t i=1;i<mass.size();i++) {
	  double diff = std::abs(mass[i]-m);
	  if (diff < mindiff) {
	    mindiff = diff;
	    best = i;
	  }
	}

	return br[best];
      }

      double TotalWidth(double m) const {

	double me = Constants::Instance().elec_mass;
	double mmu = Constants::Instance().muon_mass;

	if (m < 0.29) {
	  return GammaLL(m, me) + GammaLL(m, mmu);
	}

	double br_mu = InterpBR(fMassBR, fBRMu, m);
	double gamma_mu = GammaLL(m, mmu);

	if (br_mu > 0)
	  return gamma_mu / br_mu;

	return GammaLL(m, me);
      }

    };

    DPMakeDecay::DPMakeDecay(fhicl::ParameterSet const &pset)
      : IMeVPrtlStage("DPMakeDecay")
    {
      configure(pset);
    }

    void DPMakeDecay::configure(fhicl::ParameterSet const &pset)
    {
      fVerbose = pset.get<bool>("Verbose", false);
      fEpsilon = pset.get<double>("Epsilon");

      std::string file_ee = "/pnfs/sbnd/persistent/users/rohanr/dp/br/bfrac_dark_photon_e_e.txt";
      std::string file_had = "/pnfs/sbnd/persistent/users/rohanr/dp/br/bfrac_dark_photon_hadrons.txt";
      std::string file_mu  = "/pnfs/sbnd/persistent/users/rohanr/dp/br/bfrac_dark_photon_mu_mu.txt";

      std::ifstream eefile(file_ee);
      std::ifstream hadfile(file_had);
      std::ifstream mufile(file_mu);

      double m, br;

      while (eefile >> m >> br) {
	fMassBR.push_back(m);
	fBRElec.push_back(br);
      }

      while (hadfile >> m >> br) {
	fMassBR.push_back(m);
	fBRHad.push_back(br);
      }

      while (mufile >> m >> br) {
	fBRMu.push_back(br);
      }

      if (fVerbose)
	std::cout << "Loaded Dark Photon BR tables." << std::endl;
    }

    bool DPMakeDecay::Decay(const MeVPrtlFlux &flux,
			    const TVector3 &in,
			    const TVector3 &out,
			    MeVPrtlDecay &decay,
			    double &weight)
    {
      double m = flux.mass;

      double total_width = TotalWidth(m);
      if (total_width == 0) return false;

      double lifetime_ns = Constants::Instance().hbar / total_width;
      double mean_dist = lifetime_ns *
	flux.mom.Gamma() *
	flux.mom.Beta() *
	Constants::Instance().c_cm_per_ns;

      double in_dist  = (flux.pos.Vect() - in).Mag();
      double out_dist = (flux.pos.Vect() - out).Mag();

      double forced_weight =
	std::exp(-in_dist/mean_dist) -
	std::exp(-out_dist/mean_dist);

      weight = forced_weight;

      if (weight <= 0) return false;

      double me = Constants::Instance().elec_mass;
      //      double mmu = Constants::Instance().muon_mass;

      /*
      double gamma_ee = GammaLL(m, me);
      double gamma_mumu = GammaLL(m, mmu);

      double r = CLHEP::RandFlat::shoot(fEngine);

      bool decay_to_mu = false;
      if (gamma_mumu > 0 &&
	  r < gamma_mumu / total_width)
	decay_to_mu = true;
       
      double mlep = decay_to_mu ? mmu : me;
      int pdg = decay_to_mu ? 13 : 11;
      */

      double gamma_ee = GammaLL(m, me);

      double br_ee = gamma_ee / total_width;
      weight *= br_ee;

      double mlep = me;
      int pdg = 11;
      double p = std::sqrt(m*m/4. - mlep*mlep);

      TVector3 dir = RandomUnitVector();

      TLorentzVector lplus( p*dir, std::sqrt(p*p + mlep*mlep) );
      TLorentzVector lminus(-p*dir, std::sqrt(p*p + mlep*mlep) );

      lplus.Boost(flux.mom.BoostVector());
      lminus.Boost(flux.mom.BoostVector());

      double flat = CLHEP::RandFlat::shoot(fEngine,0,1);
      double decay_dist = -mean_dist*std::log(1-flat);

      TVector3 decay_pos =
	flux.pos.Vect() +
	decay_dist * (in - flux.pos.Vect()).Unit();

      decay.pos = TLorentzVector(decay_pos,
				 TimeOfFlight(flux, decay_pos));

      decay.daughter_mom.push_back(lplus.Vect());
      decay.daughter_e.push_back(lplus.E());
      decay.daughter_pdg.push_back(-pdg);

      decay.daughter_mom.push_back(lminus.Vect());
      decay.daughter_e.push_back(lminus.E());
      decay.daughter_pdg.push_back(pdg);

      decay.total_decay_width = total_width;
      decay.total_mean_lifetime = lifetime_ns;
      decay.total_mean_distance = mean_dist;
      decay.allowed_decay_fraction = 1.0;

      return true;
    }

    DEFINE_ART_CLASS_TOOL(DPMakeDecay)

  } // namespace ldm
} // namespace evgen
