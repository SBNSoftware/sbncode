#include <iostream>

#include "TDatabasePDG.h"
#include "TRandom3.h"

#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"
#include "core/Event.hh"

#include "Utilities.h"

#include "ubcore/LLBasicTool/GeoAlgo/GeoAABox.h"
#include "ubcore/LLBasicTool/GeoAlgo/GeoAlgo.h"
#include "ubcore/LLBasicTool/GeoAlgo/GeoLineSegment.h"

#include <TMath.h>

namespace ana {
  namespace SBNOsc {

void hello() {
  std::cout << "Hello SBNOsc!" << std::endl;
}


event::Interaction TruthReco(const simb::MCTruth& mctruth) {
  event::Interaction interaction;

  // Neutrino
  const simb::MCNeutrino& nu = mctruth.GetNeutrino();
  interaction.neutrino.energy = nu.Nu().EndMomentum().Energy();

  // Primary lepton
  const simb::MCParticle& lepton = nu.Lepton();
  interaction.lepton.pdg = lepton.PdgCode();
  interaction.lepton.energy = lepton.Momentum(0).Energy();
  interaction.lepton.momentum = lepton.Momentum(0).Vect();

  // Hadronic system
  for (int iparticle=0; iparticle<interaction.finalstate.size(); iparticle++) {
    event::FinalStateParticle fsp;
    const simb::MCParticle& particle = mctruth.GetParticle(iparticle);

    if (particle.Process() != "primary") {
      continue;
    }

    fsp.pdg = particle.PdgCode();
    fsp.energy = particle.Momentum(0).Energy();
    fsp.momentum = particle.Momentum(0).Vect();

    interaction.finalstate.push_back(fsp);
  }

  return interaction;
}


double ECCQE(const TVector3& l_momentum, double l_energy,
             double energy_distortion, double angle_distortion) {
  // Based on D. Kaleko, LowEnergyExcess LArLite module ECCQECalculator
  double M_n = 0.939565; // GeV/c^2
  double M_p = 0.938272; // GeV/c^2
  double M_e = 0.000511; // GeV/c^2
  double bindingE = 0.0300; // GeV

  double mp2 = M_p * M_p;
  double me2 = M_e * M_e;
  double mnb = M_n - bindingE;

  // mess around with lorentz vector
  TVector3 v(l_momentum);
  v.SetTheta( v.Theta() + angle_distortion);
  l_energy = l_energy + energy_distortion;
  double l_mom = sqrt(l_energy*l_energy - me2);
  double l_theta = \
    acos(v.Pz() / sqrt(v.Px()*v.Px() + v.Py()*v.Py() + v.Pz()*v.Pz()));
  double enu_top = mp2 - mnb*mnb - me2 + 2.0 * mnb * l_energy;
  double enu_bot = 2.0 * (mnb - l_energy + l_mom * cos(l_theta));
  return enu_top / enu_bot;
}


double NuMuOscillation(double numu_energy, double numu_dist,
                       double osc_dm2, double osc_angle) {
  double overlap = sin(2*osc_angle);
  double energy_factor = sin(1.27 * osc_dm2 * numu_dist / numu_energy);
  return 1 - overlap * overlap * energy_factor * energy_factor;
}


double containedLength(const TVector3 &v0, const TVector3 &v1,
                       const std::vector<geoalgo::AABox> &boxes) {
  static const geoalgo::GeoAlgo algo;

  // if points are the same, return 0
  if ((v0 - v1).Mag() < 1e-6) return 0;

  // construct individual points
  geoalgo::Point_t p0(v0);
  geoalgo::Point_t p1(v1);

  // construct line segment
  geoalgo::LineSegment line(p0, p1);

  double length = 0;

  // total contained length is sum of lengths in all boxes
  // assuming they are non-overlapping
  for (auto const &box: boxes) {
    int n_contained = box.Contain(p0) + box.Contain(p1);
    // both points contained -- length is total length (also can break out of loop)
    if (n_contained == 2) {
      length = (v1 - v0).Mag();
      break;
    }
    // one contained -- have to find intersection point (which must exist)
    if (n_contained == 1) {
      auto intersections = algo.Intersection(line, box);
      // Because of floating point errors, it can sometimes happen
      // that there is 1 contained point but no "Intersections"
      // if one of the points is right on the edge
      if (intersections.size() == 0) {
        // determine which point is on the edge
        double tol = 1e-5;
        bool p0_edge = algo.SqDist(p0, box) < tol;
        bool p1_edge = algo.SqDist(p1, box) < tol;
        assert(p0_edge || p1_edge);
        // contained one is on edge -- can treat both as not contained
        //
        // In this case, no length
        if ((p0_edge && box.Contain(p0)) || (box.Contain(p1) && p1_edge))
          continue;
        // un-contaned one is on edge -- treat both as contained
        else if ((p0_edge && box.Contain(p1)) || (box.Contain(p0) && p1_edge)) {
	  length = (v1 - v0).Mag();
	  break;
        }
        else {
          assert(false); // bad
        }
      }
      // floating point errors can also falsely cause 2 intersection points
      //
      // in this case, one of the intersections must be very close to the 
      // "contained" point, so the total contained length will be about
      // the same as the distance between the two intersection points
      else if (intersections.size() == 2) {
        length += (intersections.at(0).ToTLorentzVector().Vect() - intersections.at(1).ToTLorentzVector().Vect()).Mag();
        continue;
      }
      // "Correct"/ideal case -- 1 intersection point
      else if (intersections.size() == 1) {
        // get TVector at intersection point
        TVector3 int_tv(intersections.at(0).ToTLorentzVector().Vect());
        length += ( box.Contain(p0) ? (v0 - int_tv).Mag() : (v1 - int_tv).Mag() ); 
      }
      else assert(false); // bad
    }
    // none contained -- either must have zero or two intersections
    if (n_contained == 0) {
      auto intersections = algo.Intersection(line, box);
      assert(intersections.size() == 0 || intersections.size() == 2);
      if (intersections.size() == 2) {
        TVector3 start(intersections.at(0).ToTLorentzVector().Vect());
        TVector3 end(intersections.at(1).ToTLorentzVector().Vect());
        length += (start - end).Mag();
      }
    }
  }

  return length;
}


double visibleEnergy(TRandom3& rand, const simb::MCTruth &mctruth, const std::vector<sim::MCTrack> &mctrack_list, const std::vector<sim::MCShower> &mcshower_list, 
		     const VisibleEnergyCalculator &calculator, bool include_showers) {
  double visible_E = 0;

  // primary leptron track
  const sim::MCTrack *lepton_track = NULL;
  bool lepton_track_exists = false;

  // total up visible energy from tracks...
  unsigned ind = 0;
  for (auto const &mct: mctrack_list) {
    // ignore particles not from nu vertex, non primary particles, and uncharged particles
    if (!isFromNuVertex(mctruth, mct) || abs(PDGCharge(mct.PdgCode())) < 1e-4 || mct.Process() != "primary")
       continue;
    // account for primary lepton later
    if ((abs(mct.PdgCode()) == 13 || abs(mct.PdgCode()) == 11) && calculator.lepton_index == ind) {
      continue;
    }

    double mass = PDGMass(mct.PdgCode());
    double this_visible_energy = (mct.Start().E() - mass) / 1000. /* MeV to GeV */;
    if (calculator.track_energy_distortion > 1e-4) {
      this_visible_energy = rand.Gaus(this_visible_energy, calculator.track_energy_distortion*this_visible_energy);
      // clamp to 0
      this_visible_energy = std::max(this_visible_energy, 0.);
    }
    if (this_visible_energy > calculator.track_threshold) {
      visible_E += this_visible_energy;
    }
    ind ++;
  }

  // ...and showers
  if (include_showers) {
    for (auto const &mcs: mcshower_list) {
      // ignore particles not from nu vertex, non primary particles, and uncharged particles
      if (!isFromNuVertex(mctruth, mcs) || abs(PDGCharge(mcs.PdgCode())) < 1e-4 || mcs.Process() != "primary")
        continue; 
      // account for primary lepton later
      if ((abs(mcs.PdgCode()) == 13 || abs(mcs.PdgCode()) == 11) && isFromNuVertex(mctruth, mcs))
        continue; 

      double mass = PDGMass(mcs.PdgCode());
      double this_visible_energy = (mcs.Start().E() - mass) / 1000. /* MeV to GeV */;
      if (calculator.shower_energy_distortion > 1e-4) {
        this_visible_energy = rand.Gaus(this_visible_energy, calculator.shower_energy_distortion*this_visible_energy);
        // clamp to 0
        this_visible_energy = std::max(this_visible_energy, 0.);
      }
      if (this_visible_energy > calculator.shower_threshold) {
        visible_E += this_visible_energy;
      }
    }
  }

  // ...and primary lepton energy (for CC events)
  // only add in extra here if identified "lepton" is actually a lepton
  if (calculator.lepton_index >= 0 && (abs(mctrack_list[calculator.lepton_index].PdgCode()) == 13 || abs(mctrack_list[calculator.lepton_index].PdgCode()) == 11)) {
    visible_E += smearLeptonEnergy(rand,mctrack_list[calculator.lepton_index], calculator);
  }

  return visible_E;
}

//Function returns the Hadronic visibleEnergy[0], shower visibleEnergy[1] and the neutrino lepton energy visibleEnerg[2]  energy at the vertex 
    std::vector<double> FlavourEnergyDeposition(TRandom3& rand, const simb::MCTruth &mctruth, std::map<int,const simb::MCParticle*>& mcparticles,std::map<int,double>& mcvisibleparticles, std::vector<geoalgo::AABox>& Volumes, const VisibleEnergyCalculator &calculator){

      std::vector<int> used_particles;

      double Hadronic_E = 0;
      double Other_E   = 0;
      double Leptonic_E = 0;
      //      for(auto const& mcvisibleparticle: mcvisibleparticles){
      
      for(auto const& mcparticleit:  mcparticles){

	const simb::MCParticle*  mcparticle = mcparticleit.second;
	//double this_visible_energy = mcvisibleparticle.second;
	//const simb::MCParticle*  mcparticle =  mcparticles[mcvisibleparticle.first];

	double mass = PDGMass(mcparticle->PdgCode());
	
	double KEStart = (mcparticle->E()*1000. - mass)/1000.;
	double KEEnd   = (mcparticle->E(mcparticle->NumberTrajectoryPoints()-1)*1000 - mass)/1000;
	double KE = KEStart - KEEnd; 
	double this_visible_energy = KEStart - KEEnd;

	//	double this_visible_energy = (mcparticle->E()*1000. - mass) / 1000.;

	//       	std::cout << "mc particle with track id: " << mcparticle->TrackId() << " pdg: " << mcparticle->PdgCode() << " this_visible_energy: " << this_visible_energy << " start KE: " << KEStart << std::endl;

	//Ignore nuclear fragments
	if(abs(mcparticle->PdgCode()) > 5000){
	  continue;
	}

	//has to be charged 
	if(abs(PDGCharge(mcparticle->PdgCode())) < 1e-4){
	  continue;
	}

	//Has to be in the visible list
	if(mcvisibleparticles.find(mcparticle->TrackId()) == mcvisibleparticles.end()){continue;}
	
	//Anything less than an 8.3keV will be not be seen for the ADC resolution 
        if(mcvisibleparticles[mcparticle->TrackId()] < 8.3e-6){
	  continue;
        }
	

	//Consider only particles in the active volume.
	geoalgo::Point_t p(mcparticle->Position().Vect());
	bool contained = false;
	for (auto const& AV: Volumes){
	  if(AV.Contain(p)){contained = true;}
	}
	if(!contained){continue;}


	//       	Std::cout << "passed first three" << std::endl;

	//Check if it comes from the neutrino vertex
	//Is it a daugther of the the lepton candidate?
	int track_id =  mcparticle->TrackId();
        //Find the initial  mother
	int mother_id = track_id;
        while(mother_id != 0){
          if(mcparticles.find(mother_id) != mcparticles.end()){
            	   
	    //Don't account for shower will be done later - Nue Only 
	    if((track_id != mother_id) && 
	       (abs(mcparticles[track_id]->PdgCode()) == 11 || abs(mcparticles[track_id]->PdgCode()) == 22) && 
	       (abs(mcparticles[mother_id]->PdgCode()) == 11 || abs(mcparticles[mother_id]->PdgCode()) == 22)){break;}

	    track_id = mother_id;
            if(mother_id == calculator.lepton_index){break;}

	    //Try to roll up energy from delta rays into the hadrons
	    //	    if(std::find(used_particles.begin(), used_particles.end(),track_id) != used_particles.end()){break;} 
            
	    mother_id = mcparticles[mother_id]->Mother();
          }
          else{
            break;
          }
	}

	//leptons energy will be accounted for 
	if(track_id == calculator.lepton_index){continue;}

	//Ignore the showering particles 
	if((abs(mcparticle->PdgCode()) == 22 || abs(mcparticle->PdgCode()) == 11) && 
	   (abs(mcparticles[track_id]->PdgCode()) == 11 || abs(mcparticles[track_id]->PdgCode()) == 22)){continue;}
	
	//Check the particles energy has not aready been accounted for 
	//       	if(std::find(used_particles.begin(), used_particles.end(),track_id) != used_particles.end()){
	//  this_visible_energy = rand.Gaus(this_visible_energy, calculator.track_energy_distortion*this_visible_energy);
	//  this_visible_energy = std::max(this_visible_energy, 0.);
	//  Hadronic_E += this_visible_energy; 
	//  continue;
	//}


	if(!isFromNuVertex(mctruth,mcparticles[track_id])){continue;}
	
	//Hadronic Energy
	if(abs(mcparticle->PdgCode()) == 2212){
	  float track_threshold = calculator.track_threshold;
	  if(KE < track_threshold){continue;}
	  this_visible_energy = rand.Gaus(this_visible_energy, calculator.track_energy_distortion*this_visible_energy);
	  this_visible_energy = std::max(this_visible_energy, 0.);
	  Hadronic_E += this_visible_energy; 
	  used_particles.push_back(mcparticle->TrackId());
	  continue;
	}
      
	if(abs(mcparticle->PdgCode()) == 211){
	  //	  std::cout << "calculator.track_energy_distortion : " << calculator.track_energy_distortion << " this_visible_energy: " << this_visible_energy << std::endl;
	  double energy = mcparticle->E();
	  this_visible_energy = rand.Gaus(energy, calculator.track_energy_distortion*energy);
	  this_visible_energy = std::max(this_visible_energy, 0.);
	  Hadronic_E += this_visible_energy; 
	  //	  std::cout << "  hadron with track ID: " << mcparticle->TrackId() << " E: " << this_visible_energy << std::endl;
	  used_particles.push_back(mcparticle->TrackId());
	  continue;
	}

	if(abs(mcparticle->PdgCode()) == 321){
	  double energy = mcparticle->E();
	  //	  std::cout << "calculator.track_energy_distortion : " << calculator.track_energy_distortion << " this_visible_energy: " << this_visible_energy << std::endl;
	  this_visible_energy = rand.Gaus(energy, calculator.track_energy_distortion*energy);
	  this_visible_energy = std::max(this_visible_energy, 0.);
	  Hadronic_E += this_visible_energy; 
	  //	  std::cout << "  hadron with track ID: " << mcparticle->TrackId() << " E: " << this_visible_energy << std::endl;
	  used_particles.push_back(mcparticle->TrackId());
	  continue;
	}

      		
	//if did come from the lepton then add it energy
	if(track_id == calculator.lepton_index){
	  //Leptonic_E += this_visible_energy;
	  //	  used_particles.push_back(mcparticle->TrackId());
	  continue; 
	}
      
	//	this_visible_energy = (mcparticle->E());
	//	std::cout << "calculator.track_energy_distortion : " << calculator.track_energy_distortion << " this_visible_energy: " << this_visible_energy << std::endl;
	this_visible_energy = rand.Gaus(this_visible_energy, calculator.track_energy_distortion*this_visible_energy);
	this_visible_energy = std::max(this_visible_energy, 0.);
	Other_E += this_visible_energy;
	//	used_particles.push_back(mcparticle->TrackId());
	//	std::cout << "   other with track ID: " << mcparticle->TrackId() << " E1: " << mcparticle->E() << " E: " << this_visible_energy <<  " pdg: " << mcparticle->PdgCode() << std::endl;
	
      }
      
      //Calculate the lepton smearing 
      if (calculator.lepton_index >= 0 && (abs(mcparticles[calculator.lepton_index]->PdgCode()) == 13 || abs(mcparticles[calculator.lepton_index]->PdgCode()) == 11)) {
	if(mcparticles.find(calculator.lepton_index) != mcparticles.end()){
	  Leptonic_E = smearLeptonEnergy(rand,mcparticles[calculator.lepton_index], calculator);         
	}
      }
      std::vector<double>  visibleEnergy = {Hadronic_E, Other_E, Leptonic_E};
      return  visibleEnergy;
    }
    
  // std::vector<double>  visibleEnergy;
  // double visible_E = 0;
    
  // // total up visible energy from tracks...
  // for(auto const &mct: mcparticles){
    
  //   unsigned ind = mct.first;

  //   const simb::MCParticle*  mcparticle = mct.second;


  //   // ignore particles not from nu vertex, non primary particles, and uncharged particles
  //   if(abs(PDGCharge(mcparticle->PdgCode())) < 1e-4){
  //     continue;
  //   }

  //   //Ignore nuclear fragments
  //   if(abs(mcparticle->PdgCode()) > 5000){
  //     continue;
  //   }

  //   // account for primary lepton later                                               
  //   if (abs(mcparticle->PdgCode()) == 11) {
  //     continue;
  //   }

  //   if (abs(mcparticle->PdgCode()) == 13 && mcparticle->Mother() == 0){
  //     continue;
  //   }

  //   //Check particle comes from the vertex and is the first charge particle.
  //   if (!isChargedPrimary(mctruth, mcparticles, (int) mcparticle->TrackId())){
  //     continue;
  //   }

  //   //make sure the particle deposited energy 
  //   if(std::find(mcvisibleparticles.begin(),mcvisibleparticles.end(),(int) mcparticle->TrackId()) == mcvisibleparticles.end()){
  //     continue;
  //   }
    
  //   float track_threshold = calculator.track_threshold;

  //   if(abs(mcparticle->PdgCode()) != 2212){
  //     track_threshold = 0;
  //   }

  //   //std::cout << "Track passed charge + position  with id: " << mcparticle->TrackId() << " has pdgcode: " << mcparticle->PdgCode() <<  " charge: " << abs(PDGCharge(mcparticle->PdgCode())) << std::endl;
    
  //   double mass = PDGMass((mct.second)->PdgCode());
  //   double this_visible_energy = (mcparticle->E()*1000. - mass) / 1000. /* MeV to GeV */;

  //   if (this_visible_energy >= track_threshold){
  //     if(abs(mcparticle->PdgCode()) == 211 || abs(mcparticle->PdgCode()) == 321){ 
  // 	this_visible_energy = rand.Gaus(this_visible_energy, calculator.track_energy_distortion*mcparticle->E());
  //     }
  //     else if(abs(mcparticle->PdgCode()) == 2212){
	
  // 	this_visible_energy = rand.Gaus(this_visible_energy, calculator.track_energy_distortion*this_visible_energy);
  //     }
  //     else{
  // 	this_visible_energy = rand.Gaus(this_visible_energy, calculator.track_energy_distortion*mcparticle->E());
  // 	//	std::cout << "a particle with pdgcode: " <<  mcparticle->PdgCode()  << " and track id: " << mcparticle->TrackId() << " this_visible_energy: " << this_visible_energy  << std::endl;
  //     }
      
  //     this_visible_energy = std::max(this_visible_energy, 0.);
  //     visible_E += this_visible_energy;
      
  //     std::cout << "first this_visible_energy: " << this_visible_energy << " mass: " << mass << " E: " << mcparticle->E() << " threshold: " << track_threshold  << " track ID: " << mcparticle->TrackId() << " pdgcode: " << mcparticle->PdgCode() << std::endl;
  //   //if (this_visible_energy >= track_threshold) { //Move to before the distorton ala proposal
  //     //if (calculator.track_energy_distortion > 1e-4) {
  // 	//this_visible_energy = rand.Gaus(this_visible_energy, calculator.track_energy_distortion*this_visible_energy);
  // 	// clamp to 0
  // 	//this_visible_energy = std::max(this_visible_energy, 0.);
  // 	//	std::cout << " dis this_visible_energy: " << this_visible_energy << " mass: " << mass << " (mct.second)a.Start().E(): " << mcparticle->E() << std::endl;
  // 	//}
  //     // std::cout << "used" << std::endl;
  //     // visible_E += this_visible_energy;
  //   }
  // }
  
  // visibleEnergy.push_back(visible_E);
  // visible_E = 0;

  // //  ...and showers. This doesn't really matter if there is another shower the event will be cut. 
  // for(auto const &mcs: mcparticles){

  //   unsigned ind = mcs.first;

  //   const simb::MCParticle*  mcparticle = mcs.second;

  //   // ignore non shower particles and daughter showers.
  //   if(abs(mcparticle->PdgCode()) != 11 || abs(mcparticle->PdgCode()) != 22 || !isFromNuVertex(mctruth,  mcparticle))
  //     continue;

  //   // account for primary lepton later
  //   if ((abs(mcparticle->PdgCode()) == 13 || abs(mcparticle->PdgCode()) == 11) && isFromNuVertex(mctruth,  mcparticle)  && calculator.lepton_index == ind)
  //     continue; 

  //   //make sure the particle deposited energy 
  //   if(std::find(mcvisibleparticles.begin(),mcvisibleparticles.end(),(int) mcparticle->TrackId()) == mcvisibleparticles.end()){
  //     continue;
  //   }


  //   double mass = PDGMass(mcparticle->PdgCode());
  //   double this_visible_energy = (mcparticle->E()*1000. - mass) / 1000.  /* MeV to GeV */;
  //   if (calculator.shower_energy_distortion > 1e-4) {
  //     this_visible_energy = rand.Gaus(this_visible_energy, calculator.shower_energy_distortion*mcparticle->E());
  //     // clamp to 0
  //     this_visible_energy = std::max(this_visible_energy, 0.);
  //   }
  //   if (this_visible_energy > calculator.shower_threshold) {
  //     //      visible_E += this_visible_energy;
  //   }
  // }

  // visibleEnergy.push_back(visible_E);
  // visible_E = 0;

  // //  ...and primary lepton energy (for CC events)
  // //  only add in extra here if identified "lepton" is actually a lepton
  // if (calculator.lepton_index >= 0 && (abs(mcparticles[calculator.lepton_index]->PdgCode()) == 13 || abs(mcparticles[calculator.lepton_index]->PdgCode()) == 11)) {
  //      if(mcparticles.find(calculator.lepton_index) != mcparticles.end()){
  //        visible_E += smearLeptonEnergy(rand,mcparticles[calculator.lepton_index], calculator);
  //    }
  // }


  // visibleEnergy.push_back(visible_E);
  
  // return visibleEnergy;

    double smearLeptonEnergy(TRandom3& rand, const sim::MCTrack &mct, const VisibleEnergyCalculator &calculator) {
      // setup distortion
      
      double smearing_percentage;
      if (calculator.lepton_contained) {
	smearing_percentage = calculator.lepton_energy_distortion_contained;
      } else {
	double A = calculator.lepton_energy_distortion_leaving_A;
	double B = calculator.lepton_energy_distortion_leaving_B;
	smearing_percentage = -A * TMath::Log(B * calculator.lepton_contained_length);
      }
      // smear visible energy
      double lepton_visible_energy = (mct.Start().E()) / 1000.; /* MeV to GeV */
      double smeared_lepton_visible_energy = rand.Gaus(lepton_visible_energy, smearing_percentage * lepton_visible_energy);
      // clamp to 0
      smeared_lepton_visible_energy = std::max(smeared_lepton_visible_energy, 0.);
      return smeared_lepton_visible_energy;
    }
    
    
double smearLeptonEnergy(TRandom3& rand, const simb::MCParticle* &lepton, const VisibleEnergyCalculator &calculator){
  
  //Like proposal but not sure why.
  //  rand.SetSeed(0);

  
  if(lepton->PdgCode() == 13){
  double smearing_percentage;
  if (calculator.lepton_contained) {
    smearing_percentage = calculator.lepton_energy_distortion_contained;
  } else {
    double A = calculator.lepton_energy_distortion_leaving_A;
    double B = calculator.lepton_energy_distortion_leaving_B;
    smearing_percentage = -A * TMath::Log(B * calculator.lepton_contained_length);
  }
  // smear visible energy
  double lepton_visible_energy = (lepton->E());
  double smeared_lepton_visible_energy = rand.Gaus(lepton_visible_energy, smearing_percentage * lepton_visible_energy);

  // clamp to 0
  smeared_lepton_visible_energy = std::max(smeared_lepton_visible_energy, 0.);
  return smeared_lepton_visible_energy;
  }

  if(TMath::Abs(lepton->PdgCode()) == 11 || lepton->PdgCode() == 22){

    //Smear the electron energy 
    double smearing_percentage = calculator.lepton_energy_distortion_contained;
    double lepton_visible_energy = (lepton->E());
    smearing_percentage /= TMath::Sqrt(lepton_visible_energy);
    double smeared_lepton_visible_energy = rand.Gaus(lepton_visible_energy, smearing_percentage * lepton_visible_energy);
    //    std::cout << "smeared_lepton_visible_energy new: " << smeared_lepton_visible_energy <<  "lepton_visible_energy: " << lepton_visible_energy    << " smearing_percentage: " << smearing_percentage  << "sigma: " << (smeared_lepton_visible_energy-lepton_visible_energy)/(smearing_percentage * lepton_visible_energy)<< std::endl; 
    smeared_lepton_visible_energy = std::max(smeared_lepton_visible_energy, 0.);
    return smeared_lepton_visible_energy;
  }
}

double smearLeptonEnergy(TRandom3& rand, double& LeptonE, int& PdgCode,  const VisibleEnergyCalculator &calculator){

  //Like proposal but not sure why.
  rand.SetSeed(0);
  
  if(PdgCode == 13){
    double smearing_percentage;
    if (calculator.lepton_contained) {
      smearing_percentage = calculator.lepton_energy_distortion_contained;
    } 
    else {
      double A = calculator.lepton_energy_distortion_leaving_A;
      double B = calculator.lepton_energy_distortion_leaving_B;
      smearing_percentage = -A * TMath::Log(B * calculator.lepton_contained_length);
    }
    // smear visible energy
    double lepton_visible_energy = LeptonE;
    double smeared_lepton_visible_energy = rand.Gaus(lepton_visible_energy, smearing_percentage * lepton_visible_energy);
    
    // clamp to 0
    smeared_lepton_visible_energy = std::max(smeared_lepton_visible_energy, 0.);
    return smeared_lepton_visible_energy;
  }

  if(TMath::Abs(PdgCode) == 11 || PdgCode == 22){
    
    //Smear the electron energy 
    double smearing_percentage = calculator.lepton_energy_distortion_contained;
    double lepton_visible_energy = LeptonE;
    smearing_percentage /= TMath::Sqrt(lepton_visible_energy);
    double smeared_lepton_visible_energy = rand.Gaus(lepton_visible_energy, smearing_percentage * lepton_visible_energy);
    smeared_lepton_visible_energy = std::max(smeared_lepton_visible_energy, 0.);
    return smeared_lepton_visible_energy;
  }
}
  


// define global static const PDGTable to be used by helper functions
static const TDatabasePDG *PDGTable(new TDatabasePDG);


double PDGMass(int pdg) {
  // regular particle
  if (pdg < 1000000000) {
    TParticlePDG* ple = PDGTable->GetParticle(pdg);
    return ple->Mass() * 1000.0;
  }
  // ion
  else {
    int p = (pdg % 10000000) / 10000;
    int n = (pdg % 10000) / 10 - p;
    return (PDGTable->GetParticle(2212)->Mass() * p +
            PDGTable->GetParticle(2112)->Mass() * n) * 1000.0;
  }
}


double PDGCharge(int pdg) {
  // regular particle
  if (pdg < 1000000000) {
    TParticlePDG* ple = PDGTable->GetParticle(pdg);
    return ple->Charge();
  }
  // ion
  else {
    int p = (pdg % 10000000) / 10000;
    return p * 3;
  }
}


bool isFromNuVertex(const simb::MCTruth& mc, const sim::MCShower& show,
                    float distance)  {
  TLorentzVector nuVtx = mc.GetNeutrino().Nu().Trajectory().Position(0);
  TLorentzVector showStart = show.Start().Position();
  return TMath::Abs((showStart - nuVtx).Mag()) < distance;
}


bool isFromNuVertex(const simb::MCTruth& mc, const sim::MCTrack& track,
                    float distance) {
  TLorentzVector nuVtx = mc.GetNeutrino().Nu().Trajectory().Position(0);
  TLorentzVector trkStart = track.Start().Position();
  //  std::cout << "TMath::Abs((trkStart - nuVtx).Mag()): " << TMath::Abs((trkStart - nuVtx).Mag()) << " distance: " << distance << std::endl;
  return TMath::Abs((trkStart - nuVtx).Mag()) < distance;
}

    bool isFromNuVertex(const simb::MCTruth& mc, const simb::MCParticle* &particle, float distance){
  
      TLorentzVector nuVtx     = mc.GetNeutrino().Nu().Trajectory().Position(0);
      TLorentzVector partstart = particle->Position();
      return TMath::Abs((partstart - nuVtx).Mag()) < distance;
    }
    

    bool isFromNuVertex(std::vector<simb::MCTruth>& mcs,std::map<int, const simb::MCParticle*>& mcparticles, int particle_id, float distance){
      
      //Find the mother
      int track_id =  particle_id;
      if(mcparticles.find(track_id) == mcparticles.end()){ return false;}
      //Find the initial  mother 
      int mother_id = track_id;
      while(mother_id != 0){
      	if(mcparticles.find(mother_id) != mcparticles.end()){
      	  track_id = mother_id;
      	  mother_id = mcparticles[mother_id]->Mother();
      	}
      	else{
      	  mother_id = track_id;
      	  break;
      	}
      }
      
      //Check it is the same origin 
      // for(auto const mc: mcs){
      //      	bool isNuVertex = isFromNuVertex(mc,mcparticles[track_id]); 
      //      	if(isNuVertex){return true;}
      //      }
      return false;
    }
  


    // bool isChargedPrimary(const simb::MCTruth& mc, std::map<int, const simb::MCParticle*>& mcparticles, int particle_id){
      
    //   int track_id =  particle_id;
    //   if(mcparticles.find(track_id) == mcparticles.end()){ return false;}
    //   //Find the initial  mother 
    //   int mother_id = track_id;
    //   while(mother_id != 0){
    // 	if(mcparticles.find(mother_id) != mcparticles.end()){
    // 	  //If the mother is charge the energy is accounted for.
    // 	  if(abs(PDGCharge(mcparticles[mother_id]->PdgCode())) > 1e-4 && track_id != mother_id){return false;}
    // 	  track_id = mother_id;
    // 	  mother_id = mcparticles[mother_id]->Mother();
    // 	}
    // 	else{
    // 	  return false;
    // 	}
    //   }
    //   //Check it is the same origin 
    //   if(mother_id == 0){
    // 	bool isNuVertex = isFromNuVertex(mc,mcparticles[track_id]); 
    // 	return isNuVertex;
    //   }
    //   return false;
	
    // }

    bool isChargedPrimary(const simb::MCTruth& mc, std::map<int, const simb::MCParticle*>& mcparticles, int particle_id){
      
      int track_id =  particle_id;
      if(mcparticles.find(track_id) == mcparticles.end()){ return false;}
      //Find the initial  mother 
      int mother_id = track_id;
      while(mother_id != 0){
    	if(mcparticles.find(mother_id) != mcparticles.end()){
    	  //If the mother is charge the energy is accounted for.
    	  if(abs(PDGCharge(mcparticles[mother_id]->PdgCode())) > 1e-4 && track_id != mother_id && mcparticles[mother_id]->PdgCode() < 5000){return false;}
    	  track_id = mother_id;
    	  mother_id = mcparticles[mother_id]->Mother();
    	}
    	else{
    	  mother_id = track_id;
    	  break;
    	}
      }
      //Check it is the same origin 
    	bool isNuVertex = isFromNuVertex(mc,mcparticles[track_id]); 
    	return isNuVertex;
    }



  }
}  // namespace ana
