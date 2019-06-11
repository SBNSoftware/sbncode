#include "../include/CC1piAnalysisHelper.h"

namespace selection{
  
  //-----------------------------------------------------------------------------------------------
  
  ParticleMatrix CC1piAnalysisHelper::LengthBasedStatistics(const Event &e, const TopologyMap signal_map_topology, ParticleMatrix &count_longest, ParticleMatrix &count_second_longest){
/*
    if(e.CheckMCTopology(signal_map_topology) == 1){
      //LONGEST
      if(GeneralAnalysisHelper::GetMCLengthWithPdg( 13 ) > GeneralAnalysisHelper::GetMCLengthWithPdg( 2212 ) && GeneralAnalysisHelper::GetMCLengthWithPdg( 13 ) > GeneralAnalysisHelper::GetMCLengthWithPdg( 211 ) ) {count_longest[0][0]++;}
      else if(GeneralAnalysisHelper::GetMCLengthWithPdg( 211 ) > GeneralAnalysisHelper::GetMCLengthWithPdg( 2212 ) && GeneralAnalysisHelper::GetMCLengthWithPdg( 211 ) > GeneralAnalysisHelper::GetMCLengthWithPdg(13) ) {count_longest[0][1]++;}
      else if(GeneralAnalysisHelper::GetMCLengthWithPdg( 2212 ) > GeneralAnalysisHelper::GetMCLengthWithPdg( 211 ) && GeneralAnalysisHelper::GetMCLengthWithPdg( 2212 ) > GeneralAnalysisHelper::GetMCLengthWithPdg(13) ) {count_longest[0][2]++;}

      if( e.CheckRecoTopology(GeneralAnalysisHelper::GetCC1PiTopologyMap()) == 1 ){
        if( GeneralAnalysisHelper::GetRecoLengthWithPdg( 13 ) > GeneralAnalysisHelper::GetRecoLengthWithPdg( 2212 ) && GeneralAnalysisHelper::GetRecoLengthWithPdg( 13 ) > GeneralAnalysisHelper::GetRecoLengthWithPdg( 211 ) ){count_longest[1][0]++ ;}
        else if( GeneralAnalysisHelper::GetRecoLengthWithPdg( 211 ) > GeneralAnalysisHelper::GetRecoLengthWithPdg( 2212 ) && GeneralAnalysisHelper::GetRecoLengthWithPdg( 211 ) > GeneralAnalysisHelper::GetRecoLengthWithPdg( 13 ) ) {count_longest[1][1]++;} 
        else if( GeneralAnalysisHelper::GetRecoLengthWithPdg( 2212 ) > GeneralAnalysisHelper::GetRecoLengthWithPdg( 211 ) && GeneralAnalysisHelper::GetRecoLengthWithPdg( 2212 ) > GeneralAnalysisHelper::GetRecoLengthWithPdg( 13 ) ) {count_longest[1][2]++;}
      }

      if( e.CheckRecoTopology( GeneralAnalysisHelper::GetCC1PiTopologyMap() ) == 1 ) {
        if( GeneralAnalysisHelper::GetRecoLengthWithPdg( 13 ) > GeneralAnalysisHelper::GetRecoLengthWithPdg( 2212 ) && GeneralAnalysisHelper::GetRecoLengthWithPdg( 13 ) > GeneralAnalysisHelper::GetRecoLengthWithPdg( 211 ) ){count_longest[2][0]++ ;}
        else if( GeneralAnalysisHelper::GetRecoLengthWithPdg( 211 ) > GeneralAnalysisHelper::GetRecoLengthWithPdg( 2212 ) && GeneralAnalysisHelper::GetRecoLengthWithPdg( 211 ) > GeneralAnalysisHelper::GetRecoLengthWithPdg( 13 ) ) {count_longest[2][1]++;} 
        else if( GeneralAnalysisHelper::GetRecoLengthWithPdg( 2212 ) > GeneralAnalysisHelper::GetRecoLengthWithPdg( 211 ) && GeneralAnalysisHelper::GetRecoLengthWithPdg( 2212 ) > GeneralAnalysisHelper::GetRecoLengthWithPdg( 13 ) ) {count_longest[2][2]++;}
      }

      //SECOND LONGEST
      if( GeneralAnalysisHelper::GetMCLengthWithPdg( 13 ) > GeneralAnalysisHelper::GetMCLengthWithPdg( 211 ) && GeneralAnalysisHelper::GetMCLengthWithPdg( 211 ) > GeneralAnalysisHelper::GetMCLengthWithPdg( 2212 ) ) {count_second_longest[0][1]++;} 
      else if( GeneralAnalysisHelper::GetMCLengthWithPdg( 2212 ) > GeneralAnalysisHelper::GetMCLengthWithPdg( 211 ) && GeneralAnalysisHelper::GetMCLengthWithPdg( 211 ) > GeneralAnalysisHelper::GetMCLengthWithPdg( 13 ) ) {count_second_longest[0][1]++;} 
      else if( GeneralAnalysisHelper::GetMCLengthWithPdg( 211 ) > GeneralAnalysisHelper::GetMCLengthWithPdg( 13 ) && GeneralAnalysisHelper::GetMCLengthWithPdg( 13 ) > GeneralAnalysisHelper::GetMCLengthWithPdg( 2212 ) ) {count_second_longest[0][0]++;} 
      else if( GeneralAnalysisHelper::GetMCLengthWithPdg( 2212 ) > GeneralAnalysisHelper::GetMCLengthWithPdg( 13 ) && GeneralAnalysisHelper::GetMCLengthWithPdg( 13 ) > GeneralAnalysisHelper::GetMCLengthWithPdg( 211 ) ) {count_second_longest[0][0]++;} 
      else if( GeneralAnalysisHelper::GetMCLengthWithPdg( 13 ) > GeneralAnalysisHelper::GetMCLengthWithPdg( 2212 ) && GeneralAnalysisHelper::GetMCLengthWithPdg( 2212 ) > GeneralAnalysisHelper::GetMCLengthWithPdg( 211 ) ) {count_second_longest[0][2]++; std::cout<<"Hello"<<std::endl;} 
      else if( GeneralAnalysisHelper::GetMCLengthWithPdg( 211 ) > GeneralAnalysisHelper::GetMCLengthWithPdg( 2212 ) && GeneralAnalysisHelper::GetMCLengthWithPdg( 2212 ) > GeneralAnalysisHelper::GetMCLengthWithPdg( 13 ) ) {count_second_longest[0][2]++; std::cout<<"Hello"<<std::endl;}
	
      if( e.CheckRecoTopology(GeneralAnalysisHelper::GetCC1PiTopologyMap()) == 1 ){
        if( GeneralAnalysisHelper::GetRecoLengthWithPdg( 13 ) > GeneralAnalysisHelper::GetRecoLengthWithPdg( 211 ) && GeneralAnalysisHelper::GetRecoLengthWithPdg( 211 ) > GeneralAnalysisHelper::GetRecoLengthWithPdg( 2212 ) ) {count_second_longest[1][1]++;}
        else if( GeneralAnalysisHelper::GetRecoLengthWithPdg( 2212 ) > GeneralAnalysisHelper::GetRecoLengthWithPdg( 211 ) && GeneralAnalysisHelper::GetRecoLengthWithPdg( 211 ) > GeneralAnalysisHelper::GetRecoLengthWithPdg( 13 ) ) {count_second_longest[1][1]++;}
        else if( GeneralAnalysisHelper::GetRecoLengthWithPdg( 211 ) > GeneralAnalysisHelper::GetRecoLengthWithPdg( 13 ) && GeneralAnalysisHelper::GetRecoLengthWithPdg( 13 ) > GeneralAnalysisHelper::GetRecoLengthWithPdg( 2212 ) ) {count_second_longest[1][0]++;}
        else if( GeneralAnalysisHelper::GetRecoLengthWithPdg( 2212 ) > GeneralAnalysisHelper::GetRecoLengthWithPdg( 13 ) && GeneralAnalysisHelper::GetRecoLengthWithPdg( 13 ) > GeneralAnalysisHelper::GetRecoLengthWithPdg( 211 ) ) {count_second_longest[1][0]++;}
        else if( GeneralAnalysisHelper::GetRecoLengthWithPdg( 13 ) > GeneralAnalysisHelper::GetRecoLengthWithPdg( 2212 ) && GeneralAnalysisHelper::GetRecoLengthWithPdg( 2212 ) > GeneralAnalysisHelper::GetRecoLengthWithPdg( 211 ) ) {count_second_longest[1][2]++;} 
        else if( GeneralAnalysisHelper::GetRecoLengthWithPdg( 211 ) > GeneralAnalysisHelper::GetRecoLengthWithPdg( 2212 ) && GeneralAnalysisHelper::GetRecoLengthWithPdg( 2212 ) > GeneralAnalysisHelper::GetRecoLengthWithPdg( 13 ) ) {count_second_longest[1][2]++;}
          
      }
      if( e.CheckRecoTopology( GeneralAnalysisHelper::GetCC1PiTopologyMap() ) == 1 ) {
        if( GeneralAnalysisHelper::GetRecoLengthWithPdg( 13 ) > GeneralAnalysisHelper::GetRecoLengthWithPdg( 211 ) && GeneralAnalysisHelper::GetRecoLengthWithPdg( 211 ) > GeneralAnalysisHelper::GetRecoLengthWithPdg( 2212 ) ) {count_second_longest[2][1]++;}
        else if( GeneralAnalysisHelper::GetRecoLengthWithPdg( 2212 ) > GeneralAnalysisHelper::GetRecoLengthWithPdg( 211 ) && GeneralAnalysisHelper::GetRecoLengthWithPdg( 211 ) > GeneralAnalysisHelper::GetRecoLengthWithPdg( 13 ) ) {count_second_longest[2][1]++;}
        else if( GeneralAnalysisHelper::GetRecoLengthWithPdg( 211 ) > GeneralAnalysisHelper::GetRecoLengthWithPdg( 13 ) && GeneralAnalysisHelper::GetRecoLengthWithPdg( 13 ) > GeneralAnalysisHelper::GetRecoLengthWithPdg( 2212 ) ) {count_second_longest[2][0]++;}
        else if( GeneralAnalysisHelper::GetRecoLengthWithPdg( 2212 ) > GeneralAnalysisHelper::GetRecoLengthWithPdg( 13 ) && GeneralAnalysisHelper::GetRecoLengthWithPdg( 13 ) > GeneralAnalysisHelper::GetRecoLengthWithPdg( 211 ) ) {count_second_longest[2][0]++;}
        else if( GeneralAnalysisHelper::GetRecoLengthWithPdg( 13 ) > GeneralAnalysisHelper::GetRecoLengthWithPdg( 2212 ) && GeneralAnalysisHelper::GetRecoLengthWithPdg( 2212 ) > GeneralAnalysisHelper::GetRecoLengthWithPdg( 211 ) ) {count_second_longest[2][2]++;} 
        else if( GeneralAnalysisHelper::GetRecoLengthWithPdg( 211 ) > GeneralAnalysisHelper::GetRecoLengthWithPdg( 2212 ) && GeneralAnalysisHelper::GetRecoLengthWithPdg( 2212 ) > GeneralAnalysisHelper::GetRecoLengthWithPdg( 13 ) ) {count_second_longest[2][2]++;} 
      }
    }*/
    return count_longest;
    
  }
  
  //-----------------------------------------------------------------------------------------------

  float CC1piAnalysisHelper::GetMCQ2WithPdg(const Event &e, const int pdg) {
    //return -( GeneralAnalysisHelper::GetMCModulusMomentumWithPdg( 13 )-GeneralAnalysisHelper::GetMCModulusMomentumWithPdg( 211 ))*(GeneralAnalysisHelper::GetMCModulusMomentumWithPdg( 13 )-GeneralAnalysisHelper::GetMCModulusMomentumWithPdg( 211 ));
    return 0.;
  }

  //------------------------------------------------------------------------------------------------


  float CC1piAnalysisHelper::GetRecoQ2WithPdg(const Event &e, const int pdg) {
    //return -( GeneralAnalysisHelper::GetRecoModulusMomentumWithPdg( 13 )-GeneralAnalysisHelper::GetRecoModulusMomentumWithPdg( 211 ))*(GeneralAnalysisHelper::GetRecoModulusMomentumWithPdg( 13 )-GeneralAnalysisHelper::GetRecoModulusMomentumWithPdg( 211 ));
    return 0.;
  }

  //-----------------------------------------------------------------------------------------------

  float CC1piAnalysisHelper::GetEnergyLongest(const Event &e, const ParticleList &particle_list) {
 
    float MaxEnergy = 0.;
    float MaxLength = 0.;
    if( e.CheckMCTopology( GeneralAnalysisHelper::GetCC1PiTopologyMap() ) == 1  ){
      for(unsigned int i = 0; i < particle_list.size(); ++i) {
        if(MaxLength < particle_list[i].GetLength() && particle_list[i].GetPdgCode() != 111 ){
          MaxEnergy = particle_list[i].GetEnergy();
          MaxLength = particle_list[i].GetLength();
        }
      }
    }
    return MaxEnergy;
  }

  //-----------------------------------------------------------------------------------------------

  float CC1piAnalysisHelper::GetMCEnergyLongest(const Event &e) {
    return GetEnergyLongest(e, e.GetMCParticleList());
  }

  //-----------------------------------------------------------------------------------------------

  float CC1piAnalysisHelper::GetRecoEnergyLongest(const Event &e) {
    return GetEnergyLongest(e, e.GetRecoParticleList());
  }

  //-----------------------------------------------------------------------------------------------

  float CC1piAnalysisHelper::GetEnergySecondLongest(const Event &e, const ParticleList &particle_list) {
    float SecondMaxEnergy = 0.;
    float SecondMaxLength = 0;
    float MaxLength = 0.;
    if( e.CheckMCTopology( GeneralAnalysisHelper::GetCC1PiTopologyMap() ) == 1 ){
      for(unsigned int i = 0; i < particle_list.size(); ++i) {
        if(MaxLength < particle_list[i].GetLength() && particle_list[i].GetPdgCode() != 111 ){
          MaxLength = particle_list[i].GetLength();
        }
        else if (MaxLength > particle_list[i].GetLength() && SecondMaxLength < particle_list[i].GetLength() && particle_list[i].GetPdgCode() != 111){
          SecondMaxEnergy = particle_list[i].GetEnergy();
          SecondMaxLength = particle_list[i].GetLength();
        }
      }
      return SecondMaxEnergy;
    }
    return 0.;
  }

  //-----------------------------------------------------------------------------------------------

  float CC1piAnalysisHelper::GetMCEnergySecondLongest(const Event &e) {
    return GetEnergySecondLongest(e, e.GetMCParticleList());
  }

  //-----------------------------------------------------------------------------------------------

  float CC1piAnalysisHelper::GetRecoEnergySecondLongest(const Event &e) {
    return GetEnergySecondLongest(e, e.GetRecoParticleList());
  }
  
  //-----------------------------------------------------------------------------------------------

  float CC1piAnalysisHelper::GetKineticEnergyLongest(const Event &e, const ParticleList & particle_list) {
    float MaxKineticEnergy = 0.;
    float MaxLength = 0.;
    if( e.CheckMCTopology( GeneralAnalysisHelper::GetCC1PiTopologyMap() ) == 1  ){
      for(unsigned int i = 0; i < particle_list.size(); ++i) {
        if(MaxLength < particle_list[i].GetLength() && particle_list[i].GetPdgCode() != 111 ){
          MaxKineticEnergy = particle_list[i].GetKineticEnergy();
          MaxLength = particle_list[i].GetLength();
        }
      }
    }
    return MaxKineticEnergy;
  }

  //-----------------------------------------------------------------------------------------------

  float CC1piAnalysisHelper::GetMCKineticEnergyLongest(const Event &e) {
    return GetKineticEnergyLongest(e, e.GetMCParticleList());
  }

  //-----------------------------------------------------------------------------------------------

  float CC1piAnalysisHelper::GetRecoKineticEnergyLongest(const Event &e) {
    return GetKineticEnergyLongest(e, e.GetMCParticleList());
  }

  //-----------------------------------------------------------------------------------------------

  float CC1piAnalysisHelper::GetKineticEnergySecondLongest(const Event &e, const ParticleList & particle_list) {
    float SecondMaxKineticEnergy = 0.;
    float SecondMaxLength = 0;
    float MaxLength = 0.;
    if( e.CheckMCTopology( GeneralAnalysisHelper::GetCC1PiTopologyMap() ) == 1 ){
      for(unsigned int i = 0; i < particle_list.size(); ++i) {
        if(MaxLength < particle_list[i].GetLength() && particle_list[i].GetPdgCode() != 111 ){
          MaxLength = particle_list[i].GetLength();
        }
        else if ( MaxLength > particle_list[i].GetLength() && SecondMaxLength < particle_list[i].GetLength() && particle_list[i].GetPdgCode() != 111  ){
          SecondMaxKineticEnergy = particle_list[i].GetKineticEnergy();
          SecondMaxLength = particle_list[i].GetLength();
        }
      }
      return SecondMaxKineticEnergy;
    }
    return 0.;
  }

  //-----------------------------------------------------------------------------------------------

  float CC1piAnalysisHelper::GetMCKineticEnergySecondLongest(const Event &e) {
    return GetKineticEnergySecondLongest(e, e.GetMCParticleList());
  }

  //-----------------------------------------------------------------------------------------------

  float CC1piAnalysisHelper::GetRecoKineticEnergySecondLongest(const Event &e) {
    return GetKineticEnergySecondLongest(e, e.GetRecoParticleList());
  }

  //-----------------------------------------------------------------------------------------------

  float CC1piAnalysisHelper::GetModulusMomentumLongest(const Event &e, const ParticleList & particle_list) {
    float MaxMomentum = 0.;
    float MaxLength = 0.;
    if( e.CheckMCTopology( GeneralAnalysisHelper::GetCC1PiTopologyMap() ) == 1  ){
      for(unsigned int i = 0; i < particle_list.size(); ++i) {
        if(MaxLength < particle_list[i].GetLength() && particle_list[i].GetPdgCode() != 111 ){
          MaxMomentum = particle_list[i].GetModulusMomentum();
          MaxLength = particle_list[i].GetLength();
        }
      }
    }
    return MaxMomentum;
  }

  //-----------------------------------------------------------------------------------------------

  float CC1piAnalysisHelper::GetMCModulusMomentumLongest(const Event &e) {
    return GetModulusMomentumLongest(e, e.GetMCParticleList());
  }


  //-----------------------------------------------------------------------------------------------

  float CC1piAnalysisHelper::GetRecoModulusMomentumLongest(const Event &e) {
    return GetModulusMomentumLongest(e, e.GetRecoParticleList());
  }
  
  //-----------------------------------------------------------------------------------------------

  float CC1piAnalysisHelper::GetModulusMomentumSecondLongest(const Event &e, const ParticleList & particle_list) {
    float SecondMaxModulusMomentum = 0.;
    float SecondMaxLength = 0;
    float MaxLength = 0.;
    if( e.CheckMCTopology( GeneralAnalysisHelper::GetCC1PiTopologyMap() ) == 1 ){
      for(unsigned int i = 0; i < particle_list.size(); ++i) {
        if(MaxLength < particle_list[i].GetLength() && particle_list[i].GetPdgCode() != 111 ){
          MaxLength = particle_list[i].GetLength();
        }
        else if ( MaxLength > particle_list[i].GetLength() && SecondMaxLength < particle_list[i].GetLength() && particle_list[i].GetPdgCode() != 111  ){
          SecondMaxModulusMomentum = particle_list[i].GetModulusMomentum();
          SecondMaxLength = particle_list[i].GetLength();
        }
      }
      return SecondMaxModulusMomentum;
    }
    return 0.;
  }

  //-----------------------------------------------------------------------------------------------

  float CC1piAnalysisHelper::GetMCModulusMomentumSecondLongest(const Event &e) {
    return GetModulusMomentumSecondLongest(e, e.GetMCParticleList());
  }


  //-----------------------------------------------------------------------------------------------

  float CC1piAnalysisHelper::GetRecoModulusMomentumSecondLongest(const Event &e) {
    return GetModulusMomentumSecondLongest(e, e.GetRecoParticleList());
  }

  //-----------------------------------------------------------------------------------------------

  float CC1piAnalysisHelper::GetDeltaEnergy(const Event &e, const ParticleList &particle_list) {
    /*
     *
     * NUANCE CODE NO LONGER EXISTS 
     * SEE SCATTER CODES INSTEAD
     *
     *
     *
     * TLorentzVector pion, p2;
    for(unsigned int i = 0; i < particle_list.size(); ++i) {
	    if( e.GetNuanceCode() == 1003 || e.GetNuanceCode() == 1009 ){
	      if( particle_list[i].GetPdgCode() ==  211 ){
	        pion.SetE( particle_list[i].GetEnergy() );
	        pion.SetVect( particle_list[i].GetMomentum() );
	      }
        else if( particle_list[i].GetPdgCode() == 2212 ) {
	        p2.SetE( particle_list[i].GetEnergy() );
	        p2.SetVect( particle_list[i].GetMomentum() );
	      }
      } 
      else if( e.GetNuanceCode() == 1005 ){
	      if( particle_list[i].GetPdgCode() == 211 ){
	        pion.SetE( particle_list[i].GetEnergy() );
	        pion.SetVect( particle_list[i].GetMomentum() );
	      }
        else if( particle_list[i].GetPdgCode() == 2112 ) {
	        p2.SetE( particle_list[i].GetEnergy() );
	        p2.SetVect( particle_list[i].GetMomentum() );
	      }
      }
    }
    return (pion + p2).M();
    */
  }
  
  //---------------------------------------------------------------------------------------------                                                     

  float CC1piAnalysisHelper::GetDeltaEnergy_p(const ParticleList &particle_list) {
    TLorentzVector pion, p2;                                                                                          
    for(unsigned int i = 0; i < particle_list.size(); ++i) {
      if( particle_list[i].GetPdgCode() ==  211 ){
        pion.SetE( particle_list[i].GetEnergy() );
        pion.SetVect( particle_list[i].GetMomentum() );
      }
      else if( particle_list[i].GetPdgCode() == 2212 ) {
        p2.SetE( particle_list[i].GetEnergy() );
        p2.SetVect( particle_list[i].GetMomentum() );
      }
    }
    return (pion + p2).M();
  }

  //-----------------------------------------------------------------------------------------------
  
  float CC1piAnalysisHelper::GetMCDeltaEnergy(const Event &e) {
     return GetDeltaEnergy(e, e.GetMCParticleList());
  }

  //-----------------------------------------------------------------------------------------------

  float CC1piAnalysisHelper::GetCC1piNeutrinoEnergy(const ParticleList &particle_list) {
    // JLAB measured V on Argon - 0.0295 GeV
    // The variables from the branches and get the leaves
    float m_n   = 0.93957;   // Neutron mass, GeV
    float m_p   = 0.93828;   // Neutron mass, GeV
    float m_D   =   1.232;   // Delta mass, GeV
    float V     = 0.02950;   // Nucleon removal energy, GeV
    float m_mu  = 0.10566;   // Muon mass, GeV
    float reco=0, e, p, cth;   // track variables
    // Assuming D+ production ( xsec(n) bigger in Ar ). Modified m_n to consider the D++ production.
    m_n=(22*m_n + 18*m_p)/40.;
    // Vector of z direction
    TVector3 z;
    z[0] = 0;
    z[1] = 0;
    z[2] = 1;

    for(unsigned int i = 0; i < particle_list.size(); ++i) {
      if( particle_list[i].GetPdgCode() == 13 ){
        // Get the values needed
        e    = particle_list[i].GetEnergy();
        p    = particle_list[i].GetMomentum().Mag();
        cth  = (1/p) * (particle_list[i].GetMomentum()).Dot(z);
        
        reco = (1/(m_n - V - e + p*cth))*((m_n - V)*e - (m_mu*m_mu*0.5) + m_n*V - (V*V*0.5) + (m_D*m_D - m_n*m_n)*0.5);
        return reco;
      }
    }
    return reco;
  }

  //-------------------------------------------------------------------------------------------------

  float CC1piAnalysisHelper::GetCC1piNeutrinoEnergyMethod2(const ParticleList &particle_list) {
    // JLAB measured V on Argon - 0.0295 GeV
    // The variables from the branches and get the leaves
    float m_n   = 0.93957;   // Neutron mass, GeV
    float m_p   = 0.93828;   // Neutron mass, GeV
    float m_mu  = 0.10566;   // Muon mass, GeV
    float m_pi  = 0.13957;   // Pion mass, GeV
    float reco=0,p_mu, p_pi, e_mu, e_pi;   // track variables

    // Vector of z direction
    TVector3  z;
    z[0] = 0;
    z[1] = 0;
    z[2] = 1;
    double angle_mu, angle_pi; //cos angle_mu_pi;

    for(unsigned int i = 0; i < particle_list.size(); ++i) {
      for(unsigned int j = 0; j < particle_list.size(); ++j) {
        if( particle_list[i].GetPdgCode() == 13 ){
          e_mu    = particle_list[i].GetEnergy();
          p_mu    = particle_list[i].GetMomentum().Mag();
          angle_mu = particle_list[i].GetCosTheta();
          if( particle_list[j].GetPdgCode() == 211 ){
            e_pi    = particle_list[j].GetEnergy();
            p_pi    = particle_list[j].GetMomentum().Mag();
            angle_pi =  particle_list[j].GetCosTheta();
            reco = (m_mu*m_mu+m_pi*m_pi-2*m_p*(e_mu+e_pi)+2*p_mu*p_pi)/(2*(e_mu+e_pi-p_mu*angle_mu-p_pi*angle_pi-m_p));
            return reco;
          }
        }
      }
    }
    return reco;
  }

  //-------------------------------------------------------------------------------------------------

  float CC1piAnalysisHelper::GetMCCC1piNeutrinoEnergy(const Event &e) {
    return GetCC1piNeutrinoEnergy(e.GetMCParticleList());   
  }

  //-------------------------------------------------------------------------------------------------

  float CC1piAnalysisHelper::GetRecoCC1piNeutrinoEnergy(const Event &e) {
    return GetCC1piNeutrinoEnergy(e.GetRecoParticleList());   
  }

  //-------------------------------------------------------------------------------------------------

  float CC1piAnalysisHelper::GetRecoCC1piNeutrinoEnergyMethod2(const Event &e) {
    return GetCC1piNeutrinoEnergyMethod2(e.GetRecoParticleList());   
  }
  
  } // selection
