#include "SinglePhoton_module.h"

namespace single_photon
{
    void SinglePhoton::ClearMCTruths(){
        m_mctruth_num = 0;
        m_mctruth_origin = -99;
        m_mctruth_mode = -99;
        m_mctruth_interaction_type = -99;
        m_mctruth_nu_vertex_x = -9999;
        m_mctruth_nu_vertex_y = -9999;
        m_mctruth_nu_vertex_z = -9999;
        m_mctruth_ccnc = -99;
        m_mctruth_qsqr = -99;
        m_mctruth_nu_E = -99;
        m_mctruth_nu_pdg = 0;
        m_mctruth_lepton_pdg = 0;
        m_mctruth_num_daughter_particles = -99;
        m_mctruth_daughters_pdg.clear();
        m_mctruth_daughters_E.clear();

        m_mctruth_is_delta_radiative = 0;
        m_mctruth_delta_radiative_1g1p_or_1g1n = -999;

        m_mctruth_delta_photon_energy=-999;
        m_mctruth_delta_proton_energy=-999;
        m_mctruth_delta_neutron_energy=-999;

        m_mctruth_num_exiting_photons =0;
        m_mctruth_num_exiting_protons =0;
        m_mctruth_num_exiting_pi0 =0;
        m_mctruth_num_exiting_pipm =0;
        m_mctruth_num_exiting_neutrons=0;
        m_mctruth_num_exiting_delta0=0;
        m_mctruth_num_exiting_deltapm=0;
        m_mctruth_num_exiting_deltapp=0;

        m_mctruth_num_reconstructable_protons = 0;

        m_mctruth_is_reconstructable_1g1p = 0;
        m_mctruth_is_reconstructable_1g0p = 0;

        m_mctruth_leading_exiting_proton_energy = -9999;

        m_mctruth_exiting_pi0_E.clear();
        m_mctruth_exiting_pi0_px.clear();
        m_mctruth_exiting_pi0_py.clear();
        m_mctruth_exiting_pi0_pz.clear();

        m_mctruth_pi0_leading_photon_energy = -9999;
        m_mctruth_pi0_subleading_photon_energy = -9999;

        m_mctruth_exiting_delta0_num_daughters.clear();

        m_mctruth_exiting_photon_mother_trackID.clear();
        m_mctruth_exiting_photon_trackID.clear();
        m_mctruth_exiting_photon_from_delta_decay.clear();
        m_mctruth_exiting_photon_energy.clear();

        m_mctruth_exiting_proton_mother_trackID.clear();
        m_mctruth_exiting_proton_trackID.clear();
        m_mctruth_exiting_proton_from_delta_decay.clear();
        m_mctruth_exiting_proton_energy.clear();

    }

    void SinglePhoton::ResizeMCTruths(size_t size){
        m_mctruth_daughters_pdg.resize(size);
        m_mctruth_daughters_E.resize(size);

    }

    //add into vertex tree?
    void SinglePhoton::CreateMCTruthBranches(){
        vertex_tree->Branch("mctruth_num",&m_mctruth_num);
        vertex_tree->Branch("mctruth_origin",&m_mctruth_origin);
        vertex_tree->Branch("mctruth_nu_pdg",&m_mctruth_nu_pdg);
        vertex_tree->Branch("mctruth_nu_E",&m_mctruth_nu_E);

        vertex_tree->Branch("mctruth_nu_vertex_x",&m_mctruth_nu_vertex_x);
        vertex_tree->Branch("mctruth_nu_vertex_y",&m_mctruth_nu_vertex_y);
        vertex_tree->Branch("mctruth_nu_vertex_z",&m_mctruth_nu_vertex_z);

        vertex_tree->Branch("mctruth_lepton_pdg",&m_mctruth_lepton_pdg);
        vertex_tree->Branch("mctruth_lepton_E",&m_mctruth_lepton_E);
        vertex_tree->Branch("mctruth_mode",&m_mctruth_mode);
        vertex_tree->Branch("mctruth_qsqr",&m_mctruth_qsqr);
        vertex_tree->Branch("mctruth_cc_or_nc",&m_mctruth_ccnc);
        vertex_tree->Branch("mctruth_interaction_type",&m_mctruth_interaction_type);

        vertex_tree->Branch("mctruth_num_daughter_particles",&m_mctruth_num_daughter_particles);
        vertex_tree->Branch("mctruth_daughters_pdg",&m_mctruth_daughters_pdg);
        vertex_tree->Branch("mctruth_daughters_E",&m_mctruth_daughters_E);


        vertex_tree->Branch("mctruth_num_exiting_protons",&m_mctruth_num_exiting_protons);
        vertex_tree->Branch("mctruth_num_exiting_photons",&m_mctruth_num_exiting_photons);
        vertex_tree->Branch("mctruth_num_exiting_neutrons",&m_mctruth_num_exiting_neutrons);
        vertex_tree->Branch("mctruth_num_exiting_pi0",&m_mctruth_num_exiting_pi0);
        vertex_tree->Branch("mctruth_num_exiting_pipm",&m_mctruth_num_exiting_pipm);
        vertex_tree->Branch("mctruth_num_exiting_delta0",&m_mctruth_num_exiting_delta0);
        vertex_tree->Branch("mctruth_num_exiting_deltapm",&m_mctruth_num_exiting_deltapm);
        vertex_tree->Branch("mctruth_num_exiting_deltapp",&m_mctruth_num_exiting_deltapp);

        vertex_tree->Branch("mctruth_leading_exiting_proton_energy",&m_mctruth_leading_exiting_proton_energy);
        vertex_tree->Branch("mctruth_is_delta_radiative",&m_mctruth_is_delta_radiative);
        vertex_tree->Branch("mctruth_delta_radiative_1g1p_or_1g1n",&m_mctruth_delta_radiative_1g1p_or_1g1n);
        vertex_tree->Branch("mctruth_delta_photon_energy",&m_mctruth_delta_photon_energy);
        vertex_tree->Branch("mctruth_delta_proton_energy",&m_mctruth_delta_proton_energy);
        vertex_tree->Branch("mctruth_delta_neutron_energy",&m_mctruth_delta_neutron_energy);
        vertex_tree->Branch("mctruth_exiting_delta0_num_daughters",&m_mctruth_exiting_delta0_num_daughters);

        vertex_tree->Branch("mctruth_exiting_photon_trackID",&m_mctruth_exiting_photon_trackID);
        vertex_tree->Branch("mctruth_exiting_photon_mother_trackID",&m_mctruth_exiting_photon_mother_trackID);
        vertex_tree->Branch("mctruth_exiting_photon_from_delta_decay",&m_mctruth_exiting_photon_from_delta_decay);
        vertex_tree->Branch("mctruth_exiting_photon_energy",&m_mctruth_exiting_photon_energy);

        vertex_tree->Branch("mctruth_exiting_proton_trackID",&m_mctruth_exiting_proton_trackID);
        vertex_tree->Branch("mctruth_exiting_proton_mother_trackID",&m_mctruth_exiting_proton_mother_trackID);
        vertex_tree->Branch("mctruth_exiting_proton_from_delta_decay",&m_mctruth_exiting_proton_from_delta_decay);
        vertex_tree->Branch("mctruth_exiting_proton_energy",&m_mctruth_exiting_proton_energy);

        vertex_tree->Branch("mctruth_is_reconstructable_1g1p",&m_mctruth_is_reconstructable_1g1p);

        vertex_tree->Branch("mctruth_is_reconstructable_1g1p",&m_mctruth_is_reconstructable_1g1p);
        vertex_tree->Branch("mctruth_num_reconstructable_protons",&m_mctruth_num_reconstructable_protons);

        vertex_tree->Branch("mctruth_pi0_leading_photon_energy",&m_mctruth_pi0_leading_photon_energy);
        vertex_tree->Branch("mctruth_pi0_subleading_photon_energy",&m_mctruth_pi0_subleading_photon_energy);

        vertex_tree->Branch("mctruth_exiting_pi0_E",&m_mctruth_exiting_pi0_E);
        vertex_tree->Branch("mctruth_exiting_pi0_px",&m_mctruth_exiting_pi0_px);
        vertex_tree->Branch("mctruth_exiting_pi0_py",&m_mctruth_exiting_pi0_py);
        vertex_tree->Branch("mctruth_exiting_pi0_pz",&m_mctruth_exiting_pi0_pz);
    }

    void SinglePhoton::AnalyzeMCTruths(std::vector<art::Ptr<simb::MCTruth>> & mcTruthVector , std::vector<art::Ptr<simb::MCParticle>> & mcParticleVector){
        m_mctruth_num = mcTruthVector.size();
        this->ResizeMCTruths(m_mctruth_num);

        if(m_is_verbose) std::cout<<"SinglePhoton::AnalyzeMCTruths()\t||\t Starting to analyze "<<m_mctruth_num<<" simb::MCTruth's."<<std::endl;
        if(m_mctruth_num >1){
            std::cout<<"SinglePhoton::AnalyzeMCTruths()\t||\t WARNING There is more than 1 MCTruth neutrino interaction. Just running over the first simb::MCTruth."<<std::endl;
        }

        //one mctruth per event.  contains list of all particles 
        for(int i=0; i<std::min(1,m_mctruth_num); i++){
            const art::Ptr<simb::MCTruth> truth = mcTruthVector[i];

            m_mctruth_origin = truth->Origin();
            m_mctruth_ccnc = truth->GetNeutrino().CCNC();
            m_mctruth_mode = truth->GetNeutrino().Mode();
            m_mctruth_interaction_type = truth->GetNeutrino().InteractionType();
            m_mctruth_qsqr = truth->GetNeutrino().QSqr();
            m_mctruth_nu_pdg = truth->GetNeutrino().Nu().PdgCode();
            m_mctruth_nu_E = truth->GetNeutrino().Nu().E();
            m_mctruth_lepton_pdg = truth->GetNeutrino().Lepton().PdgCode();
            m_mctruth_lepton_E = truth->GetNeutrino().Lepton().E();

            std::vector<double> corrected(3);
            this->spacecharge_correction( truth->GetNeutrino().Lepton(),corrected);

            m_mctruth_nu_vertex_x = corrected[0];
            m_mctruth_nu_vertex_y = corrected[1];
            m_mctruth_nu_vertex_z = corrected[2];


            m_mctruth_num_daughter_particles = truth->NParticles(); //MCTruth_NParticles
            this->ResizeMCTruths(m_mctruth_num_daughter_particles);


            //some temp variables to see if its 1g1p or 1g1n
            int tmp_n_photons_from_delta = 0; 
            int tmp_n_protons_from_delta = 0; 
            int tmp_n_neutrons_from_delta = 0; 


            m_mctruth_leading_exiting_proton_energy = -9999;

            for(int j=0; j< m_mctruth_num_daughter_particles; j++){

                const simb::MCParticle par = truth->GetParticle(j);
                m_mctruth_daughters_pdg[j] = par.PdgCode();
                m_mctruth_daughters_E[j] = par.E();



                switch(m_mctruth_daughters_pdg[j]){
                    case(22):
                        {
                            if(par.StatusCode() == 1){
                                m_mctruth_num_exiting_photons++;
                                m_mctruth_exiting_photon_mother_trackID.push_back(par.Mother());
                                m_mctruth_exiting_photon_trackID.push_back(par.TrackId());
                                m_mctruth_exiting_photon_energy.push_back(par.E());

                            }
                         if(m_is_verbose)   std::cout<<"SinglePhoton::AnalyzeMCTruths()\t||\t Photon "<<par.PdgCode()<<" (id: "<<par.TrackId()<<") with mother trackID: "<<par.Mother()<<". Status Code: "<<par.StatusCode()<<" and photon energy "<<par.E()<<std::endl;

                            //if its mother is a delta with statuscode 3, and it has status code 14, then its the internal product of the delta decay.
                            const  simb::MCParticle mother = truth->GetParticle(par.Mother());
                            if((par.StatusCode()==1 || par.StatusCode()==14 ) && is_delta_map.count(mother.PdgCode())>0 && mother.StatusCode()==3){
                                m_mctruth_delta_photon_energy = par.E();
                                tmp_n_photons_from_delta ++;
                                m_mctruth_is_delta_radiative++;
                            }
                        }
                        break;
                    case(111):
                        {
                            // Make sure the pi0 actually exits the nucleus
                            if (par.StatusCode() == 1) {
                                m_mctruth_exiting_pi0_E.push_back(par.E());
                                m_mctruth_exiting_pi0_px.push_back(par.Px());
                                m_mctruth_exiting_pi0_py.push_back(par.Py());
                                m_mctruth_exiting_pi0_pz.push_back(par.Pz());
                                m_mctruth_num_exiting_pi0++;
                            }
                            break;
                        }
                    case(211):
                    case(-211):
                        if (par.StatusCode() == 1) {
                            m_mctruth_num_exiting_pipm++;
                        }
                        break;
                    case(2212):
                        {
                            if(par.StatusCode() == 1){
                                m_mctruth_num_exiting_protons++;
                                m_mctruth_exiting_proton_mother_trackID.push_back(par.Mother());
                                m_mctruth_exiting_proton_trackID.push_back(par.TrackId());
                                m_mctruth_exiting_proton_energy.push_back(par.E());
                            }


                         if(m_is_verbose)   std::cout<<"SingleProton::AnalyzeMCTruths()\t||\t Proton "<<par.PdgCode()<<" (id: "<<par.TrackId()<<") with mother trackID: "<<par.Mother()<<". Status Code: "<<par.StatusCode()<<" and proton energy "<<par.E()<<std::endl;


                            //if its mother is a delta with statuscode 3, and it has status code 14, then its the internal product of the delta decay.
                            const  simb::MCParticle mother = truth->GetParticle(par.Mother());
                            if(par.StatusCode()==14 && is_delta_map.count(mother.PdgCode())>0 && mother.StatusCode()==3){
                                m_mctruth_delta_proton_energy = par.E();
                                tmp_n_protons_from_delta ++;
                            }


                            break;
                        }
                    case(2112):
                        {

                            m_mctruth_num_exiting_neutrons++;
                      if(m_is_verbose)      std::cout<<"SingleProton::AnalyzeMCTruths()\t||\t Neutron "<<par.PdgCode()<<" (id: "<<par.TrackId()<<") with mother trackID: "<<par.Mother()<<". Status Code: "<<par.StatusCode()<<" and neutron energy "<<par.E()<<std::endl;

                            //if its mother is a delta with statuscode 3, and it has status code 14, then its the internal product of the delta decay.
                            const  simb::MCParticle mother = truth->GetParticle(par.Mother());
                            if(par.StatusCode()==14 && is_delta_map.count(mother.PdgCode())>0 && mother.StatusCode()==3){
                                m_mctruth_delta_neutron_energy = par.E();
                                tmp_n_neutrons_from_delta ++;
                            }
                        }

                        break;
                    case(-2224):
                    case(2224):
                        if(par.StatusCode() == 1){  m_mctruth_num_exiting_deltapp++; }
                        break;
                    case(-2214):
                    case(2214):
                    case(-1114):
                    case(1114):
                        if(par.StatusCode() == 1){ m_mctruth_num_exiting_deltapm++; }
                        break;
                    case(-2114):
                    case(2114):
                        if(par.StatusCode() == 1){ 
                            m_mctruth_num_exiting_delta0++;
                            m_mctruth_exiting_delta0_num_daughters.push_back(par.NumberDaughters());
                         if(m_is_verbose)   std::cout<<"SinglePhoton::AnalyzeMCTruths()\t||\t Delta0 "<<par.PdgCode()<<" (id: "<<par.TrackId()<<") with "<<m_mctruth_exiting_delta0_num_daughters.back()<<" daughters. StatusCode "<<par.StatusCode()<<std::endl;
                        }
                        break;
                    default:
                        break;
                }
            }


            for(size_t p=0; p< m_mctruth_exiting_proton_energy.size(); p++){
                if( m_mctruth_exiting_proton_energy[p] > m_mctruth_leading_exiting_proton_energy ){
                    m_mctruth_leading_exiting_proton_energy = m_mctruth_exiting_proton_energy[p];
                }
            }



              std::cout<<"SinglePhoton::AnalyzeMCTruths()\t||\t This event is ";
            if(tmp_n_photons_from_delta==1 && tmp_n_protons_from_delta==1){
                m_mctruth_delta_radiative_1g1p_or_1g1n = 1;
                std::cout<<"a 1g1p delta radiative event"<<std::endl;
            }else if(tmp_n_photons_from_delta==1 && tmp_n_neutrons_from_delta==1){
                m_mctruth_delta_radiative_1g1p_or_1g1n = 0;
                std::cout<<"a 1g1n delta radiative event"<<std::endl;
                std::cout<<"SinglePhoton::AnalyzeMCTruths()\t||\t This event is  ";
            }else{
                std::cout<<"NOT a 1g1p or 1g1n delta radiative decay"<<std::endl;;
            }

            //Now for FSI exiting particles!
            m_mctruth_exiting_photon_from_delta_decay.resize(m_mctruth_num_exiting_photons,0);
            m_mctruth_exiting_proton_from_delta_decay.resize(m_mctruth_num_exiting_protons,0);


            //second loop for some dauhter info
            // status codes!
            // 0 initial state
            // 1 stable final state
            // 2 intermediate state
            // 3 decayed state
            // 11 Nucleon target
            // 14 hadron in the nucleas 

            // So if a  final_state_particle has a status(3) delta in its history its "from" a delta.
            //first we loop over all 14's to see which have a direct mother delta. [done above]
            //so first we loop over all state 1 (exiting) to see what a LArTPC sees (post FSI)
                for (unsigned int p = 0; p <  m_mctruth_exiting_photon_energy.size(); p++){

            if ( m_mctruth_exiting_photon_energy[p] > m_exiting_photon_energy_threshold){
                        m_mctruth_num_reconstructable_protons++;

            }//if g above threshold
                }

            //if it's a true delta radiative event, check the energies



            if (m_mctruth_is_delta_radiative==true){
                for (unsigned int p = 0; p <  m_mctruth_exiting_photon_energy.size(); p++){
                    std::cout<<"SinglePhoton::AnalyzeMCTruths()\t||\tLooking at exiting photon with energy "<<m_mctruth_exiting_photon_energy[p]<<std::endl;
                    if ( m_mctruth_exiting_photon_energy[p] > m_exiting_photon_energy_threshold){
                        m_mctruth_is_reconstructable_1g0p = true;

                    }//if g above threshold
                }//for all exiting g
                for(unsigned int pr = 0; pr <  m_mctruth_exiting_proton_energy.size(); pr++){
                   if ( m_mctruth_exiting_proton_energy[pr]> m_exiting_proton_energy_threshold){
                        //if it's already 1g1p then we've found a 1g2p which we aren't counting
                        if( m_mctruth_is_reconstructable_1g1p == true && m_mctruth_is_reconstructable_1g0p == false){
                            m_mctruth_is_reconstructable_1g1p = false;
                        }
                        //if there's a photon then it's actually a 1g1p
                        if( m_mctruth_is_reconstructable_1g0p == true &&  m_mctruth_is_reconstructable_1g1p == false){
                            m_mctruth_is_reconstructable_1g1p = true;
                            m_mctruth_is_reconstructable_1g0p = false;
                        } 
                       std::cout<<"SinglePhoton::AnalyzeMCTruths()\t||\tChecking proton with energy "<<m_mctruth_exiting_proton_energy[pr]<<", is 1g1p/1g0p= "<< m_mctruth_is_reconstructable_1g1p<<"/"<< m_mctruth_is_reconstructable_1g0p<<std::endl;
                         }//if p above threshold
                }//for all exiting p

            }//if ncdelta


            //So for all photons that have status code 1 i.e all exiting ones...
            for(int p =0; p < m_mctruth_num_exiting_photons; ++p){
                const simb::MCParticle mother = truth->GetParticle(m_mctruth_exiting_photon_mother_trackID[p]);

                std::cout<<"SinglePhoton::AnalyzeMCTruths()\t||\t -- gamma ("<<m_mctruth_exiting_photon_trackID[p]<<") of status_code 1.. "<<std::endl;
                std::cout<<"SinglePhoton::AnalyzeMCTruths()\t||\t ---- with mother "<<mother.PdgCode()<<" ("<<m_mctruth_exiting_photon_mother_trackID[p]<<") status_code "<<mother.StatusCode()<<std::endl;
                simb::MCParticle nth_mother = mother;
                int n_generation = 2;

                while(nth_mother.StatusCode() != 0){

                    nth_mother = truth->GetParticle(nth_mother.Mother()); 
                    std::cout<<"SinglePhoton::AnalyzeMCTruths()\t||\t ---- and "<<n_generation<<"-mother "<<nth_mother.PdgCode()<<" ("<<nth_mother.TrackId()<<") and status_code "<<nth_mother.StatusCode()<<std::endl;
                    if( is_delta_map.count(nth_mother.PdgCode())>0 && nth_mother.StatusCode()==3){
                        std::cout<<"SinglePhoton::AnalyzeMCTruths()\t||\t ------ Defintely From a Delta Decay! : "<<is_delta_map[nth_mother.PdgCode()]<<std::endl;
                        m_mctruth_exiting_photon_from_delta_decay[p] = 1;
                    }
                    n_generation++;
                }
            }

            //So for all protons that have status code 1 i.e all exiting ones...
            for(int p =0; p < m_mctruth_num_exiting_protons; ++p){
                const simb::MCParticle mother = truth->GetParticle(m_mctruth_exiting_proton_mother_trackID[p]);

                if(m_is_verbose){
                    std::cout<<"SinglePhoton::AnalyzeMCTruths()\t||\t -- proton ("<<m_mctruth_exiting_proton_trackID[p]<<") of status_code 1.. "<<std::endl;
                    std::cout<<"SinglePhoton::AnalyzeMCTruths()\t||\t ---- with mother "<<mother.PdgCode()<<" ("<<m_mctruth_exiting_proton_mother_trackID[p]<<") status_code "<<mother.StatusCode()<<std::endl;
                }
                simb::MCParticle nth_mother = mother;
                int n_generation = 2;

                while(nth_mother.StatusCode() != 0){

                    nth_mother = truth->GetParticle(nth_mother.Mother()); 
                    if(m_is_verbose) std::cout<<"SinglePhoton::AnalyzeMCTruths()\t||\t ---- and "<<n_generation<<"-mother "<<nth_mother.PdgCode()<<" ("<<nth_mother.TrackId()<<") and status_code "<<nth_mother.StatusCode()<<std::endl;
                    if(is_delta_map.count(nth_mother.PdgCode())>0 && nth_mother.StatusCode()==3){
                       if(m_is_verbose) std::cout<<"SinglePhoton::AnalyzeMCTruths()\t||\t ------ Defintely From a Delta Decay! : "<<is_delta_map[nth_mother.PdgCode()]<<std::endl;
                        m_mctruth_exiting_proton_from_delta_decay[p] = 1;
                    } 
                    n_generation++;
                }


            }





            if(m_is_verbose){
                std::cout<<"SinglePhoton::AnalyzeMCTruths()\t||\t This is a CCNC: "<<m_mctruth_ccnc<<" event with a nu_pdg: "<<m_mctruth_nu_pdg<<" and "<<m_mctruth_num_daughter_particles<<" exiting particles."<<std::endl;
                std::cout<<"SinglePhoton::AnalyzeMCTruths()\t||\t With  "<<m_mctruth_num_exiting_pi0<<" Pi0, "<<m_mctruth_num_exiting_pipm<<" Pi+/-, "<<m_mctruth_num_exiting_protons<<" Protons, "<<m_mctruth_num_exiting_neutrons<<" neutrons and "<<m_mctruth_num_exiting_delta0<<" delta0, "<<m_mctruth_num_exiting_deltapm<<" deltapm, "<<m_mctruth_num_exiting_deltapp<<" Deltas++"<<std::endl;
            }

        }// end of MCtruth loo

        //make a stupid temp map
        std::map<size_t,size_t> mymap;
        for(size_t k = 0; k < mcParticleVector.size(); k++){
            const art::Ptr<simb::MCParticle> mcp = mcParticleVector[k];
            mymap[mcp->TrackId()]       = k;
        }


        //Just some VERY hacky pi^0 photon stuff
        int npi0check = 0;
        for(size_t k = 0; k < mcParticleVector.size(); k++){
            const art::Ptr<simb::MCParticle> mcp = mcParticleVector[k];

            if(false) std::cout << k << " Mother:"<< mcp->Mother() << " pdgcode: " << mcp->PdgCode() << " trkid: " << mcp->TrackId() << " statuscode: " << mcp->StatusCode() << std::endl;


            if(mcp->PdgCode() == 111 && mcp->Mother() == 0 && mcp->NumberDaughters()==2 ){
                npi0check++;
                const art::Ptr<simb::MCParticle> dau1 = mcParticleVector[mymap[mcp->Daughter(0)]];
                const art::Ptr<simb::MCParticle> dau2 = mcParticleVector[mymap[mcp->Daughter(1)]];

                if(false)  std::cout<<"On Dau1: "<<" Mother:"<< dau1->Mother()<<" pdgcode: "<<dau1->PdgCode()<<" trkid: "<<dau1->TrackId()<<" statuscode: "<<dau1->StatusCode()<<std::endl;
                if(false)  std::cout<<"On Dau2: "<<" Mother:"<< dau2->Mother()<<" pdgcode: "<<dau2->PdgCode()<<" trkid: "<<dau2->TrackId()<<" statuscode: "<<dau2->StatusCode()<<std::endl;

                double e1 = dau1->E();
                double e2 = dau2->E();
                if(e2<e1){
                    m_mctruth_pi0_leading_photon_energy = e1;
                    m_mctruth_pi0_subleading_photon_energy = e2;
                }else{
                    m_mctruth_pi0_leading_photon_energy = e2;
                    m_mctruth_pi0_subleading_photon_energy = e1;
                }

            }
        }

        if(npi0check>1) std::cout<<"WARNING WARNING!!!! there are "<<npi0check<<" Pi0's in this event in geant4 that come from the nucleas"<<std::endl;

    }//end of analyze this


}
