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
        m_mctruth_reco_vertex_dist = -9999;
        m_mctruth_ccnc = -99;
        m_mctruth_qsqr = -99;
        m_mctruth_nu_E = -99;
        m_mctruth_nu_pdg = 0;
        m_mctruth_lepton_pdg = 0;
        m_mctruth_num_daughter_particles = -99;
        m_mctruth_daughters_pdg.clear();
        m_mctruth_daughters_E.clear();

        m_mctruth_daughters_status_code.clear();
        m_mctruth_daughters_trackID.clear();
        m_mctruth_daughters_mother_trackID.clear();
        m_mctruth_daughters_px.clear();
        m_mctruth_daughters_py.clear();
        m_mctruth_daughters_pz.clear();
        m_mctruth_daughters_startx.clear();
        m_mctruth_daughters_starty.clear();
        m_mctruth_daughters_startz.clear();
        m_mctruth_daughters_time.clear();
        m_mctruth_daughters_endx.clear();
        m_mctruth_daughters_endy.clear();
        m_mctruth_daughters_endz.clear();
        m_mctruth_daughters_endtime.clear();
        m_mctruth_daughters_process.clear();
        m_mctruth_daughters_end_process.clear();


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
        m_mctruth_exiting_pi0_mom.clear();
        m_mctruth_exiting_pi0_px.clear();
        m_mctruth_exiting_pi0_py.clear();
        m_mctruth_exiting_pi0_pz.clear();

        m_mctruth_pi0_leading_photon_energy = -9999;
        m_mctruth_pi0_subleading_photon_energy = -9999;
        m_mctruth_pi0_leading_photon_end_process = "none";
        m_mctruth_pi0_subleading_photon_end_process = "none";
        m_mctruth_pi0_leading_photon_end = {-9999,-9999,-9999};
        m_mctruth_pi0_leading_photon_start = {-9999,-9999,-9999};
        m_mctruth_pi0_subleading_photon_end = {-9999,-9999,-9999};
        m_mctruth_pi0_subleading_photon_start = {-9999,-9999,-9999};
        m_mctruth_pi0_leading_photon_exiting_TPC = -999;
        m_mctruth_pi0_subleading_photon_exiting_TPC = -999;
        m_mctruth_pi0_leading_photon_mom = {-9999,-9999,-9999};
        m_mctruth_pi0_subleading_photon_mom = {-9999,-9999,-9999};

        m_mctruth_exiting_delta0_num_daughters.clear();

        m_mctruth_exiting_photon_mother_trackID.clear();
        m_mctruth_exiting_photon_trackID.clear();
        m_mctruth_exiting_photon_from_delta_decay.clear();
        m_mctruth_exiting_photon_energy.clear();
        m_mctruth_exiting_photon_px.clear();
        m_mctruth_exiting_photon_py.clear();
        m_mctruth_exiting_photon_pz.clear();

        m_mctruth_exiting_proton_mother_trackID.clear();
        m_mctruth_exiting_proton_trackID.clear();
        m_mctruth_exiting_proton_from_delta_decay.clear();
        m_mctruth_exiting_proton_energy.clear();
        m_mctruth_exiting_proton_px.clear();
        m_mctruth_exiting_proton_py.clear();
        m_mctruth_exiting_proton_pz.clear();

        m_mctruth_exiting_neutron_mother_trackID.clear();
        m_mctruth_exiting_neutron_trackID.clear();
        m_mctruth_exiting_neutron_from_delta_decay.clear();
        m_mctruth_exiting_neutron_energy.clear();
        m_mctruth_exiting_neutron_px.clear();
        m_mctruth_exiting_neutron_py.clear();
        m_mctruth_exiting_neutron_pz.clear();

    }

    void SinglePhoton::ResizeMCTruths(size_t size){
        m_mctruth_daughters_pdg.resize(size);
        m_mctruth_daughters_E.resize(size);
        m_mctruth_daughters_status_code.resize(size);
        m_mctruth_daughters_trackID.resize(size);
        m_mctruth_daughters_mother_trackID.resize(size);
        m_mctruth_daughters_px.resize(size);
        m_mctruth_daughters_py.resize(size);
        m_mctruth_daughters_pz.resize(size);
        m_mctruth_daughters_startx.resize(size);
        m_mctruth_daughters_starty.resize(size);
        m_mctruth_daughters_startz.resize(size);
        m_mctruth_daughters_time.resize(size);
        m_mctruth_daughters_endx.resize(size);
        m_mctruth_daughters_endy.resize(size);
        m_mctruth_daughters_endz.resize(size);
        m_mctruth_daughters_endtime.resize(size);
        m_mctruth_daughters_end_process.resize(size);
        m_mctruth_daughters_process.resize(size);

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
        vertex_tree->Branch("mctruth_reco_vertex_dist",&m_mctruth_reco_vertex_dist);

        vertex_tree->Branch("mctruth_lepton_pdg",&m_mctruth_lepton_pdg);
        vertex_tree->Branch("mctruth_lepton_E",&m_mctruth_lepton_E);
        vertex_tree->Branch("mctruth_mode",&m_mctruth_mode);
        vertex_tree->Branch("mctruth_qsqr",&m_mctruth_qsqr);
        vertex_tree->Branch("mctruth_cc_or_nc",&m_mctruth_ccnc);
        vertex_tree->Branch("mctruth_interaction_type",&m_mctruth_interaction_type);

        vertex_tree->Branch("mctruth_num_daughter_particles",&m_mctruth_num_daughter_particles);
        vertex_tree->Branch("mctruth_daughters_pdg",&m_mctruth_daughters_pdg);
        vertex_tree->Branch("mctruth_daughters_E",&m_mctruth_daughters_E);
        vertex_tree->Branch("mctruth_daughters_status_code",&m_mctruth_daughters_status_code);
        vertex_tree->Branch("mctruth_daughters_trackID",&m_mctruth_daughters_trackID);
        vertex_tree->Branch("mctruth_daughters_mother_trackID",&m_mctruth_daughters_mother_trackID);
        vertex_tree->Branch("mctruth_daughters_px",&m_mctruth_daughters_px);
        vertex_tree->Branch("mctruth_daughters_py",&m_mctruth_daughters_py);
        vertex_tree->Branch("mctruth_daughters_pz",&m_mctruth_daughters_pz);
        vertex_tree->Branch("mctruth_daughters_startx",&m_mctruth_daughters_startx);
        vertex_tree->Branch("mctruth_daughters_starty",&m_mctruth_daughters_starty);
        vertex_tree->Branch("mctruth_daughters_startz",&m_mctruth_daughters_startz);
        vertex_tree->Branch("mctruth_daughters_time",&m_mctruth_daughters_time);
        vertex_tree->Branch("mctruth_daughters_endx",&m_mctruth_daughters_endx);
        vertex_tree->Branch("mctruth_daughters_endy",&m_mctruth_daughters_endy);
        vertex_tree->Branch("mctruth_daughters_endz",&m_mctruth_daughters_endz);
        vertex_tree->Branch("mctruth_daughters_endtime",&m_mctruth_daughters_endtime);
        vertex_tree->Branch("mctruth_daughters_process",&m_mctruth_daughters_process);
        vertex_tree->Branch("mctruth_daughters_end_process",&m_mctruth_daughters_end_process);




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
        vertex_tree->Branch("mctruth_exiting_photon_px",&m_mctruth_exiting_photon_px);
        vertex_tree->Branch("mctruth_exiting_photon_py",&m_mctruth_exiting_photon_py);
        vertex_tree->Branch("mctruth_exiting_photon_pz",&m_mctruth_exiting_photon_pz);

        vertex_tree->Branch("mctruth_exiting_proton_trackID",&m_mctruth_exiting_proton_trackID);
        vertex_tree->Branch("mctruth_exiting_proton_mother_trackID",&m_mctruth_exiting_proton_mother_trackID);
        vertex_tree->Branch("mctruth_exiting_proton_from_delta_decay",&m_mctruth_exiting_proton_from_delta_decay);
        vertex_tree->Branch("mctruth_exiting_proton_energy",&m_mctruth_exiting_proton_energy);
        vertex_tree->Branch("mctruth_exiting_proton_px",&m_mctruth_exiting_proton_px);
        vertex_tree->Branch("mctruth_exiting_proton_py",&m_mctruth_exiting_proton_py);
        vertex_tree->Branch("mctruth_exiting_proton_pz",&m_mctruth_exiting_proton_pz);

        vertex_tree->Branch("mctruth_exiting_neutron_trackID",&m_mctruth_exiting_neutron_trackID);
        vertex_tree->Branch("mctruth_exiting_neutron_mother_trackID",&m_mctruth_exiting_neutron_mother_trackID);
        vertex_tree->Branch("mctruth_exiting_neutron_from_delta_decay",&m_mctruth_exiting_neutron_from_delta_decay);
        vertex_tree->Branch("mctruth_exiting_neutron_energy",&m_mctruth_exiting_neutron_energy);
        vertex_tree->Branch("mctruth_exiting_neutron_px",&m_mctruth_exiting_neutron_px);
        vertex_tree->Branch("mctruth_exiting_neutron_py",&m_mctruth_exiting_neutron_py);
        vertex_tree->Branch("mctruth_exiting_neutron_pz",&m_mctruth_exiting_neutron_pz);


        vertex_tree->Branch("mctruth_is_reconstructable_1g1p",&m_mctruth_is_reconstructable_1g1p);

        vertex_tree->Branch("mctruth_is_reconstructable_1g1p",&m_mctruth_is_reconstructable_1g1p);
        vertex_tree->Branch("mctruth_num_reconstructable_protons",&m_mctruth_num_reconstructable_protons);

        vertex_tree->Branch("mctruth_pi0_leading_photon_energy",&m_mctruth_pi0_leading_photon_energy);
        vertex_tree->Branch("mctruth_pi0_leading_photon_mom",&m_mctruth_pi0_leading_photon_mom);
        vertex_tree->Branch("mctruth_pi0_leading_photon_start",&m_mctruth_pi0_leading_photon_start);
        vertex_tree->Branch("mctruth_pi0_leading_photon_end",&m_mctruth_pi0_leading_photon_end);
        vertex_tree->Branch("mctruth_pi0_leading_photon_exiting_TPC",&m_mctruth_pi0_leading_photon_exiting_TPC);
        vertex_tree->Branch("mctruth_pi0_subleading_photon_energy",&m_mctruth_pi0_subleading_photon_energy);
        vertex_tree->Branch("mctruth_pi0_subleading_photon_mom",&m_mctruth_pi0_subleading_photon_mom);
        vertex_tree->Branch("mctruth_pi0_subleading_photon_end_process",&m_mctruth_pi0_subleading_photon_end_process);
        vertex_tree->Branch("mctruth_pi0_subleading_photon_start",&m_mctruth_pi0_subleading_photon_start);
        vertex_tree->Branch("mctruth_pi0_subleading_photon_end",&m_mctruth_pi0_subleading_photon_end);
        vertex_tree->Branch("mctruth_pi0_subleading_photon_exiting_TPC",&m_mctruth_pi0_subleading_photon_exiting_TPC);


        vertex_tree->Branch("mctruth_exiting_pi0_E",&m_mctruth_exiting_pi0_E);
        vertex_tree->Branch("mctruth_exiting_pi0_mom",&m_mctruth_exiting_pi0_mom);
        vertex_tree->Branch("mctruth_exiting_pi0_px",&m_mctruth_exiting_pi0_px);
        vertex_tree->Branch("mctruth_exiting_pi0_py",&m_mctruth_exiting_pi0_py);
        vertex_tree->Branch("mctruth_exiting_pi0_pz",&m_mctruth_exiting_pi0_pz);
    }

	//CHECK how is this different from MCReco matching?
    void SinglePhoton::AnalyzeMCTruths(std::vector<art::Ptr<simb::MCTruth>> & mcTruthVector , std::vector<art::Ptr<simb::MCParticle>> & mcParticleVector){
        m_mctruth_num = mcTruthVector.size();
        if(m_is_verbose) std::cout<<"SinglePhoton::AnalyzeMCTruths()\t||\tTotal  "<<m_mctruth_num<<" simb::MCTruth's."<<std::endl;
        if(m_mctruth_num >1){
            std::cout<<"SinglePhoton::AnalyzeMCTruths()\t||\t WARNING There is more than 1 MCTruth neutrino interaction. Just running over the first simb::MCTruth."<<std::endl;
        }else if(m_mctruth_num==0){
            std::cout<<"SinglePhoton::AnalyzeMCTruths()\t||\t WARNING There is 0 MCTruth neutrino interaction. Break simb::MCTruth."<<std::endl;
        }

        //one mctruth per event.  contains list of all particles 
        for(int i=0; i<std::min(1,m_mctruth_num); i++){
            const art::Ptr<simb::MCTruth> truth = mcTruthVector[i];
//			std::cout<<"\nCHECK THIS MCTruth!! "<<*truth<<std::endl;


            m_mctruth_origin = truth->Origin();
//            if(m_is_verbose) std::cout<<"Getting origin "<<truth->Origin()<<std::endl;

            if(!truth->NeutrinoSet()){
                if(m_is_verbose) std::cout<<"Warning, no neutrino set skipping. "<<std::endl;
            }else{
//                if(m_is_verbose) std::cout<<"Getting origin "<<truth->Origin()<<std::endl;
                m_mctruth_ccnc = truth->GetNeutrino().CCNC();
//                if(m_is_verbose) std::cout<<"Getting ccnc "<<truth->GetNeutrino().CCNC()<<std::endl;
                m_mctruth_mode = truth->GetNeutrino().Mode();
//                if(m_is_verbose) std::cout<<"Getting Mode"<<std::endl;
                m_mctruth_interaction_type = truth->GetNeutrino().InteractionType();
//                if(m_is_verbose) std::cout<<"Getting Type"<<std::endl;
                m_mctruth_qsqr = truth->GetNeutrino().QSqr();
//                if(m_is_verbose) std::cout<<"Getting Q"<<std::endl;
                m_mctruth_nu_pdg = truth->GetNeutrino().Nu().PdgCode();
//                if(m_is_verbose) std::cout<<"Getting E"<<std::endl;
                m_mctruth_nu_E = truth->GetNeutrino().Nu().E();
//                if(m_is_verbose) std::cout<<"Getting pdg"<<std::endl;
                m_mctruth_lepton_pdg = truth->GetNeutrino().Lepton().PdgCode();
//                if(m_is_verbose) std::cout<<"Getting pdg lepton"<<std::endl;
                m_mctruth_lepton_E = truth->GetNeutrino().Lepton().E();
//                if(m_is_verbose) std::cout<<"Getting lepton E"<<std::endl;
            
//                if(m_is_verbose) std::cout<<"Getting SC corrected vertex position"<<std::endl;
                std::vector<double> corrected(3);
	        // get corrected lepton position
			// CHECK, turn simb::mcparticle to art::Ptr<simb::MCParticle, then it can remove one unnecessary spacecharge_correction() function.
                this->spacecharge_correction( truth->GetNeutrino().Lepton(),corrected);

                m_mctruth_nu_vertex_x = corrected[0];
                m_mctruth_nu_vertex_y = corrected[1];
                m_mctruth_nu_vertex_z = corrected[2];
                m_mctruth_reco_vertex_dist = sqrt(pow (m_mctruth_nu_vertex_x-m_vertex_pos_x,2)+pow (m_mctruth_nu_vertex_y-m_vertex_pos_y,2)+pow (m_mctruth_nu_vertex_z-m_vertex_pos_z,2));
            
            }



            m_mctruth_num_daughter_particles = truth->NParticles(); //MCTruth_NParticles
			if(m_is_verbose){
				std::cout<<"We are working with :"<<m_mctruth_num_daughter_particles<<" m_mctruth_num_daugher_particles "<<std::endl;
			}
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

                m_mctruth_daughters_status_code[j] = par.StatusCode();
                m_mctruth_daughters_trackID[j] = par.TrackId();
                m_mctruth_daughters_mother_trackID[j] = par.Mother();
                m_mctruth_daughters_px[j] = par.Px();
                m_mctruth_daughters_py[j] = par.Py();
                m_mctruth_daughters_pz[j] = par.Pz();
                m_mctruth_daughters_startx[j] = par.Vx();
                m_mctruth_daughters_starty[j] = par.Vy();
                m_mctruth_daughters_startz[j] = par.Vz();
                m_mctruth_daughters_time[j] = par.T();
                m_mctruth_daughters_endx[j] = par.EndX();
                m_mctruth_daughters_endy[j] = par.EndY();
                m_mctruth_daughters_endz[j] = par.EndZ();
                m_mctruth_daughters_endtime[j] = par.EndT();
                m_mctruth_daughters_process[j] = par.Process();  //Process() and EndProcess() return string
                m_mctruth_daughters_end_process[j] = par.EndProcess();

                if(m_is_textgen) continue; //quick hack, fix in files

                switch(m_mctruth_daughters_pdg[j]){
                    case(22): // if it's a gamma
                        {
                            if(par.StatusCode() == 1){
                                m_mctruth_num_exiting_photons++;
                                m_mctruth_exiting_photon_mother_trackID.push_back(par.Mother());
                                m_mctruth_exiting_photon_trackID.push_back(par.TrackId());
                                m_mctruth_exiting_photon_energy.push_back(par.E());
                                m_mctruth_exiting_photon_px.push_back(par.Px());
                                m_mctruth_exiting_photon_py.push_back(par.Py());
                                m_mctruth_exiting_photon_pz.push_back(par.Pz());
                            }
                         if(m_is_verbose)   std::cout<<"SinglePhoton::AnalyzeMCTruths()\t||\t Photon "<<par.PdgCode()<<" (id: "<<par.TrackId()<<") with mother trackID: "<<par.Mother()<<". Status Code: "<<par.StatusCode()<<" and photon energy "<<par.E()<<std::endl;

                            //if its mother is a delta with statuscode 3, and it has status code 14, then its the internal product of the delta decay.
                            if((par.StatusCode()==1 || par.StatusCode()==14 )){
                                const  simb::MCParticle mother = truth->GetParticle(par.Mother());
                                
                                if(is_delta_map.count(mother.PdgCode())>0 && mother.StatusCode()==3){
                                  m_mctruth_delta_photon_energy = par.E();
                                     tmp_n_photons_from_delta ++;
                                    m_mctruth_is_delta_radiative++;
                                }
                            }
                        }
                        break;
                    case(111): // if it's a pi0
                        {
                            // Make sure the pi0 actually exits the nucleus
                            if (par.StatusCode() == 1) {
                                m_mctruth_exiting_pi0_E.push_back(par.E());
                                m_mctruth_exiting_pi0_mom.push_back(sqrt(pow(par.Px(),2)+pow(par.Py(),2)+pow(par.Pz(),2)));
                                m_mctruth_exiting_pi0_px.push_back(par.Px());
                                m_mctruth_exiting_pi0_py.push_back(par.Py());
                                m_mctruth_exiting_pi0_pz.push_back(par.Pz());
                                m_mctruth_num_exiting_pi0++;
                            }
                            break;
                        }
                    case(211):
                    case(-211):  // it's pi+ or pi-
                        if (par.StatusCode() == 1) {
                            m_mctruth_num_exiting_pipm++;
                        }
                        break;
                    case(2212):  // if it's a proton
                        {
                            if(par.StatusCode() == 1){
                                m_mctruth_num_exiting_protons++;
                                m_mctruth_exiting_proton_mother_trackID.push_back(par.Mother());
                                m_mctruth_exiting_proton_trackID.push_back(par.TrackId());
                                m_mctruth_exiting_proton_energy.push_back(par.E());
                                m_mctruth_exiting_proton_px.push_back(par.Px());
                                m_mctruth_exiting_proton_py.push_back(par.Py());
                                m_mctruth_exiting_proton_pz.push_back(par.Pz());
                            }


                         if(m_is_verbose)   std::cout<<"SingleProton::AnalyzeMCTruths()\t||\t Proton "<<par.PdgCode()<<" (id: "<<par.TrackId()<<") with mother trackID: "<<par.Mother()<<". Status Code: "<<par.StatusCode()<<" and proton energy "<<par.E()<<std::endl;


                            //if its mother is a delta with statuscode 3, and it has status code 14, then its the internal product of the delta decay.
                            if(par.StatusCode()==14 ){
                                
                                const  simb::MCParticle mother = truth->GetParticle(par.Mother());
                                if(is_delta_map.count(mother.PdgCode())>0 && mother.StatusCode()==3){
                                    m_mctruth_delta_proton_energy = par.E();
                                    tmp_n_protons_from_delta ++;
                                }
                            }


                            break;
                        }
                    case(2112): // if it's a neutron
                        {

                            m_mctruth_num_exiting_neutrons++;  // Guanqun: neutron always exits the nucleus? should check it
                                m_mctruth_exiting_neutron_mother_trackID.push_back(par.Mother());
                                m_mctruth_exiting_neutron_trackID.push_back(par.TrackId());
                                m_mctruth_exiting_neutron_energy.push_back(par.E());
                                m_mctruth_exiting_neutron_px.push_back(par.Px());
                                m_mctruth_exiting_neutron_py.push_back(par.Py());
                                m_mctruth_exiting_neutron_pz.push_back(par.Pz());

                      if(m_is_verbose)      std::cout<<"SingleProton::AnalyzeMCTruths()\t||\t Neutron "<<par.PdgCode()<<" (id: "<<par.TrackId()<<") with mother trackID: "<<par.Mother()<<". Status Code: "<<par.StatusCode()<<" and neutron energy "<<par.E()<<std::endl;

                            //if its mother is a delta with statuscode 3, and it has status code 14, then its the internal product of the delta decay.
                            if(par.StatusCode()==14){
                                const  simb::MCParticle mother = truth->GetParticle(par.Mother());
                                if(is_delta_map.count(mother.PdgCode())>0 && mother.StatusCode()==3){
                                m_mctruth_delta_neutron_energy = par.E();
                                tmp_n_neutrons_from_delta ++;
                                }
                            }
                        }

                        break;
                    case(-2224):
                    case(2224):
                        if(par.StatusCode() == 1){  m_mctruth_num_exiting_deltapp++; }
                        break;
                    case(-2214):
                    case(2214)://delta +
                    case(-1114):
                    case(1114): // if it's delta-
                        if(par.StatusCode() == 1){ m_mctruth_num_exiting_deltapm++; }
                        break;
                    case(-2114):
                    case(2114): // if it's delta0
                        if(par.StatusCode() == 1){ 
                            m_mctruth_num_exiting_delta0++;
                            m_mctruth_exiting_delta0_num_daughters.push_back(par.NumberDaughters());
                         if(m_is_verbose)   std::cout<<"SinglePhoton::AnalyzeMCTruths()\t||\t Delta0 "<<par.PdgCode()<<" (id: "<<par.TrackId()<<") with "<<m_mctruth_exiting_delta0_num_daughters.back()<<" daughters. StatusCode "<<par.StatusCode()<<std::endl;
                        }
                        break;
                    default:
                        break;
                }
            } // end of m_mctruth_num_daughter_particles loop

            if(m_is_textgen) continue; //quick hack, fix in files

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
		    // m_exiting_photon_energy_threshold is read from pset
		    if ( m_mctruth_exiting_photon_energy[p] > m_exiting_photon_energy_threshold){
				m_mctruth_num_reconstructable_protons++;

		    }//if g above threshold
                }

            //if it's a true delta radiative event, check the energies



            if (m_mctruth_is_delta_radiative==true){
                for (unsigned int p = 0; p <  m_mctruth_exiting_photon_energy.size(); p++){
                    std::cout<<"SinglePhoton::AnalyzeMCTruths()\t||\tLooking at exiting photon with energy "<<m_mctruth_exiting_photon_energy[p]<<std::endl;
                    if ( m_mctruth_exiting_photon_energy[p] > m_exiting_photon_energy_threshold){
                        m_mctruth_is_reconstructable_1g0p = true; // Guanqun: just means now we have a reconstructable shower, but we don't know if there is a reconstructed nucleon or not yet.

                    }//if g above threshold
                }//for all exiting g
                for(unsigned int pr = 0; pr <  m_mctruth_exiting_proton_energy.size(); pr++){
                   if ( m_mctruth_exiting_proton_energy[pr]> m_exiting_proton_energy_threshold){
                        //if it's already 1g1p then we've found a 1g2p which we aren't counting
                        // Guanqun: limit to only 1 reconstructable proton?
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

		// Guanqun: why not consider its first-generation mother?
		// for a photon exiting nucleus, its first mother is always also a photon (photon exits the nucleus, it becomes another photon..)
                while(nth_mother.StatusCode() != 0 || n_generation < 4){

                    if(nth_mother.Mother()<0) break;
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

                while(nth_mother.StatusCode() != 0 && n_generation < 4){
                    if(nth_mother.Mother()<0) break;
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

	    // if it's a pi0, its mother trackID is 0 and it has two daughters
            if(mcp->PdgCode() == 111 && mcp->Mother() == 0 && mcp->NumberDaughters()==2 ){
                npi0check++;
		// get its two daughters
                const art::Ptr<simb::MCParticle> dau1 = mcParticleVector[mymap[mcp->Daughter(0)]];
                const art::Ptr<simb::MCParticle> dau2 = mcParticleVector[mymap[mcp->Daughter(1)]];

                if(false)  std::cout<<"On Dau1: "<<" Mother:"<< dau1->Mother()<<" pdgcode: "<<dau1->PdgCode()<<" trkid: "<<dau1->TrackId()<<" statuscode: "<<dau1->StatusCode()<<std::endl;
                if(false)  std::cout<<"On Dau2: "<<" Mother:"<< dau2->Mother()<<" pdgcode: "<<dau2->PdgCode()<<" trkid: "<<dau2->TrackId()<<" statuscode: "<<dau2->StatusCode()<<std::endl;

                double e1 = dau1->E();
                double e2 = dau2->E();

                std::vector<double> raw_1_End  ={dau1->EndX(), dau1->EndY(), dau1->EndZ()};
                std::vector<double> raw_1_Start  ={dau1->Vx(), dau1->Vy(), dau1->Vz()};
            
                std::vector<double> raw_2_End  ={dau2->EndX(), dau2->EndY(), dau2->EndZ()};
                std::vector<double> raw_2_Start  ={dau2->Vx(), dau2->Vy(), dau2->Vz()};

                std::vector<double> corrected_1_start(3), corrected_2_start(3);
                std::vector<double> corrected_1_end(3), corrected_2_end(3);

                this->spacecharge_correction(dau1, corrected_1_start, raw_1_Start);
                this->spacecharge_correction(dau1, corrected_1_end, raw_1_End);
                
                this->spacecharge_correction(dau2, corrected_2_start, raw_2_Start);
                this->spacecharge_correction(dau2, corrected_2_end, raw_2_End);

                for(int p1=0; p1<dau1->NumberDaughters();p1++){
                        auto dd = mcParticleVector[mymap[dau1->Daughter(p1)]];
                        std::cout<<"Post1 "<<dd->PdgCode()<<" "<<dd->TrackId()<<" "<<dd->StatusCode()<<" "<<dd->EndProcess()<<std::endl;
                }
 
                for(int p1=0; p1<dau2->NumberDaughters();p1++){
                        auto dd = mcParticleVector[mymap[dau2->Daughter(p1)]];
                        std::cout<<"Post2 "<<dd->PdgCode()<<" "<<dd->TrackId()<<" "<<dd->StatusCode()<<" "<<dd->EndProcess()<<" "<<dd->E()<<std::endl;
                }
                
                int exit1 = this->isInTPCActive(corrected_1_end);
                int exit2 = this->isInTPCActive(corrected_2_end);

                if(e2<e1){
                    m_mctruth_pi0_leading_photon_energy = e1;
                    m_mctruth_pi0_subleading_photon_energy = e2;
                    m_mctruth_pi0_leading_photon_end_process = dau1->EndProcess();
                    m_mctruth_pi0_subleading_photon_end_process = dau2->EndProcess();
                    m_mctruth_pi0_leading_photon_start = corrected_1_start;
                    m_mctruth_pi0_leading_photon_end = corrected_1_end;
                    m_mctruth_pi0_subleading_photon_start = corrected_2_start;
                    m_mctruth_pi0_subleading_photon_end = corrected_2_end;
         		    //note: the order of subleading/leading photon is reversed// Fixed as of 2022 reprocess!
                    m_mctruth_pi0_leading_photon_exiting_TPC =exit1; 
                    m_mctruth_pi0_subleading_photon_exiting_TPC = exit2;
                    m_mctruth_pi0_leading_photon_mom = {dau1->Px(),dau1->Py(),dau1->Pz()};
                    m_mctruth_pi0_subleading_photon_mom = {dau2->Px(),dau2->Py(),dau2->Pz()};

                }else{
                    m_mctruth_pi0_leading_photon_energy = e2;
                    m_mctruth_pi0_subleading_photon_energy = e1;
                    m_mctruth_pi0_leading_photon_end_process = dau2->EndProcess();
                    m_mctruth_pi0_subleading_photon_end_process = dau1->EndProcess();
                    m_mctruth_pi0_leading_photon_start = corrected_2_start;
                    m_mctruth_pi0_leading_photon_end = corrected_2_end;
                    m_mctruth_pi0_subleading_photon_start = corrected_1_start;
                    m_mctruth_pi0_subleading_photon_end = corrected_1_end;
                    m_mctruth_pi0_leading_photon_exiting_TPC = exit2;
                    m_mctruth_pi0_subleading_photon_exiting_TPC = exit1;
                    m_mctruth_pi0_subleading_photon_mom = {dau1->Px(),dau1->Py(),dau1->Pz()};
                    m_mctruth_pi0_leading_photon_mom = {dau2->Px(),dau2->Py(),dau2->Pz()};

                }

            }
        }

        if(npi0check>1) std::cout<<"WARNING WARNING!!!! there are "<<npi0check<<" Pi0's in this event in geant4 that come from the nucleas"<<std::endl;

    }//end of analyze this


}
