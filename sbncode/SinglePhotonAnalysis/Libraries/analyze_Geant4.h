namespace single_photon
{


    void SinglePhoton::AnalyzeGeant4( const    std::vector<art::Ptr<simb::MCParticle>> &mcParticleVector){    


		std::vector<int> spacers = Printer_header({"#MCP","   pdg", " Status"," trkID"," Mother"," Process", "            Process_End","     Energy", "      Vertex(x,  ","       y,      ",",      z     )"});
        for(size_t j=0;j< mcParticleVector.size();j++){

            const art::Ptr<simb::MCParticle> mcp = mcParticleVector[j];
//            std::cout<<"PARG: "<<j<<" PDG "<<mcp->PdgCode()<<" Status "<<mcp->StatusCode()<<" trackid: "<<mcp->TrackId()<<" Mothe "<<mcp->Mother()<<" Process "<<mcp->Process()<<" EndProcess "<<mcp->EndProcess()<<" Energy "<<mcp->E()<<" start ("<<mcp->Vx()<<","<<mcp->Vy()<<","<<mcp->Vz()<<")"<<std::endl;

			Printer_content({
					std::to_string(j),
					std::to_string(mcp->PdgCode()),
					std::to_string(mcp->StatusCode()),
					std::to_string(mcp->TrackId()),
					std::to_string(mcp->Mother()),
					mcp->Process(),
					mcp->EndProcess(),
					std::to_string(mcp->E()),
					std::to_string(mcp->Vx()),
					std::to_string(mcp->Vy()),
					std::to_string(mcp->Vz())
					},
					spacers);

            m_geant4_pdg.push_back(mcp->PdgCode());
            m_geant4_trackid.push_back(mcp->TrackId());
            m_geant4_statuscode.push_back(mcp->StatusCode());
            m_geant4_mother.push_back(mcp->Mother());
            m_geant4_E.push_back(mcp->E());
            m_geant4_mass.push_back(mcp->Mass());
            m_geant4_px.push_back(mcp->Px());
            m_geant4_py.push_back(mcp->Py());
            m_geant4_pz.push_back(mcp->Pz());
            m_geant4_vx.push_back(mcp->Vx());
            m_geant4_vy.push_back(mcp->Vy());
            m_geant4_vz.push_back(mcp->Vz());
            m_geant4_end_process.push_back(mcp->EndProcess());
            m_geant4_process.push_back(mcp->Process());
            m_geant4_costheta.push_back(m_geant4_pz.back()/sqrt(pow(m_geant4_pz.back(),2)+pow(m_geant4_px.back(),2)+pow(m_geant4_py.back(),2)));
            m_geant4_dx.push_back(mcp->Px()/sqrt(pow(m_geant4_pz.back(),2)+pow(m_geant4_px.back(),2)+pow(m_geant4_py.back(),2)));
            m_geant4_dy.push_back(mcp->Py()/sqrt(pow(m_geant4_pz.back(),2)+pow(m_geant4_px.back(),2)+pow(m_geant4_py.back(),2)));
            m_geant4_dz.push_back(mcp->Pz()/sqrt(pow(m_geant4_pz.back(),2)+pow(m_geant4_px.back(),2)+pow(m_geant4_py.back(),2)));



            if(j>2)break;

        }


    }
}
