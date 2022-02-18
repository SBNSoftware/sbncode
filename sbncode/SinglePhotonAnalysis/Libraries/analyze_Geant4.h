#include "SinglePhoton_module.h"

namespace single_photon
{
    void SinglePhoton::ClearGeant4Branches(){

        m_geant4_pdg.clear();
        m_geant4_trackid.clear();
        m_geant4_mother.clear();
        m_geant4_statuscode.clear();
        m_geant4_E.clear();
        m_geant4_mass.clear();
        m_geant4_px.clear();
        m_geant4_py.clear();
        m_geant4_pz.clear();
        m_geant4_dx.clear();
        m_geant4_dy.clear();
        m_geant4_dz.clear();

        m_geant4_vx.clear();
        m_geant4_vy.clear();
        m_geant4_vz.clear();
        m_geant4_process.clear();
        m_geant4_end_process.clear();


        m_geant4_costheta.clear();


    }

    void SinglePhoton::CreateGeant4Branches(){
        geant4_tree->Branch("geant4_pdg",&m_geant4_pdg);
        geant4_tree->Branch("geant4_trackid",&m_geant4_trackid);
        geant4_tree->Branch("geant4_mother",&m_geant4_mother);
        geant4_tree->Branch("geant4_statuscode",&m_geant4_statuscode);
        geant4_tree->Branch("geant4_E",&m_geant4_E);
        geant4_tree->Branch("geant4_mass",&m_geant4_mass);
        geant4_tree->Branch("geant4_px", &m_geant4_px);
        geant4_tree->Branch("geant4_py", &m_geant4_py);
        geant4_tree->Branch("geant4_pz", &m_geant4_pz);

        geant4_tree->Branch("geant4_dx", &m_geant4_dx);
        geant4_tree->Branch("geant4_dy", &m_geant4_dy);
        geant4_tree->Branch("geant4_dz", &m_geant4_dz);

        geant4_tree->Branch("geant4_vx", &m_geant4_vx);
        geant4_tree->Branch("geant4_vy", &m_geant4_vy);
        geant4_tree->Branch("geant4_vz", &m_geant4_vz);
        geant4_tree->Branch("geant4_costheta",&m_geant4_costheta);

        geant4_tree->Branch("geant4_end_process", &m_geant4_end_process);
        geant4_tree->Branch("geant4_process", &m_geant4_process);


    }

    void SinglePhoton::AnalyzeGeant4( const    std::vector<art::Ptr<simb::MCParticle>> &mcParticleVector){    

        if(m_is_verbose) std::cout<<"SinglePhoton::AnalyzeGeant4s()\t||\t Begininning recob::Geant4 analysis suite"<<std::endl;;


        for(size_t j=0;j< mcParticleVector.size();j++){

            const art::Ptr<simb::MCParticle> mcp = mcParticleVector[j];
            std::cout<<"PARG: "<<j<<" PDG "<<mcp->PdgCode()<<" Status "<<mcp->StatusCode()<<" trackid: "<<mcp->TrackId()<<" Mothe "<<mcp->Mother()<<" Process "<<mcp->Process()<<" EndProcess "<<mcp->EndProcess()<<" Energy "<<mcp->E()<<" start ("<<mcp->Vx()<<","<<mcp->Vy()<<","<<mcp->Vz()<<")"<<std::endl;
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
