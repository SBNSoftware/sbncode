#include "SinglePhoton_module.h"

namespace single_photon
{
    void SinglePhoton::ClearTemplates(){
        m_reco_num_templates = 0;
        m_reco_template.clear();
    }

    void SinglePhoton::ResizeTemplates(size_t size){
        m_reco_template.resize(size);
    }


    void SinglePhoton::CreateTemplateBranches(){
        vertex_tree->Branch("reco_template",&m_reco_template);
    }

    void SinglePhoton::AnalyzeTemplates(){
        m_reco_num_templates = 1;
        this->ResizeTemplates(m_reco_num_templates);

        m_reco_template[0]=99;

    }
}
