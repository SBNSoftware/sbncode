namespace single_photon
{

	//set the vertex for now;
	void SinglePhoton::Output_PFParticleInfo( std::vector<PandoraPFParticle> PPFPs){

		int pfp_size = PPFPs.size();

		for(int index = 0; index < pfp_size; index++){

			PandoraPFParticle* temp_p = &PPFPs[index];
			if(!(pfp_w_bestnuID == temp_p->get_SliceID() && temp_p->get_IsNeutrino()) ) continue;
			m_vertex_pos_x = temp_p->get_Vertex_pos()[0];
			m_vertex_pos_y = temp_p->get_Vertex_pos()[1];
			m_vertex_pos_z = temp_p->get_Vertex_pos()[2];
			std::cout<<"Best NuScore is found, define the vertice as: ("<<temp_p->get_Vertex_pos()[0]<<","<<temp_p->get_Vertex_pos()[1]<<","<<temp_p->get_Vertex_pos()[2]<<")"<<std::endl;
		
			std::vector<double> tmp = {m_vertex_pos_x, m_vertex_pos_y, m_vertex_pos_z};
			m_reco_vertex_in_SCB = this->distToSCB(m_reco_vertex_dist_to_SCB,tmp);
			m_reco_vertex_dist_to_active_TPC =  this->distToTPCActive(tmp);
			m_reco_vertex_dist_to_CPA =  this->distToCPA(tmp);
		}
	}

}
