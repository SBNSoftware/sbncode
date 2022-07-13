namespace single_photon
{
    void SinglePhoton::AnalyzeFlashes(const std::vector<art::Ptr<recob::OpFlash>>& flashes, art::Handle<std::vector<sbn::crt::CRTHit>> crthit_h, double evt_timeGPS_nsec,  std::map<art::Ptr<recob::OpFlash>, std::vector< art::Ptr<sbn::crt::CRTHit>>> crtvetoToFlashMap){

        //  void SinglePhoton::AnalyzeFlashes(const std::vector<art::Ptr<recob::OpFlash>>& flashes, art::Handle<std::vector<sbn::crt::CRTHit>> crthit_h){}

      
        for(auto pair: crtvetoToFlashMap){
            std::cout<<"for flash at time "<< pair.first->Time()<<" has "<< pair.second.size() << " associated  CRT hits "<<std::endl;
            if(pair.second.size() > 0){
                for (auto hit: pair.second){
                    std::cout<<"---- associated CRT hit at time "<<hit->ts0_ns/1000. <<" with PE "<<hit->peshit<<std::endl;
                    m_CRT_veto_hit_PE.push_back(hit->peshit);
                }

            }
            m_CRT_veto_nhits =  pair.second.size();//save the number of associated CRT veto hits
        }

        size_t flash_size = flashes.size();
        m_reco_num_flashes = flash_size;

        if(m_is_verbose) std::cout<<"SinglePhoton::AnalyzeFlashes()\t||\t Beginning analysis of recob::OpFlash\n";

        this->ResizeFlashes(flash_size);

        for(size_t i = 0; i < flash_size; ++i) {


            art::Ptr<recob::OpFlash> const & flash = flashes[i];

            m_reco_flash_total_pe[i]=(flash->TotalPE());
            m_reco_flash_time[i]=(flash->Time());
            m_reco_flash_time_width[i]=flash->TimeWidth();
            m_reco_flash_abs_time[i]=flash->AbsTime();
            m_reco_flash_frame[i]=flash->Frame();
            m_reco_flash_ycenter[i]=flash->YCenter();
            m_reco_flash_ywidth[i]=flash->YWidth();
            m_reco_flash_zcenter[i]=flash->ZCenter();
            m_reco_flash_zwidth[i]=flash->ZWidth();

      // m_beamgate_flash_end/m_beamgate_flash_start are read from pset
            if(m_reco_flash_time[i] <= m_beamgate_flash_end && m_reco_flash_time[i] >= m_beamgate_flash_start){
                m_reco_num_flashes_in_beamgate++;
                m_reco_flash_total_pe_in_beamgate[i]=(flash->TotalPE());
                m_reco_flash_time_in_beamgate[i]=(flash->Time());
                m_reco_flash_ycenter_in_beamgate[i] = flash->YCenter();
                m_reco_flash_zcenter_in_beamgate[i] = flash->ZCenter();
            }



        }

        if(m_is_verbose) std::cout<<"SinglePhoton::AnalyzeFlashes()\t||\t Finished. There was "<<flash_size<<" flashes with: "<<m_reco_num_flashes_in_beamgate<<" in the beamgate defined by: "<<m_beamgate_flash_start<<" <-> "<<m_beamgate_flash_end<<std::endl;

        //fill these values only for events that have CRT information - run3 G and later
        //code taken from ubcrt/UBCRTCosmicFilter/UBCRTCosmicFilter_module.cc
        if(m_runCRT){
            if (m_reco_num_flashes_in_beamgate == 1){ //fill only if there's a flash in the beamgate

                int  _nCRThits_in_event = crthit_h->size();

                double _dt_abs   = 100000.0;
                //  double  _within_resolution = 0;
                double _beam_flash_time  =  m_reco_flash_time_in_beamgate[0];  // Guanqun: why use index 0?

                // Loop over the CRT hits.
                for (int j = 0; j < _nCRThits_in_event; j++)
                {
                    /*
                       if (verbose)
                       std::cout << "\t Time of the CRT Hit wrt the event timestamp = " << ((crthit_h->at(j).ts0_ns - evt_timeGPS_nsec + fDTOffset) / 1000.) << " us." << std::endl;
                       */
                    double _crt_time_temp = ((crthit_h->at(j).ts0_ns - evt_timeGPS_nsec + m_DTOffset) / 1000.);

                    // Fill the vector variables.
                    m_CRT_hits_time.push_back(_crt_time_temp);
                    m_CRT_hits_PE.push_back(crthit_h->at(j).peshit);
                    m_CRT_hits_x.push_back(crthit_h->at(j).x_pos);
                    m_CRT_hits_y.push_back(crthit_h->at(j).y_pos);
                    m_CRT_hits_z.push_back(crthit_h->at(j).z_pos);

                    if (fabs(_beam_flash_time - _crt_time_temp) < _dt_abs)
                    {
                        _dt_abs = fabs(_beam_flash_time - _crt_time_temp);
                        m_CRT_dt = _beam_flash_time - _crt_time_temp;
                        m_CRT_min_hit_time = _crt_time_temp;
                        // set 'within_resolution' to 'true' and break the loop if 'closest_crt_diff' is less than fResolution.
                        if (_dt_abs < m_Resolution)
                        {
                            //_within_resolution = 1;
                            // Set the position information and the intensity of the CRT hit.
                            m_CRT_min_hit_PE = crthit_h->at(j).peshit;
                            m_CRT_min_hit_x = crthit_h->at(j).x_pos;
                            m_CRT_min_hit_y = crthit_h->at(j).y_pos;
                            m_CRT_min_hit_z = crthit_h->at(j).z_pos;


                            // if (verbose)
                            // {
                            std::cout << "CRT hit PE = " << m_CRT_min_hit_PE << " PEs." << std::endl;
                            std::cout << "CRT hit x = " << m_CRT_min_hit_x << " cm." << std::endl;
                            std::cout << "CRT hit y = " << m_CRT_min_hit_y << " cm." << std::endl;
                            std::cout << "CRT hit z = " << m_CRT_min_hit_z << " cm." << std::endl;
                            // }
                            break;
                        }
                    } // End of conditional for closest CRT hit time.
                } // End of loop over CRT hits.
            } //if there is 1 flash in beamgate
        }//if runCRT


    }//analyze flashes


    }
