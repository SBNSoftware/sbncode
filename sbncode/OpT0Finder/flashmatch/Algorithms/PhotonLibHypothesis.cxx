#ifndef PHOTONLIBHYPOTHESIS_CXX
#define PHOTONLIBHYPOTHESIS_CXX

#include "PhotonLibHypothesis.h"
#include <assert.h>

#ifndef USING_LARSOFT
#define USING_LARSOFT 1
#endif

#if USING_LARSOFT == 0
#include <omp.h>
#define NUM_THREADS 12
#endif

//using namespace std::chrono;

namespace flashmatch {

    static PhotonLibHypothesisFactory __global_PhotonLibHypothesisFactory__;

    PhotonLibHypothesis::PhotonLibHypothesis(const std::string name)
    : BaseFlashHypothesis(name)
    {
        #if USING_LARSOFT == 0
        omp_set_num_threads(NUM_THREADS); 
        #endif
    }

    void PhotonLibHypothesis::_Configure_(const Config_t &pset)
    {

        _global_qe = pset.get<double>("GlobalQE");
        _global_qe_refl = pset.get<double>("GlobalQERefl", -1);
        _reco_pe_calib = pset.get<double>("RecoPECalibFactor",1.0);
        _extend_tracks = pset.get<bool>("ExtendTracks", false);
        /*
        if(_extend_tracks) {
            FLASH_CRITICAL() << "ExtendTracks currently has a problem and cannot be enabled!" << std::endl
            << "If you would like to debug, here is useful info..." << std::endl
            << "1. InspectTouchingEdges function checks if a track is touching either of 2 PMT planes (within 1 cryostat)" << std::endl
            << "2. But TrackExtension function only extends to the negative x direction (and it does so disregard of the fact if that end is touching or not)" << std::endl
            << "This is completely wrong (a bug). Hence we disabled at the moment. TrackExtension may not be critical for flash matching performance" <<std::endl;
            throw OpT0FinderException();
        } 
        */
        _threshold_proximity = pset.get<double>("ExtensionProximityThreshold", 5.0);
        _threshold_track_len = pset.get<double>("ExtensionTrackLengthMaxThreshold", 20.0);
        _segment_size = pset.get<double>("SegmentSize", 0.5);

        _qe_v.clear();
        _qe_v = pset.get<std::vector<double> >("CCVCorrection",_qe_v);
        if(_qe_v.empty()) _qe_v.resize(DetectorSpecs::GetME().NOpDets(),1.0);
        if(_qe_v.size() != DetectorSpecs::GetME().NOpDets()) {
          FLASH_CRITICAL() << "CCVCorrection factor array has size " << _qe_v.size()
          << " != number of opdet (" << DetectorSpecs::GetME().NOpDets() << ")!" << std::endl;
          throw OpT0FinderException();
        }

    }

    int PhotonLibHypothesis::InspectTouchingEdges(const QCluster_t& trk) const
    {
        // Return code: 0=no touch, 1=start touching, 2=end touching, 3= both touching
        int result = 0;

        // Check if the start/end is near the edge
        auto const& start = trk.front();
        auto const& end   = trk.back();
        if (((start.x - DetectorSpecs::GetME().ActiveVolume().Min()[0]) < _threshold_proximity) ||
            ((start.y - DetectorSpecs::GetME().ActiveVolume().Min()[1]) < _threshold_proximity) ||
            ((start.z - DetectorSpecs::GetME().ActiveVolume().Min()[2]) < _threshold_proximity) ||
            ((DetectorSpecs::GetME().ActiveVolume().Max()[0] - start.x) < _threshold_proximity) ||
            ((DetectorSpecs::GetME().ActiveVolume().Max()[1] - start.y) < _threshold_proximity) ||
            ((DetectorSpecs::GetME().ActiveVolume().Max()[2] - start.z) < _threshold_proximity)) {
                //std::cout << "*** Extending track " << trk.size() << " " << min_x << " " << max_x << std::endl;
                //std::cout << trk.front().x << " " << trk.back().x << std::endl;
            result += 1;
        }
        if (((end.x - DetectorSpecs::GetME().ActiveVolume().Min()[0]) < _threshold_proximity) ||
            ((end.y - DetectorSpecs::GetME().ActiveVolume().Min()[1]) < _threshold_proximity) ||
            ((end.z - DetectorSpecs::GetME().ActiveVolume().Min()[2]) < _threshold_proximity) ||
            ((DetectorSpecs::GetME().ActiveVolume().Max()[0] - end.x) < _threshold_proximity) ||
            ((DetectorSpecs::GetME().ActiveVolume().Max()[1] - end.y) < _threshold_proximity) ||
            ((DetectorSpecs::GetME().ActiveVolume().Max()[2] - end.z) < _threshold_proximity)) {
                //std::cout << "*** Extending track " << trk.size() << " " << min_x << " " << max_x << std::endl;
                //std::cout << trk.front().x << " " << trk.back().x << std::endl;
            result += 2;
        }
        /*
        // OLD CODE looking at ANY POINT close to the edge
        double min_x = kINVALID_DOUBLE; double max_x = -kINVALID_DOUBLE;
        min_idx = max_idx = 0;
        for (size_t pt_index = 0; pt_index < trk.size(); ++pt_index) {
            if (trk[pt_index].x < min_x) {
                min_x = trk[pt_index].x;
                min_idx = pt_index;
            }
            if (trk[pt_index].x > max_x) {
                max_x = trk[pt_index].x;
                max_idx = pt_index;
            }
        }
        if (((trk[min_idx].x - DetectorSpecs::GetME().ActiveVolume().Min()[0]) < _threshold_proximity) ||
            ((trk[min_idx].y - DetectorSpecs::GetME().ActiveVolume().Min()[1]) < _threshold_proximity) ||
            ((trk[min_idx].z - DetectorSpecs::GetME().ActiveVolume().Min()[2]) < _threshold_proximity) ||
            ((DetectorSpecs::GetME().ActiveVolume().Max()[0] - trk[max_idx].x) < _threshold_proximity) ||
            ((DetectorSpecs::GetME().ActiveVolume().Max()[1] - trk[max_idx].y) < _threshold_proximity) ||
            ((DetectorSpecs::GetME().ActiveVolume().Max()[2] - trk[max_idx].z) < _threshold_proximity)) {
                //std::cout << "*** Extending track " << trk.size() << " " << min_x << " " << max_x << std::endl;
                //std::cout << trk.front().x << " " << trk.back().x << std::endl;
            return true;
        }
        */
        return result;

    }

    QCluster_t PhotonLibHypothesis::ComputeExtension(const geoalgo::Vector& A, const geoalgo::Vector& B) const 
    {
        QCluster_t extension;

        // direction to extend should be A=>B ... so B should be closer to the edge
        geoalgo::Vector pt = B;
        geoalgo::Vector AB = B - A;
        auto length = AB.Length();
        if(length < 1.e-9) {
            FLASH_WARNING() << "Extension calculation halted as the direction estimate is unreliable for a short segment (" 
            << length << " cm)" << std::endl;
            return extension;
        }
        AB.Normalize();
        const auto& box0 = DetectorSpecs::GetME().PhotonLibraryVolume();
        const auto& box1 = DetectorSpecs::GetME().ActiveVolume();

        double dist_inside = 0;
        extension.reserve(100);
        pt += (AB * _segment_size/2.);
        while(box0.Contain(pt) && dist_inside < _threshold_proximity) {
          if(!box1.Contain(pt)) {
            QPoint_t qpt;
            qpt.x = pt[0];
            qpt.y = pt[1];
            qpt.z = pt[2];
            //std::cout << qpt.x << " " << qpt.y << " " << qpt.z << std::endl;
            qpt.q = _segment_size * DetectorSpecs::GetME().LightYield() * DetectorSpecs::GetME().MIPdEdx();
            extension.emplace_back(qpt);
          }
          else{ dist_inside += _segment_size; }

          pt += (AB * _segment_size);
        }

        return extension;
    }

    QCluster_t PhotonLibHypothesis::TrackExtension(const QCluster_t& in_trk, const int touch) const
    {
        auto start_touch = bool(touch & 0x1);
        auto end_touch   = bool(touch & 0x2);

        QCluster_t trk;
        if(start_touch) {
            geoalgo::Vector B(in_trk[0].x, in_trk[0].y, in_trk[0].z);
            geoalgo::Vector A(3);
            // search for a canddiate point to define a direction (and avoid using the same point)
            for(size_t i=1; i<in_trk.size(); ++i) {
                A[0] = in_trk[i].x;
                A[1] = in_trk[i].y;
                A[2] = in_trk[i].z;
                if(A.Dist(B) > _segment_size) break;
            }
            trk = this->ComputeExtension(A,B);
        }
        trk.reserve(trk.size() + in_trk.size());
        for(auto const& pt : in_trk) trk.push_back(pt);
        if(end_touch) {
            geoalgo::Vector B(in_trk[in_trk.size()-1].x, in_trk[in_trk.size()-1].y, in_trk[in_trk.size()-1].z);
            geoalgo::Vector A = B;
            // search for a canddiate point to define a direction (and avoid using the same point)
            for(int i=in_trk.size()-2; i>=0; --i) {
                A[0] = in_trk[i].x;
                A[1] = in_trk[i].y;
                A[2] = in_trk[i].z;
                if(A.Dist(B) > _segment_size) break;
            }
            auto end_trk = this->ComputeExtension(A,B);
            trk.reserve(trk.size()+end_trk.size());
            for(auto const& pt : end_trk)
                trk.push_back(pt);
        }

        /*

        QCluster_t trk = in_trk;
        // Extend the track
        // Compute coordinates of final point first
        geoalgo::Vector A(trk[max_idx].x, trk[max_idx].y, trk[max_idx].z);
        geoalgo::Vector B(trk[min_idx].x, trk[min_idx].y, trk[min_idx].z);
        geoalgo::Vector AB = B - A;

        double x_C = DetectorSpecs::GetME().PhotonLibraryVolume().Min()[0];
        double lengthBC = (x_C - B[0])/(B[0] - A[0]) * AB.Length();
        geoalgo::Vector C = B + AB / AB.Length() * lengthBC;

                // Add to _var_trk the part betwen boundary and C
                //_custom_algo->MakeQCluster(C, B, _var_trk, -1);
        geoalgo::Vector unit = (C - B).Dir();

        QPoint_t q_pt;
        geoalgo::Vector current = B;
        int num_pts = int(lengthBC / _segment_size);
        trk.reserve(trk.size() + num_pts);
        for (int i = 0; i < num_pts+1; i++) {
            double current_segment_size = (i < num_pts ? _segment_size : (lengthBC - _segment_size*num_pts));
            current = current + unit * current_segment_size/2.0;
            q_pt.x = current[0];
            q_pt.y = current[1];
            q_pt.z = current[2];
            q_pt.q = current_segment_size * DetectorSpecs::GetME().LightYield() * DetectorSpecs::GetME().MIPdEdx();
            if (trk.front().x < trk.back().x) {
                trk.insert(trk.begin(), q_pt);
            }
            else {
                trk.emplace_back(q_pt);
            }
            current = current + unit * current_segment_size/2.0;
            //std::cout << "Adding point " << current  << " " << i << " " << num_pts << std::endl;
        }
        //std::cout << " done " << trk.size() << std::endl;
        //std::cout << trk.front().x << " " << trk.back().x << std::endl;
        */
        return trk;
    }


    void PhotonLibHypothesis::FillEstimate(const QCluster_t& tpc_trk, Flash_t &flash) const
    {
        size_t n_pmt = DetectorSpecs::GetME().NOpDets();//n_pmt returns 0 now, needs to be fixed
        if(flash.pe_v.empty()) flash.pe_v.resize(n_pmt);
        if(flash.pe_err_v.empty()) flash.pe_err_v.resize(n_pmt);
        if(flash.pe_true_v.empty()) flash.pe_true_v.resize(n_pmt);

        assert(flash.pe_v.size()     == n_pmt);
        assert(flash.pe_true_v.size() == n_pmt);
        assert(flash.pe_err_v.size() == n_pmt);

        for (auto& v : flash.pe_v      ) {v = 0;}
        for (auto& v : flash.pe_err_v  ) {v = 0;}
        for (auto& v : flash.pe_true_v ) {v = 0;}

        double track_length = tpc_trk.front().dist(tpc_trk.back());
        int touch = this->InspectTouchingEdges(tpc_trk);
        bool extend_tracks = (_extend_tracks && touch && track_length < _threshold_track_len);
        if(_extend_tracks) {
            FLASH_DEBUG() << "Extend? " << extend_tracks << " ... track length " << track_length << " touch-or-not " << touch << std::endl;
        }


        if(!extend_tracks) 
            this->BuildHypothesis(tpc_trk,flash);
        else {
            auto trk = this->TrackExtension(tpc_trk,touch);
            this->BuildHypothesis(trk,flash);
        }

    }

    void PhotonLibHypothesis::BuildHypothesis(const QCluster_t& trk, Flash_t &flash) const
    {
        //bool speak=false;
        //std::cout<<trk[0].x<<" "<<trk[0].y<<" "<<trk[0].z<<std::endl;
        //if(std::fabs(trk[0].x + 211.635)<1.0) speak=true;
        size_t n_pmt = DetectorSpecs::GetME().NOpDets();
	#if USING_LARSOFT == 0
        #pragma omp parallel
        #endif
        {
            size_t num_pts  = trk.size();
            size_t start_pt = 0;

            #if USING_LARSOFT == 0
            size_t thread_id = omp_get_thread_num();
            size_t num_threads = omp_get_num_threads();
            if(num_threads == 1 || num_threads>num_pts) {
                if(thread_id > 0) {
                    start_pt = 0;
                }
            }else{
                num_pts = trk.size() / (num_threads-1);
                start_pt = num_pts * thread_id;
                if(thread_id+1 == num_threads) num_pts = (trk.size() % (num_threads-1));
                //sleep(thread_id);
                //std::cout<<"Start " <<start_pt << " num pts " << num_pts << " / " << trk.size() << std::endl;
            }
            #endif

            auto const& vox_def = DetectorSpecs::GetME().GetVoxelDef();
            std::vector<double> local_pe_v(n_pmt,0.);
            std::vector<double> local_pe_refl_v(n_pmt,0.);
            double pos[3];
            size_t nproc0,nproc1; // for debug cout
            nproc0 = nproc1 = 0;  // for debug cout

            for(size_t ipt = start_pt; ipt < start_pt+num_pts; ++ipt) {

                auto const& pt = trk[ipt];
                pos[0] = pt.x;
                pos[1] = pt.y;
                pos[2] = pt.z;
                nproc0 += 1; // for debug cout

                int vox_id = vox_def.GetVoxelID(pos);
                if (vox_id < 0) continue;

		auto const& lib_data = DetectorSpecs::GetME().GetLibraryEntries(vox_id);

		#if USING_LARSOFT == 1
		auto const& lib_data_refl = DetectorSpecs::GetME().GetLibraryEntries(vox_id,true);
		#endif

                nproc1 += 1; // for debug cout
                double qsum = 0.;
                double vsum = 0.;
                for(size_t ipmt=0; ipmt < n_pmt; ++ipmt) {

                    if(_channel_mask[ipmt]) continue;

                    if(!_uncoated_pmt_list[ipmt])
		      local_pe_v[ipmt] += pt.q * lib_data[ipmt];
            qsum += pt.q * lib_data[ipmt];
            vsum += lib_data[ipmt];
		    #if USING_LARSOFT == 1
		    if(_global_qe_refl > 0.)
		      local_pe_refl_v[ipmt] += pt.q * lib_data_refl[ipmt];
                    #endif
                    //local_pe_v[ipmt] += pt.q * vis_pmt[ipmt];
                    //std::cout << "PMT : " << ipmt << " [x,y,z] -> [q] : [" << pt.x << ", " << pt.y << ", " << pt.z << "] -> [" << q0 << "," << q1 << "]" << std::endl;
                }
                //if(speak) std::cout << pt.x << " " << pt.y << " " << pt.z << " ... " << pt.q << ", " << qsum << ", " << vsum << std::endl;
            }

            #if USING_LARSOFT == 0
            #pragma omp critical
            #endif
            {
                double qsum = 0.; // for debug cout
                for(size_t ipmt = 0; ipmt < n_pmt; ++ipmt) {
                    double q0 = (local_pe_v[ipmt] * _global_qe * _reco_pe_calib / _qe_v[ipmt]);
                    double q1 = (local_pe_refl_v[ipmt] * _global_qe_refl * _reco_pe_calib / _qe_v[ipmt]);
                    double q = q0 + q1;
                    flash.pe_v[ipmt] +=  q;
                    qsum += q; // for debug cout
                }

                // for debug cout
                //double qsum2 = 0.;
                //for(auto const& v: flash.pe_v) qsum2 += v;
                //std::cout<<"Thread ID " << thread_id << " ... " << start_pt << " => " << start_pt+num_pts << " ... " << nproc0 << "/" << nproc1 << " qsum " << qsum << " total sum " << qsum2 << std::endl<<std::flush;
                //std::cout<<qsum<<std::endl;
            }

        }
        // for debug cout
        //double charge = 0.;
        //for(auto const& v : flash.pe_v) { charge+=v; }
        //if(charge>0)
        //    std::cout << std::endl << "Track size " << trk.size() << " points, total pe " << charge << std::endl << std::flush;
        //sleep(3);
        //throw OpT0FinderException();
        return;
    }

}
#endif
