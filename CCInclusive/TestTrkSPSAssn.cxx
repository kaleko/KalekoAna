#ifndef LARLITE_TESTTRKSPSASSN_CXX
#define LARLITE_TESTTRKSPSASSN_CXX

#include "TestTrkSPSAssn.h"
#include "DataFormat/track.h"
#include "DataFormat/spacepoint.h"
#include "DataFormat/mctrack.h"

namespace larlite {

    bool TestTrkSPSAssn::initialize() {

        _tmc = new kaleko::TrackMomentumCalculator();
        _tmc->SetStepSize(_seglen);

        _tree = new TTree("tree", "tree");
        _tree->Branch("n_total_sps", &_n_total_sps, "n_total_sps/I");
        _tree->Branch("n_asstd_sps", &_n_asstd_sps, "n_asstd_sps/I");
        _tree->Branch("n_track_traj_pts", &_n_track_traj_pts, "n_track_traj_pts/I");
        _tree->Branch("dist_firstsps_trackstart", &_dist_firstsps_trackstart, "dist_firstsps_trackstart/D");
        _tree->Branch("dist_firstsps_trackend", &_dist_firstsps_trackend, "dist_firstsps_trackend/D");
        _tree->Branch("track_true_E", &_track_true_E, "track_true_E/D");
        _tree->Branch("mct_MCS_E", &_mct_MCS_E, "mct_MCS_E/D");
        _tree->Branch("track_MCS_E", &_track_MCS_E, "track_MCS_E/D");
        _tree->Branch("sps_MCS_E", &_sps_MCS_E, "sps_MCS_E/D");
        _tree->Branch("trk_len", &_trk_len, "trk_len/D");
        _tree->Branch("sps_trk_len", &_sps_trk_len, "sps_trk_len/D");

        return true;
    }

    bool TestTrkSPSAssn::analyze(storage_manager* storage) {


        auto ev_mctrack = storage->get_data<event_mctrack>("mcreco");
        if (!ev_mctrack) {
            print(larlite::msg::kERROR, __FUNCTION__, Form("Did not find specified data product, mctrack!"));
            return false;
        }
        if (ev_mctrack->size() != 1) {
            // print(larlite::msg::kERROR, __FUNCTION__, Form("# of MCTracks in this event is not 1. Run this on single muons!"));
            return false;
        }

        auto ev_track = storage->get_data<event_track>(_track_producer);
        if (!ev_track) {
            print(larlite::msg::kERROR, __FUNCTION__, Form("Did not find specified data product, track!"));
            return false;
        }
        if (!ev_track->size()) {
            // print(larlite::msg::kERROR, __FUNCTION__, Form("Zero reconstructed tracks in this event!"));
            return false;
        }

        if ( !_fidvol.Contain(ev_mctrack->at(0).front().Position().Vect()) || !_fidvol.Contain(ev_mctrack->at(0).back().Position().Vect()) ) {
            // print(larlite::msg::kERROR, __FUNCTION__, Form("MC track not fully contained in fiducial volume."));
            return false;
        }
        if ( !_fidvol.Contain(ev_track->at(0).Vertex()) || !_fidvol.Contain(ev_track->at(0).End()) ) {
            // print(larlite::msg::kERROR, __FUNCTION__, Form("Reco track not fully contained in fiducial volume."));
            return false;
        }

        _trk_len = (ev_track->at(0).End() - ev_track->at(0).Vertex()).Mag();
        if (_trk_len < 100.) return false;

        _n_track_traj_pts = ev_track->at(0).NumberTrajectoryPoints();

        event_spacepoint* ev_sps = nullptr;
        auto const& ass_sps_v = storage->find_one_ass(ev_track->id(),
                                ev_sps);

        if (!ev_sps) {
            std::cerr << "no event_sps!" << std::endl;
            return false;
        }

        if (ev_sps->empty()) {
            std::cout << "ev_sps empty" << std::endl;
            return false;
        }

        if (!ass_sps_v.size()) {
            std::cout << "ass_sps_v no size??" << std::endl;
            return false;
        }
        if (!ass_sps_v.at(0).size()) {
            std::cout << "ass_sps_v.at(0) no size??" << std::endl;
            return false;
        }
        _n_total_sps = ev_sps->size();
        // std::cout << "ev_hit size is " << _n_tothits << std::endl;
        // std::cout << "ev_hit id is std::pair<"
        //           << ev_hit->id().first << ","
        //           << ev_hit->id().second.c_str()
        //           << ">" << std::endl;
        _dist_firstsps_trackstart = (TVector3(ev_sps->at(ass_sps_v.at(0).front()).XYZ()[0],
                                              ev_sps->at(ass_sps_v.at(0).front()).XYZ()[1],
                                              ev_sps->at(ass_sps_v.at(0).front()).XYZ()[2]) -
                                     ev_track->at(0).Vertex() ).Mag();
        _dist_firstsps_trackend = (TVector3(ev_sps->at(ass_sps_v.at(0).front()).XYZ()[0],
                                            ev_sps->at(ass_sps_v.at(0).front()).XYZ()[1],
                                            ev_sps->at(ass_sps_v.at(0).front()).XYZ()[2]) -
                                   ev_track->at(0).End() ).Mag();

        bool sps_flipped = _dist_firstsps_trackend < _dist_firstsps_trackstart;

        // std::cout<<"dist first sps to track start is "<<_dist_firstsps_trackstart<<std::endl;
        // std::cout<<"dist first sps to track end is "<<_dist_firstsps_trackend<<std::endl;

        // std::cout << "sps flipped is " << sps_flipped << std::endl;
        size_t reco_track_idx = 0; //since there is only one track, otherwise choose the index in ev_track you care about

        // get the hits associated with this track
        _n_asstd_sps = ass_sps_v.at(reco_track_idx).size();

        // std::cout << "Reco Track start Z is " << ev_track->at(0).Vertex().Z() << std::endl;
        // std::cout << "Reco Track end Z is " << ev_track->at(0).End().Z() << std::endl;
        // std::cout << "Reco track has " << _n_track_traj_pts << " trajectory points." << std::endl;

        //  std::cout << "The first ten (of "<<_n_track_traj_pts<<" total) traj points on the track are:" << std::endl;
        // for (size_t i = 0; i < 10; ++i)
        //     std::cout << "(" << ev_track->at(0).LocationAtPoint(i).X()
        //               << ", " << ev_track->at(0).LocationAtPoint(i).Y()
        //               << ", " << ev_track->at(0).LocationAtPoint(i).Z() << ")" << std::endl;
        // size_t counter = 0;

        // std::cout << "The first ten associated spacepoints (of "<<_n_asstd_sps<<" total) are:" << std::endl;


        /// Loop over the spacepoints and create my own new track:
        larlite::track sps_track;
        sps_track.set_track_id(999);//trk.ID());

        if (!sps_flipped) {
            for (int i = 0; i < ass_sps_v.at(reco_track_idx).size(); ++i) {
                auto const& asstd_sps_idx = ass_sps_v.at(reco_track_idx).at(i);
                const larlite::spacepoint isps = ev_sps->at(asstd_sps_idx);
                sps_track.add_vertex(TVector3(isps.XYZ()[0], isps.XYZ()[1], isps.XYZ()[2]));
                sps_track.add_direction(TVector3(isps.XYZ()[0], isps.XYZ()[1], isps.XYZ()[2]));
            }
        }
        else {
            for (int i = ass_sps_v.at(reco_track_idx).size() - 1; i >= 0; --i) {
                auto const& asstd_sps_idx = ass_sps_v.at(reco_track_idx).at(i);
                const larlite::spacepoint isps = ev_sps->at(asstd_sps_idx);
                sps_track.add_vertex(TVector3(isps.XYZ()[0], isps.XYZ()[1], isps.XYZ()[2]));
                sps_track.add_direction(TVector3(isps.XYZ()[0], isps.XYZ()[1], isps.XYZ()[2]));
            }

        }

        _sps_trk_len = (sps_track.End() - sps_track.Vertex()).Mag();
        // std::cout << "Created Track from SPS start Z is " << ev_track->at(0).Vertex().Z() << std::endl;
        // std::cout << "Created Track from SPS end Z is " << ev_track->at(0).End().Z() << std::endl;
        // std::cout << "Created Track from SPS has " << sps_track.NumberTrajectoryPoints() << " trajectory points." << std::endl;
        //  std::cout << "The first ten (of "<<_n_track_traj_pts<<" total) traj points on the track are:" << std::endl;
        // for (size_t i = 0; i < 10; ++i)
        //     std::cout << "(" << sps_track.LocationAtPoint(i).X()
        //               << ", " << sps_track.LocationAtPoint(i).Y()
        //               << ", " << sps_track.LocationAtPoint(i).Z() << ")" << std::endl;
        // for (auto const& asstd_sps_idx : ass_sps_v.at(reco_track_idx))
        // {
        //     // counter++;

        //     // const larlite::spacepoint isps = ev_sps->at(asstd_sps_idx);
        //     // std::cout << "(" << isps.XYZ()[0]
        //     //           << ", " << isps.XYZ()[1]
        //     //           << ", " << isps.XYZ()[2] << ")" << std::endl;
        //     // if (counter >= 10) break;

        //     // std::cout << "   this sps has z position: " << isps.XYZ()[2] << std::endl;

        // }




        // _mct_depE = ev_mctrack->at(0).front().E() - ev_mctrack->at(0).back().E();

        // std::cout << "MCTrack truth start kinetic energy is roughly : " << ev_mctrack->at(0).front().E() - 106. << std::endl;
        // std::cout << "MCTrack deposited energy is : " << ev_mctrack->at(0).front().E() - ev_mctrack->at(0).back().E() << " MEV." << std::endl;
        // std::cout << "Summed ADC of associated hits is : "<<_sum_hit_ADC<<std::endl;



        // _mct_MCS_E = _tmc->GetMomentumMultiScatterLLHD(ev_mctrack->at(0));
        _track_true_E = ev_mctrack->at(0).front().E() / 1000.;
        if (_sps1track0 == 0)
            _track_MCS_E = _tmc->GetMomentumMultiScatterLLHD(ev_track->at(0), false, true);
        else
            _sps_MCS_E = _tmc->GetMomentumMultiScatterLLHD(sps_track, false, true);



        // if (_sps_trk_len - _trk_len < -100) {
        //     std::cout << "found a fucking weird SPS track. let's print out info:" << std::endl;
        //     std::cout << "Reco track has " << _n_track_traj_pts << " trajectory points." << std::endl;

        //     std::cout << "The first ten (of " << _n_track_traj_pts << " total) traj points on the track are:" << std::endl;
        //     for (size_t i = 0; i < 10; ++i)
        //         std::cout << "(" << ev_track->at(0).LocationAtPoint(i).X()
        //                   << ", " << ev_track->at(0).LocationAtPoint(i).Y()
        //                   << ", " << ev_track->at(0).LocationAtPoint(i).Z() << ")" << std::endl;

        //        std::cout << "The last ten (of " << _n_track_traj_pts << " total) traj points on the track are:" << std::endl;
        //     for (size_t i = _n_track_traj_pts - 1; i > _n_track_traj_pts - 10; --i)
        //         std::cout << "(" << ev_track->at(0).LocationAtPoint(i).X()
        //                   << ", " << ev_track->at(0).LocationAtPoint(i).Y()
        //                   << ", " << ev_track->at(0).LocationAtPoint(i).Z() << ")" << std::endl;

        //     std::cout << "The first ten (of " << sps_track.NumberTrajectoryPoints() << " total) traj points on the track are:" << std::endl;
        //     for (size_t i = 0; i < 10; ++i)
        //         std::cout << "(" << sps_track.LocationAtPoint(i).X()
        //                   << ", " << sps_track.LocationAtPoint(i).Y()
        //                   << ", " << sps_track.LocationAtPoint(i).Z() << ")" << std::endl;

        // std::cout << "The last ten (of " << sps_track.NumberTrajectoryPoints() << " total) traj points on the track are:" << std::endl;
        //     for (size_t i = sps_track.NumberTrajectoryPoints() - 1; i > sps_track.NumberTrajectoryPoints() - 10; --i)
        //         std::cout << "(" << sps_track.LocationAtPoint(i).X()
        //                   << ", " << sps_track.LocationAtPoint(i).Y()
        //                   << ", " << sps_track.LocationAtPoint(i).Z() << ")" << std::endl;
        // }



        _tree->Fill();

        return true;
    }

    bool TestTrkSPSAssn::finalize() {

        if (_fout) {
            _fout->cd();
            _tree->Write();
            _tmc->GetTree()->Write();
        }

        return true;
    }

}
#endif
