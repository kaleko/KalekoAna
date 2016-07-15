#ifndef LARLITE_TRACKDATAMCCOMP_CXX
#define LARLITE_TRACKDATAMCCOMP_CXX

#include "TrackDataMCComp.h"
#include "DataFormat/mctrack.h"
#include "DataFormat/track.h"

namespace larlite {

    bool TrackDataMCComp::initialize() {

        if (_track_producer.empty()) {
            std::cout << "Track producer empty. Quitting." << std::endl;
            return false;
        }

        if (!_ana_tree) {
            _ana_tree = new TTree("evt_tree", "evt_tree");
            _ana_tree->Branch("n_mctracks_cosmic", &_n_mctracks_cosmic, "n_mctracks_cosmic/I");
            _ana_tree->Branch("n_mctracks_neutrino", &_n_mctracks_neutrino, "n_mctracks_neutrino/I");
            _ana_tree->Branch("n_recotracks", &_n_recotracks, "n_recotracks/I");
        }

        if (!_track_tree) {
            _track_tree = new TTree("trk_tree", "trk_tree");
            _track_tree->Branch("trk_len", &_trk_len, "trk_len/D");
            _track_tree->Branch("trk_start_x", &_trk_start_x, "trk_start_x/D");
            _track_tree->Branch("trk_start_y", &_trk_start_y, "trk_start_y/D");
            _track_tree->Branch("trk_start_z", &_trk_start_z, "trk_start_z/D");
            _track_tree->Branch("trk_end_x", &_trk_end_x, "trk_end_x/D");
            _track_tree->Branch("trk_end_y", &_trk_end_y, "trk_end_y/D");
            _track_tree->Branch("trk_end_z", &_trk_end_z, "trk_end_z/D");
        }

        if (!_mctrack_tree) {
            _mctrack_tree = new TTree("mctrk_tree", "mctrk_tree");
            _mctrack_tree->Branch("mctrk_len", &_mctrk_len, "mctrk_len/D");
            _mctrack_tree->Branch("mctrk_origin", &_mctrk_origin, "mctrk_origin/I");
            _mctrack_tree->Branch("mctrk_start_x", &_mctrk_start_x, "mctrk_start_x/D");
            _mctrack_tree->Branch("mctrk_start_y", &_mctrk_start_y, "mctrk_start_y/D");
            _mctrack_tree->Branch("mctrk_start_z", &_mctrk_start_z, "mctrk_start_z/D");
            _mctrack_tree->Branch("mctrk_end_x", &_mctrk_end_x, "mctrk_end_x/D");
            _mctrack_tree->Branch("mctrk_end_y", &_mctrk_end_y, "mctrk_end_y/D");
            _mctrack_tree->Branch("mctrk_end_z", &_mctrk_end_z, "mctrk_end_z/D");
        }

        return true;
    }

    bool TrackDataMCComp::analyze(storage_manager* storage) {

        // Initialize event-by-event TTree vars
        _n_mctracks_neutrino = 0;
        _n_mctracks_cosmic = 0;
        _n_recotracks = 0;

        // MCTrack info first
        if (!_running_on_data) {

            auto ev_mctrack = storage->get_data<event_mctrack>("mcreco");
            if (!ev_mctrack) {
                print(larlite::msg::kERROR, __FUNCTION__, Form("Did not find specified data product, mctrack!"));
                return false;
            }
            if (!ev_mctrack->size()) {
                //print(larlite::msg::kERROR, __FUNCTION__, Form("Zero mctracks in this event!"));
                return false;
            }

            for (auto const &mct : *ev_mctrack) {

                if ( !mct.size() ) continue;

                _mctrk_len = (mct.back().Position().Vect() - mct.front().Position().Vect()).Mag();
                if ( _mctrk_len < 3.) continue;

                if ( mct.Origin() == 1 ) _n_mctracks_neutrino++;
                else _n_mctracks_cosmic++;

                _mctrk_origin = mct.Origin();
                _mctrk_start_x = mct.front().Position().Vect().X();
                _mctrk_start_y = mct.front().Position().Vect().Y();
                _mctrk_start_z = mct.front().Position().Vect().Z();
                _mctrk_end_x = mct.back().Position().Vect().X();
                _mctrk_end_y = mct.back().Position().Vect().Y();
                _mctrk_end_z = mct.back().Position().Vect().Z();

                _mctrack_tree->Fill();

            } // End loop over MCTracks

        } // End if not running on data

        // Now Reco Track info

        auto ev_track = storage->get_data<event_track>(_track_producer);
        if (!ev_track) {
            print(larlite::msg::kERROR, __FUNCTION__, Form("Did not find specified data product, track!"));
            return false;
        }
        if (!ev_track->size()) {
            //print(larlite::msg::kERROR, __FUNCTION__, Form("Zero reconstructed tracks in this event!"));
            return false;
        }

        for (auto const &trk : *ev_track) {

            if ( trk.NumberTrajectoryPoints() < 2 ) continue;

            _trk_len = (trk.End() - trk.Vertex()).Mag();
            if ( _trk_len < 3.) continue;

            _n_recotracks++;

            _trk_start_x = trk.Vertex().X();
            _trk_start_y = trk.Vertex().Y();
            _trk_start_z = trk.Vertex().Z();
            _trk_end_x = trk.End().X();
            _trk_end_y = trk.End().Y();
            _trk_end_z = trk.End().Z();

            _track_tree->Fill();

        }

        _ana_tree->Fill();

        return true;

    }

    bool TrackDataMCComp::finalize() {

        if (_fout) {
            _ana_tree->Write();
            _track_tree->Write();
            _mctrack_tree->Write();

        }
        return true;
    }

}
#endif
