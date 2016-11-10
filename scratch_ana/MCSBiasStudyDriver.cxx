#ifndef LARLITE_MCSBIASSTUDYDRIVER_CXX
#define LARLITE_MCSBIASSTUDYDRIVER_CXX

#include "MCSBiasStudyDriver.h"

namespace larlite {

    bool MCSBiasStudyDriver::initialize() {

        if (_ana_type == MCSBiasStudy::kANALYSIS_TYPE_MAX) {
            print(larlite::msg::kERROR, __FUNCTION__,
                  Form("Did not set what kind of analysis you are doing in MCSBiasStudyDriver!!! Bahhh!"));
            return false;
        }


        if (!_module) {
            _module = new MCSBiasStudy();
            _module->SetAnalysisType(_ana_type);
        }

        if (!_driver_tree) {
            _driver_tree = new TTree("driver_tree", "driver_tree");
            _driver_tree->Branch("MCT_PDG", &_MCT_PDG, "MCT_PDG/I");
            _driver_tree->Branch("MCT_origin", &_MCT_origin, "MCT_origin/I");
            _driver_tree->Branch("run", &_run, "run/I");
            _driver_tree->Branch("subrun", &_subrun, "subrun/I");
            _driver_tree->Branch("eventid", &_eventid, "eventid/I");
        }

        counter = 0;

        return true;
    }

    bool MCSBiasStudyDriver::analyze(storage_manager* storage) {

        larlite::event_mctrack *ev_mctrack;
        larlite::mctrack mct;
        larlite::event_track *ev_track;
        larlite::track track;
        larlite::event_vertex *ev_vertex;
        larlite::vertex vtx;
        _run = 999999;
        _subrun = 999999;
        _eventid = 999999;
        _MCT_PDG = -99999;
        _MCT_origin = -99999;
        bool flip = false;

        // Different cases for different analysis types!
        // std::cout << "_ana type is " << _ana_type << std::endl;


        ////////////////////////////////////////////////////////////////////////////////////
        /// Simulated single muons. Using MCTracks.
        ////////////////////////////////////////////////////////////////////////////////////
        if (_ana_type == MCSBiasStudy::kSingleMuonMCTrack)
        {
            // std::cout << "single muon case!" << std::endl;
            ev_mctrack = storage->get_data<event_mctrack>("mcreco");
            if (!ev_mctrack) {
                print(larlite::msg::kERROR, __FUNCTION__, Form("Did not find specified data product, mctrack!"));
                return false;
            }
            if (ev_mctrack->size() != 1)
                return false;


            mct = ev_mctrack->at(0);
            _MCT_PDG = mct.PdgCode();
            _MCT_origin = mct.Origin();

            _run = ev_mctrack->run();
            _subrun = ev_mctrack->subrun();
            _eventid = ev_mctrack->event_id();

            if (mct.size() < 3)
                return false;


            if ( (mct.back().Position().Vect() - mct.front().Position().Vect()).Mag() < 100.) return false;


            if (!_fidvol.Contain(mct.front().Position().Vect()) || !_fidvol.Contain(mct.back().Position().Vect()))
                return false;

            //There are a few MCTrack muons that decay in flight, so range energy doesn't match MCS energy
            //I found ~2 of these in ~700 events
            // let's throw these out too. they have End().E() = 105.658367
            if (mct.End().E() > 106.) return false;

            // Now we are left with the MCTracks we actually care about
            _module->AnalyzeTrack(mct, _run, _subrun, _eventid);

            ///////////////////////////////////////////
            // END kSingleMuonMCTrack case
            ///////////////////////////////////////////
        }

        ////////////////////////////////////////////////////////////////////////////////////
        /// Simulated single muons. Using reconstructed tracks (pandoraNuPMA)
        /// MCTracks are also read in to determine track direction.
        ////////////////////////////////////////////////////////////////////////////////////
        else if (_ana_type == MCSBiasStudy::kSingleMuonRecoTrack)
        {
            // std::cout << "single muon reco case!" << std::endl;
            ev_track = storage->get_data<event_track>("pandoraNuPMA");
            if (!ev_track) {
                print(larlite::msg::kERROR, __FUNCTION__, Form("Did not find specified data product, track!"));
                return false;
            }

            if (ev_track->size() != 1) {
                //std::cout << "Number of reco tracks in event != 1. skipping" << std::endl;
                //6052 / 19500 make it past this part
                return false;
            }

            ev_mctrack = storage->get_data<event_mctrack>("mcreco");
            if (!ev_mctrack) {
                print(larlite::msg::kERROR, __FUNCTION__, Form("Did not find specified data product, mctrack!"));
                return false;
            }
            if (ev_mctrack->size() != 1)
                return false;

            mct = ev_mctrack->at(0);
            _run = ev_mctrack->run();
            _subrun = ev_mctrack->subrun();
            _eventid = ev_mctrack->event_id();

            if (mct.size() < 3) return false;
            if ( (mct.back().Position().Vect() - mct.front().Position().Vect()).Mag() < 100.) return false;
            if (!_fidvol.Contain(mct.front().Position().Vect()) || !_fidvol.Contain(mct.back().Position().Vect())) return false;

            //There are a few MCTrack muons that decay in flight, so range energy doesn't match MCS energy
            //I found ~2 of these in ~700 events
            // let's throw these out too. they have End().E() = 105.658367
            if (mct.End().E() > 106.) return false;

            // 1256 / 6052 make it here
            track = ev_track->at(0);
            if ( !_fidvol.Contain(track.Vertex()) || !_fidvol.Contain(track.End()) ) return false;

            // 1237 / 1256 make it here

            // Require the track is well reconstructed (Start and end within 3cm of the real thing)
            auto const mct_start = mct.front().Position().Vect();
            auto const mct_end   = mct.back().Position().Vect();
            auto const trk_start = track.Vertex();
            auto const trk_end   = track.End();


            // Track well recod and in the correct direction
            if ( (mct_start - trk_start).Mag() < 3 && (mct_end - trk_end).Mag() < 3 )
                flip = false;
            // Track well recod and in the flipped direction
            else if ( (mct_start - trk_end).Mag() < 3 && (mct_end - trk_start).Mag() < 3 )
                flip = true;
            // Track badly recod
            else return false;

            // 928 / 1237 make it here (928 are well recod, though maybe backwards)

            _module->AnalyzeTrack(track, _run, _subrun, _eventid, flip);


            ///////////////////////////////////////////
            // END kSingleMuonRecoTrack case
            ///////////////////////////////////////////
        }

        ////////////////////////////////////////////////////////////////////////////////////
        /// Simulated BNB+Cosmic events. One reco track is saved per event to the input file
        /// Also one reconstructed vertex is saved to the input file (to determine track direction)
        ////////////////////////////////////////////////////////////////////////////////////
        else if (_ana_type == MCSBiasStudy::kMCBNBSelectedRecoTrack)
        {
            // std::cout << "mc bnb selected reco track case!" << std::endl;
            // print(larlite::msg::kERROR, __FUNCTION__, Form("MCSBiasStudy::kMCBNBSelectedRecoTrack not yet implemented!"));
            // return false;

            ev_track = storage->get_data<event_track>("savedlongesttrack");
            if (!ev_track) {
                print(larlite::msg::kERROR, __FUNCTION__, Form("Did not find specified data product, track!"));
                return false;
            }

            if (ev_track->size() != 1) {
                std::cout << "Number of reco tracks in event != 1. skipping" << std::endl;
                //6052 / 19500 make it past this part
                return false;
            }

            ev_mctrack = storage->get_data<event_mctrack>("mcreco");
            if (!ev_mctrack) {
                print(larlite::msg::kERROR, __FUNCTION__, Form("Did not find specified data product, mctrack!"));
                return false;
            }

            _run = ev_mctrack->run();
            _subrun = ev_mctrack->subrun();
            _eventid = ev_mctrack->event_id();

            track = ev_track->at(0);
            // Make sure the track is contained in the fiducial volume
            if ( !_fidvol.Contain(track.Vertex()) || !_fidvol.Contain(track.End()) ) return false;

            // Make sure the track is at least 1m long
            if ( (track.End() - track.Vertex()).Mag() < 100. ) return false;

            // Find the mctrack that matches the reco track, if there is one
            bool found_mctrack = false;

            auto const trk_start = track.Vertex();
            auto const trk_end   = track.End();
            for (size_t i = 0; i < ev_mctrack->size(); ++i) {
                mct = ev_mctrack->at(i);
                if (mct.size() < 3) continue;
                auto const mct_start = mct.front().Position().Vect();
                auto const mct_end   = mct.back().Position().Vect();

                // Track well recod and in the correct direction
                if ( (mct_start - trk_start).Mag() < 3 && (mct_end - trk_end).Mag() < 3 ) {
                    found_mctrack = true;
                    flip = false;
                    break;
                }
                // Track well recod and in the flipped direction
                else if ( (mct_start - trk_end).Mag() < 3 && (mct_end - trk_start).Mag() < 3 ) {
                    found_mctrack = true;
                    flip = true;
                    break;
                }
            }

            // If there was no good MCTrack matching the reco track
            if ( !found_mctrack ) return false;

            _MCT_PDG = mct.PdgCode();
            _MCT_origin = mct.Origin();

            _module->AnalyzeTrack(track, _run, _subrun, _eventid, flip);

            ///////////////////////////////////////////
            // END kMCBNBSelectedRecoTrack case
            ///////////////////////////////////////////
        }

        ////////////////////////////////////////////////////////////////////////////////////
        /// Actual data BNB selected events. One reco track is saved per event to the input file
        /// Also one reconstructed vertex is saved to the input file (to determine track direction)
        /// Any hand-scanning information is merged in in later analysis, not in this code.
        ////////////////////////////////////////////////////////////////////////////////////
        else if (_ana_type == MCSBiasStudy::kDataBNBSelectedRecoTrack)
        {

            // std::cout << "data bnb case!" << std::endl;
            ev_track = storage->get_data<event_track>("savedlongesttrack");
            if (!ev_track) {
                print(larlite::msg::kERROR, __FUNCTION__, Form("Did not find specified data product, track!"));
                return false;
            }
            if (!ev_track->size())
                return false;
            if (ev_track->size() != 1) {
                std::cout << "More than 1 saved reco track in this event? skipping" << std::endl;
                return false;
            }

            ev_vertex = storage->get_data<event_vertex>("savedvertex");
            if (!ev_vertex) {
                print(larlite::msg::kERROR, __FUNCTION__, Form("Did not find specified data product, vertex!"));
                return false;
            }
            if (!ev_vertex->size())
                return false;
            if (ev_vertex->size() != 1) {
                std::cout << "More than 1 saved vertex in this event? skipping" << std::endl;
                return false;
            }

            _run = ev_track->run();
            _subrun = ev_track->subrun();
            _eventid = ev_track->event_id();

            track = ev_track->at(0);

            // Determine if the track needs to be flipped (force it to start at vertex)
            vtx = ev_vertex->at(0);
            TVector3 vtx_vec = TVector3(vtx.X(), vtx.Y(), vtx.Z());
            flip = (track.Vertex() - vtx_vec).Mag2() > (track.End() - vtx_vec).Mag2();
            _module->AnalyzeTrack(track, _run, _subrun, _eventid, flip);

            ///////////////////////////////////////////
            // END kDataBNBSelectedRecoTrack case
            ///////////////////////////////////////////
        }

        ////////////////////////////////////////////////////////////////////////////////////
        /// The user did not set an analysis type. Return error message.
        ////////////////////////////////////////////////////////////////////////////////////
        else
        {
            print(larlite::msg::kERROR, __FUNCTION__,
                  Form("Did not set what kind of analysis you are doing in MCSBiasStudyDriver!!! Bahhh!"));
            return false;
        }

        _driver_tree->Fill();

        return true;
    }

    bool MCSBiasStudyDriver::finalize() {

        std::cout << "counter = " << counter << std::endl;
        if (_fout) {
            _fout->cd();
            _module->GetTree()->Write();
            // _module->GetSegTree()->Write();
            _module->GetTMCTree()->Write();
            _driver_tree->Write();
        }

        return true;
    }

}
#endif
