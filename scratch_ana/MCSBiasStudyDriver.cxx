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

        _stat_decrease_counters.clear();
        for (size_t i = 0; i < 10; ++i) _stat_decrease_counters.push_back(0);

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
            _stat_decrease_counters.at(0)++;
            // std::cout << "single muon case!" << std::endl;
            ev_mctrack = storage->get_data<event_mctrack>("mcreco");
            if (!ev_mctrack) {
                print(larlite::msg::kERROR, __FUNCTION__, Form("Did not find specified data product, mctrack!"));
                return false;
            }
            if (ev_mctrack->size() != 1)
                return false;

            _stat_decrease_counters.at(1)++;

            mct = ev_mctrack->at(0);
            _MCT_PDG = mct.PdgCode();
            _MCT_origin = mct.Origin();

            _run = ev_mctrack->run();
            _subrun = ev_mctrack->subrun();
            _eventid = ev_mctrack->event_id();

            if (mct.size() < 3)
                return false;

            _stat_decrease_counters.at(2)++;


            if ( (mct.back().Position().Vect() - mct.front().Position().Vect()).Mag() < 100.) return false;

            _stat_decrease_counters.at(3)++;

            if (!_fidvol.Contain(mct.front().Position().Vect()) || !_fidvol.Contain(mct.back().Position().Vect()))
                return false;

            _stat_decrease_counters.at(4)++;

            //There are a few MCTrack muons that decay in flight, so range energy doesn't match MCS energy
            //I found ~2 of these in ~700 events
            // let's throw these out too. they have End().E() = 105.658367
            if (mct.End().E() > 106.) return false;

            _stat_decrease_counters.at(5)++;

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
            _stat_decrease_counters.at(0)++;

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

            _stat_decrease_counters.at(1)++;

            ev_mctrack = storage->get_data<event_mctrack>("mcreco");
            if (!ev_mctrack) {
                print(larlite::msg::kERROR, __FUNCTION__, Form("Did not find specified data product, mctrack!"));
                return false;
            }
            if (ev_mctrack->size() != 1)
                return false;

            _stat_decrease_counters.at(2)++;

            mct = ev_mctrack->at(0);
            _run = ev_mctrack->run();
            _subrun = ev_mctrack->subrun();
            _eventid = ev_mctrack->event_id();

            if (mct.size() < 3) return false;
            _stat_decrease_counters.at(3)++;
            if ( (mct.back().Position().Vect() - mct.front().Position().Vect()).Mag() < 100.) return false;
            _stat_decrease_counters.at(4)++;
            if (!_fidvol.Contain(mct.front().Position().Vect()) || !_fidvol.Contain(mct.back().Position().Vect())) return false;
            _stat_decrease_counters.at(5)++;

            //There are a few MCTrack muons that decay in flight, so range energy doesn't match MCS energy
            //I found ~2 of these in ~700 events
            // let's throw these out too. they have End().E() = 105.658367
            if (mct.End().E() > 106.) return false;
            _stat_decrease_counters.at(6)++;

            // 1256 / 6052 make it here
            track = ev_track->at(0);
            if ( !_fidvol.Contain(track.Vertex()) || !_fidvol.Contain(track.End()) ) return false;
            _stat_decrease_counters.at(7)++;

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
            _stat_decrease_counters.at(8)++;

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
            _stat_decrease_counters.at(0)++;

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
            _stat_decrease_counters.at(1)++;

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
            _stat_decrease_counters.at(2)++;
            // Make sure the track is at least 1m long
            if ( (track.End() - track.Vertex()).Mag() < 100. ) return false;
            _stat_decrease_counters.at(3)++;

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
            _stat_decrease_counters.at(4)++;

            _MCT_PDG = mct.PdgCode();
            _MCT_origin = mct.Origin();

            _module->AnalyzeTrack(track, _run, _subrun, _eventid, flip);

            ///////////////////////////////////////////
            // END kMCBNBSelectedRecoTrack case
            ///////////////////////////////////////////
        }

        ////////////////////////////////////////////////////////////////////////////////////
        /// Simulated BNB only events. They are pre-selected to be numuCC interaction inside of the
        /// fiducial volume with one long muon MC track that is also fully contained
        /// Here we analyze the reconstructed track that matches the MCTrack.
        ////////////////////////////////////////////////////////////////////////////////////
        else if (_ana_type == MCSBiasStudy::kMCBNBRecoTrack)
        {
            _stat_decrease_counters.at(0)++;

            ev_track = storage->get_data<event_track>("pandoraNuPMA");
            if (!ev_track) {
                print(larlite::msg::kERROR, __FUNCTION__, Form("Did not find specified data product, track!"));
                return false;
            }

            if (!ev_track->size()) {
                // std::cout << "No reco tracks in event :( Skipping..." << std::endl;
                return false;
            }
            _stat_decrease_counters.at(1)++;

            // Grab the MCTruth
            auto ev_mctruth = storage->get_data<event_mctruth>("generator");
            if (!ev_mctruth) {
                print(larlite::msg::kERROR, __FUNCTION__, Form("Did not find specified data product, mctruth!"));
                return false;
            }
            if (ev_mctruth->size() != 1) {
                // Sometimes size is 2 if there are two neutrinos... let's just throw out these events since we
                // have high stats and I don't feel like writing the code to handle them
                if (ev_mctruth->size() == 2) return false;

                // If the size is 0 or more than 2 something is wrong.
                print(larlite::msg::kERROR, __FUNCTION__,
                      Form("MCTruth size is not 1! More than two neutrinos? Cosmics? Size is %zu!", ev_mctruth->size())
                     );
                return false;
            }
            _stat_decrease_counters.at(2)++;

            auto const &mctruth = ev_mctruth->at(0);

            // Make sure the event is numuCC inside of the fiducial volume

            //Enforce CC interaction channel
            if ( mctruth.GetNeutrino().CCNC() != 0 ) return false;
            _stat_decrease_counters.at(3)++;

            // If neutrino interacts outside of fiducial volume, skip event
            auto const nu_vtx = mctruth.GetNeutrino().Nu().Trajectory().back().Position().Vect();
            if (!_fidvol.Contain(nu_vtx)) return false;
            _stat_decrease_counters.at(4)++;

            // If neutrino is not a numu, skip event
            if (mctruth.GetNeutrino().Nu().PdgCode() != 14) return false;
            _stat_decrease_counters.at(5)++;


            ev_mctrack = storage->get_data<event_mctrack>("mcreco");
            if (!ev_mctrack) {
                print(larlite::msg::kERROR, __FUNCTION__, Form("Did not find specified data product, mctrack!"));
                return false;
            }
            if (!ev_mctrack->size()) {
                std::cout << "No MCTracks in event :( Skipping..." << std::endl;
                return false;
            }
            _stat_decrease_counters.at(6)++;

            _run = ev_mctrack->run();
            _subrun = ev_mctrack->subrun();
            _eventid = ev_mctrack->event_id();


            larlite::mctrack the_mctrack;
            // Find the mctrack that is a muon starting from the neutrino interaction
            // We also require it is 1m long at least... might want to remove this if you are using this
            // filter for some reason other than a multiple coloumb scattering analysis which requires 1m
            size_t n_found_MCTs = 0;
            for (auto const& mct : *ev_mctrack) {

                // Sometimes mctracks have zero size. No idea why. Skip them.
                if ( mct.size() < 3 ) continue;

                // origin == 1 means comes from neutrino interaction (IE not cosmic)
                if (mct.Origin() != 1 ) continue;

                // MCTrack has to be truly a muon
                if (mct.PdgCode() != 13 ) continue;

                //MCTrack should start VERY CLOSE to nu vtx
                // (we're talking like 10^-28 or smaller ... below floating point precision)
                if ( (mct.front().Position().Vect() - nu_vtx).Mag2() > 0.0001) continue;

                // Enforce the muon is fully contained in fiducial volume.
                // (note we already checked the front of it is close to the nu vtx, and we already checked
                // the nu vtx is in the fidicial volume... so we only need to check the back here)
                // But let's just check the front anyway
                if (!_fidvol.Contain(mct.back().Position().Vect()) || !_fidvol.Contain(mct.front().Position().Vect()) ) continue;

                // Require the muon is at least one meter in length.
                if ((mct.back().Position().Vect() - mct.front().Position().Vect()).Mag() < 100.) continue;

                // Make sure the muon doesn't decay in flight
                if (mct.End().E() > 106.) return false;

                // Save this mctrack!
                the_mctrack = mct;
                // At this point you have found a viable muon!
                n_found_MCTs++;

            }

            // If you didn't find any viable muons, skip the event
            if (!n_found_MCTs) return false;
            _stat_decrease_counters.at(7)++;

            // If somehow you found more than one muon, something has gone wrong. Skip the event.
            if (n_found_MCTs > 1) {
                print(larlite::msg::kWARNING, __FUNCTION__, Form("More than one viable muon in the event?!"));
                return false;
            }
            _stat_decrease_counters.at(8)++;


            // Now "the_mctrack" is the one MCTrack we care about. It should have length > 1 and be
            // the muon from a true numuCC interaction
            auto const mct_start = the_mctrack.front().Position().Vect();
            auto const mct_end   = the_mctrack.back().Position().Vect();

            bool found_track = false;

            for (auto const& trk : *ev_track) {

                // Make sure the track is contained in the fiducial volume
                if ( !_fidvol.Contain(trk.Vertex()) || !_fidvol.Contain(trk.End()) ) continue;

                // Make sure the track is at least 1m long
                if ( (trk.End() - trk.Vertex()).Mag() < 100. ) continue;


                // Check if this track matches the already-found mctrack we care about ("the_mctrack")
                auto const trk_start = trk.Vertex();
                auto const trk_end   = trk.End();

                // Track well recod and in the correct direction
                if ( (mct_start - trk_start).Mag() < 3 && (mct_end - trk_end).Mag() < 3 ) {
                    flip = false;
                    found_track = true;
                    track = trk;
                    break;
                }
                // Track well recod and in the flipped direction
                else if ( (mct_start - trk_end).Mag() < 3 && (mct_end - trk_start).Mag() < 3 ) {
                    flip = true;
                    found_track = true;
                    track = trk;
                    break;
                }
                else continue;
            }

            // If there was no good reco track matching the MCTrack
            if ( !found_track ) return false;
            _stat_decrease_counters.at(9)++;

            _MCT_PDG = the_mctrack.PdgCode();
            _MCT_origin = the_mctrack.Origin();

            _module->AnalyzeTrack(track, _run, _subrun, _eventid, flip);

            ///////////////////////////////////////////
            // END kMCBNBRecoTrack case
            ///////////////////////////////////////////
        }

        ////////////////////////////////////////////////////////////////////////////////////
        /// Simulated BNB only events. They are pre-selected to be numuCC interaction inside of the
        /// fiducial volume with one long muon MC track that is also fully contained
        /// Here we analyze the MCTrack.
        ////////////////////////////////////////////////////////////////////////////////////
        else if (_ana_type == MCSBiasStudy::kMCBNBMCTrack)
        {
            _stat_decrease_counters.at(0)++;

            // Grab the MCTruth
            auto ev_mctruth = storage->get_data<event_mctruth>("generator");
            if (!ev_mctruth) {
                print(larlite::msg::kERROR, __FUNCTION__, Form("Did not find specified data product, mctruth!"));
                return false;
            }
            if (ev_mctruth->size() != 1) {
                // Sometimes size is 2 if there are two neutrinos... let's just throw out these events since we
                // have high stats and I don't feel like writing the code to handle them
                if (ev_mctruth->size() == 2) return false;

                // If the size is 0 or more than 2 something is wrong.
                print(larlite::msg::kERROR, __FUNCTION__,
                      Form("MCTruth size is not 1! More than two neutrinos? Cosmics? Size is %zu!", ev_mctruth->size())
                     );
                return false;
            }
            _stat_decrease_counters.at(1)++;

            auto const &mctruth = ev_mctruth->at(0);

            // Make sure the event is numuCC inside of the fiducial volume

            //Enforce CC interaction channel
            if ( mctruth.GetNeutrino().CCNC() != 0 ) return false;
            _stat_decrease_counters.at(2)++;

            // If neutrino interacts outside of fiducial volume, skip event
            auto const nu_vtx = mctruth.GetNeutrino().Nu().Trajectory().back().Position().Vect();
            if (!_fidvol.Contain(nu_vtx)) return false;
            _stat_decrease_counters.at(3)++;

            // If neutrino is not a numu, skip event
            if (mctruth.GetNeutrino().Nu().PdgCode() != 14) return false;
            _stat_decrease_counters.at(4)++;


            ev_mctrack = storage->get_data<event_mctrack>("mcreco");
            if (!ev_mctrack) {
                print(larlite::msg::kERROR, __FUNCTION__, Form("Did not find specified data product, mctrack!"));
                return false;
            }
            if (!ev_mctrack->size()) {
                std::cout << "No MCTracks in event :( Skipping..." << std::endl;
                return false;
            }
            _stat_decrease_counters.at(5)++;

            _run = ev_mctrack->run();
            _subrun = ev_mctrack->subrun();
            _eventid = ev_mctrack->event_id();


            larlite::mctrack the_mctrack;
            // Find the mctrack that is a muon starting from the neutrino interaction
            // We also require it is 1m long at least... might want to remove this if you are using this
            // filter for some reason other than a multiple coloumb scattering analysis which requires 1m
            size_t n_found_MCTs = 0;
            for (auto const& mct : *ev_mctrack) {

                // Sometimes mctracks have zero size. No idea why. Skip them.
                if ( mct.size() < 3 ) continue;

                // origin == 1 means comes from neutrino interaction (IE not cosmic)
                if (mct.Origin() != 1 ) continue;

                // MCTrack has to be truly a muon
                if (mct.PdgCode() != 13 ) continue;

                //MCTrack should start VERY CLOSE to nu vtx
                // (we're talking like 10^-28 or smaller ... below floating point precision)
                if ( (mct.front().Position().Vect() - nu_vtx).Mag2() > 0.0001) continue;

                // Enforce the muon is fully contained in fiducial volume.
                // (note we already checked the front of it is close to the nu vtx, and we already checked
                // the nu vtx is in the fidicial volume... so we only need to check the back here)
                // but let's check the front antwaty
                if (!_fidvol.Contain(mct.back().Position().Vect()) || !_fidvol.Contain(mct.front().Position().Vect())) continue;

                // Require the muon is at least one meter in length.
                if ((mct.back().Position().Vect() - mct.front().Position().Vect()).Mag() < 100.) continue;

                // Make sure the muon doesn't decay in flight
                if (mct.End().E() > 106.) return false;

                // Save this mctrack!
                the_mctrack = mct;
                // At this point you have found a viable muon!
                n_found_MCTs++;

            }

            // If you didn't find any viable muons, skip the event
            if (!n_found_MCTs) return false;
            _stat_decrease_counters.at(6)++;

            // If somehow you found more than one muon, something has gone wrong. Skip the event.
            if (n_found_MCTs > 1) {
                print(larlite::msg::kWARNING, __FUNCTION__, Form("More than one viable muon in the event?!"));
                return false;
            }
            _stat_decrease_counters.at(7)++;


            // Now "the_mctrack" is the one MCTrack we care about. It should have length > 1 and be
            // the muon from a true numuCC interaction
            auto const mct_start = the_mctrack.front().Position().Vect();
            auto const mct_end   = the_mctrack.back().Position().Vect();

            _MCT_PDG = the_mctrack.PdgCode();
            _MCT_origin = the_mctrack.Origin();

            // Analyze the MCTrack!
            _module->AnalyzeTrack(the_mctrack, _run, _subrun, _eventid, flip);

            ///////////////////////////////////////////
            // END kMCBNBRecoTrack case
            ///////////////////////////////////////////
        }


        ////////////////////////////////////////////////////////////////////////////////////
        /// Simulated BNB only events. They are pre-selected to be numuCC interaction inside of the
        /// fiducial volume with one long muon MC track that EXITS
        /// Here we analyze the MCTrack.
        ////////////////////////////////////////////////////////////////////////////////////
        else if (_ana_type == MCSBiasStudy::kMCBNBMCTrackExiting)
        {
            _stat_decrease_counters.at(0)++;

            // Grab the MCTruth
            auto ev_mctruth = storage->get_data<event_mctruth>("generator");
            if (!ev_mctruth) {
                print(larlite::msg::kERROR, __FUNCTION__, Form("Did not find specified data product, mctruth!"));
                return false;
            }
            if (ev_mctruth->size() != 1) {
                // Sometimes size is 2 if there are two neutrinos... let's just throw out these events since we
                // have high stats and I don't feel like writing the code to handle them
                if (ev_mctruth->size() == 2) return false;

                // If the size is 0 or more than 2 something is wrong.
                print(larlite::msg::kERROR, __FUNCTION__,
                      Form("MCTruth size is not 1! More than two neutrinos? Cosmics? Size is %zu!", ev_mctruth->size())
                     );
                return false;
            }
            _stat_decrease_counters.at(1)++;

            auto const &mctruth = ev_mctruth->at(0);

            // Make sure the event is numuCC inside of the fiducial volume

            //Enforce CC interaction channel
            if ( mctruth.GetNeutrino().CCNC() != 0 ) return false;
            _stat_decrease_counters.at(2)++;

            // If neutrino interacts outside of fiducial volume, skip event
            auto const nu_vtx = mctruth.GetNeutrino().Nu().Trajectory().back().Position().Vect();
            if (!_fidvol.Contain(nu_vtx)) return false;
            _stat_decrease_counters.at(3)++;

            // If neutrino is not a numu, skip event
            if (mctruth.GetNeutrino().Nu().PdgCode() != 14) return false;
            _stat_decrease_counters.at(4)++;


            ev_mctrack = storage->get_data<event_mctrack>("mcreco");
            if (!ev_mctrack) {
                print(larlite::msg::kERROR, __FUNCTION__, Form("Did not find specified data product, mctrack!"));
                return false;
            }
            if (!ev_mctrack->size()) {
                std::cout << "No MCTracks in event :( Skipping..." << std::endl;
                return false;
            }
            _stat_decrease_counters.at(5)++;

            _run = ev_mctrack->run();
            _subrun = ev_mctrack->subrun();
            _eventid = ev_mctrack->event_id();


            larlite::mctrack the_mctrack;
            // Find the mctrack that is a muon starting from the neutrino interaction
            // We also require it is 1m long at least... might want to remove this if you are using this
            // filter for some reason other than a multiple coloumb scattering analysis which requires 1m
            size_t n_found_MCTs = 0;
            for (auto const& mct : *ev_mctrack) {

                // Sometimes mctracks have zero size. No idea why. Skip them.
                if ( mct.size() < 3 ) continue;

                // origin == 1 means comes from neutrino interaction (IE not cosmic)
                if (mct.Origin() != 1 ) continue;

                // MCTrack has to be truly a muon
                if (mct.PdgCode() != 13 ) continue;

                //MCTrack should start VERY CLOSE to nu vtx
                // (we're talking like 10^-28 or smaller ... below floating point precision)
                if ( (mct.front().Position().Vect() - nu_vtx).Mag2() > 0.0001) continue;

                // Enforce the muon is exits the fiducial volume.
                if ( _fidvol.Contain(mct.back().Position().Vect()) ) continue;

                // Require the muon is at least one meter in length.
                if ((mct.back().Position().Vect() - mct.front().Position().Vect()).Mag() < 100.) continue;

                // Require the muon exits the actual TPC
                // this is done by requiring the energy of the last step in the MCTrack std::vector
                // (which corresponds to the last energy deposition point inside of the TPC)
                // has more energy than the rest mass energy (105.658 MeV)
                if ((mct.back().E() < 106.)) continue;

                // Save this mctrack!
                the_mctrack = mct;
                // At this point you have found a viable muon!
                n_found_MCTs++;

            }

            // If you didn't find any viable muons, skip the event
            if (!n_found_MCTs) return false;
            _stat_decrease_counters.at(6)++;

            // If somehow you found more than one muon, something has gone wrong. Skip the event.
            if (n_found_MCTs > 1) {
                print(larlite::msg::kWARNING, __FUNCTION__, Form("More than one viable muon in the event?!"));
                return false;
            }
            _stat_decrease_counters.at(7)++;


            // Now "the_mctrack" is the one MCTrack we care about. It should have length > 1 and be
            // the muon from a true numuCC interaction
            auto const mct_start = the_mctrack.front().Position().Vect();
            auto const mct_end   = the_mctrack.back().Position().Vect();

            _MCT_PDG = the_mctrack.PdgCode();
            _MCT_origin = the_mctrack.Origin();

            // Analyze the MCTrack!
            flip = false;
            _module->AnalyzeTrack(the_mctrack, _run, _subrun, _eventid, flip);

            ///////////////////////////////////////////
            // END kMCBNBRecoTrackExiting case
            ///////////////////////////////////////////
        }


        ////////////////////////////////////////////////////////////////////////////////////
        /// Actual data BNB selected events. One reco track is saved per event to the input file
        /// Also one reconstructed vertex is saved to the input file (to determine track direction)
        /// Any hand-scanning information is merged in in later analysis, not in this code.
        ////////////////////////////////////////////////////////////////////////////////////
        else if (_ana_type == MCSBiasStudy::kDataBNBSelectedRecoTrack)
        {
            _stat_decrease_counters.at(0)++;
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
            _stat_decrease_counters.at(1)++;

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
            _stat_decrease_counters.at(2)++;

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

        std::cout << "Printing out _stat_decrease_counters vector:" << std::endl;
        for (size_t i = 0; i < _stat_decrease_counters.size(); ++i)
            std::cout << " _stat_decrease_counters.at(" << i << ") = " << _stat_decrease_counters.at(i) << std::endl;

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
