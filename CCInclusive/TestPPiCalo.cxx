#ifndef LARLITE_TESTPPICALO_CXX
#define LARLITE_TESTPPICALO_CXX

#include "TestPPiCalo.h"
#include "DataFormat/mctrack.h"
#include "DataFormat/track.h"
#include "DataFormat/calorimetry.h"

namespace larlite {

    bool TestPPiCalo::initialize() {

        _h_avgdedx_proton = new TH1F("h_avgdedx_proton", "Average dEdx In Most-Hits Plane: Single Protons", 250, 0, 15);
        _h_avgdedx_pion = new TH1F("h_avgdedx_pion", "Average dEdx In Most-Hits Plane: Single Pions", 250, 0, 15);
        _h_avgdedx_muon = new TH1F("h_avgdedx_muon", "Average dEdx In Most-Hits Plane: Single Muons", 250, 0, 15);

        n_proton_evts_tot = 0;
        n_pion_evts_tot = 0;
        n_muon_evts_tot = 0;

        return true;
    }

    bool TestPPiCalo::analyze(storage_manager* storage) {

        //Grab the MCTrack
        auto ev_mctrack = storage->get_data<event_mctrack>("mcreco");
        if (!ev_mctrack) {
            print(larlite::msg::kERROR, __FUNCTION__, Form("Did not find specified data product, mctrack!"));
            return false;
        }
        // For some reason pions seem to have multiple MCTracks... for now I'll just use the longest mctrack.

        // if (ev_mctrack->size() != 1) {
        //     print(larlite::msg::kERROR, __FUNCTION__,
        //           Form("More than 1 MCTrack in this event! There are %zu!", ev_mctrack->size()));
        //     return false;
        // }

        // Pick the longest MCTrack in the event
        mctrack themctrack;
        double longest_mctrack_len = 0;
        for (auto const& mct : *ev_mctrack) {
            if (!mct.size()) continue;
            double this_mctrack_len = (mct.front().Position().Vect() - mct.back().Position().Vect()).Mag();
            if ( this_mctrack_len > longest_mctrack_len ) {
                longest_mctrack_len = this_mctrack_len;
                themctrack = mct;
            }
        }
        if (!longest_mctrack_len) {
            print(larlite::msg::kERROR, __FUNCTION__, Form("No MCTracks with nonzero size in this event!"));
            return false;
        }

        int mctrack_pdg = themctrack.PdgCode();
        int proton_1_pion_2_muon_3 = 0;
        if ( abs(mctrack_pdg) == 2212 ) {
            proton_1_pion_2_muon_3 = 1;
            n_proton_evts_tot++;
        }
        else if (abs(mctrack_pdg) == 211) {
            proton_1_pion_2_muon_3 = 2;
            n_pion_evts_tot++;
        }
        else if ( abs(mctrack_pdg) == 13) {
            proton_1_pion_2_muon_3 = 3;
            n_muon_evts_tot++;
        }
        else {
            print(larlite::msg::kERROR, __FUNCTION__, Form("The MCTrack was not proton or pion or muon. What the heck?"));
            return false;
        }

        auto ev_track = storage->get_data<event_track>("pandoraNuPMA");
        if (!ev_track) {
            print(larlite::msg::kERROR, __FUNCTION__, Form("Did not find specified data product, track!"));
            return false;
        }
        if ( !ev_track->size() ) {
            print(larlite::msg::kWARNING, __FUNCTION__,
                  Form("# of reconstructed tracks in this event is zero!"));
            return false;
        }

        /// Choose the reconstructed track with start or end point within 1cm of mctrack start or end point
        /// If there are none, throw out this event.
        track thetrack;
        int thetrack_idx = -1;
        bool found_good_recotrack = false;
        for (auto const &trk : *ev_track) {
            thetrack_idx++;
            if ( (trk.Vertex() - themctrack.front().Position().Vect()).Mag() < 1. ||
                    (trk.End() - themctrack.front().Position().Vect()).Mag() < 1. ) {
                if (found_good_recotrack) {
                    print(larlite::msg::kWARNING, __FUNCTION__,
                          Form("Found more than one good reco track matchint the mctrack ... throwing out event."));
                    return false;
                }
                else
                    found_good_recotrack = true;
            } // End if found good track
        }//End loop over reco tracks

        if (!found_good_recotrack) {
            print(larlite::msg::kWARNING, __FUNCTION__,
                  Form("Did not find good reco track to match mctrack in this event."));
            return false;
        }
        event_calorimetry* ev_calo = nullptr;
        auto const& ass_calo_v = storage->find_one_ass(ev_track->id(),
                                 ev_calo,
                                 Form("%scalo", ev_track->name().c_str()));

        if (!ev_calo) {
            std::cerr << "no event_calorimetry!" << std::endl;
            return false;
        }

        if (ev_calo->empty()) {
            std::cout << "ev_calo empty" << std::endl;
            return false;
        }

        //size of ev calo is 3 (one per plane)... let's use the plane with the most hits
        // Choose the calo object for this track by the one
        // with the most number of hits in the dEdx vector
        // (this is how analysis tree does it)

        //ev_calo->at(ass_calo_v[long_track_idx][i]).dEdx().size()

        calorimetry thecalo;
        size_t tmp_nhits = 0;
        for (size_t i = 0; i < 3; ++i)
            if (ev_calo->at(ass_calo_v[thetrack_idx][i]).dEdx().size() > tmp_nhits) {
                thecalo = ev_calo->at(ass_calo_v[thetrack_idx][i]);
                tmp_nhits = thecalo.dEdx().size();
            }

        // Parameter of interest is average dedx over the whole track
        double dedx_avg = 0;
        for (size_t i = 0; i < thecalo.dEdx().size(); ++i)
            dedx_avg += thecalo.dEdx().at(i);
        dedx_avg /= thecalo.dEdx().size();

        if (proton_1_pion_2_muon_3 == 1) _h_avgdedx_proton->Fill(dedx_avg);
        if (proton_1_pion_2_muon_3 == 2) _h_avgdedx_pion->Fill(dedx_avg);
        if (proton_1_pion_2_muon_3 == 3) _h_avgdedx_muon->Fill(dedx_avg);

        return true;
    }

    bool TestPPiCalo::finalize() {

        std::cout << " TestPPiCalo Finalizing! " << std::endl;
        std::cout << " Of " << n_proton_evts_tot << " input proton events, "
                  << _h_avgdedx_proton->GetEntries() << " had one track, valid calo, etc. and ended up in a histogram."
                  << std::endl;
        std::cout << " Of " << n_pion_evts_tot << " input pion events, "
                  << _h_avgdedx_pion->GetEntries() << " had one track, valid calo, etc. and ended up in a histogram."
                  << std::endl;
        std::cout << " Of " << n_muon_evts_tot << " input muon events, "
                  << _h_avgdedx_muon->GetEntries() << " had one track, valid calo, etc. and ended up in a histogram."
                  << std::endl;

        if (_fout) {
            _h_avgdedx_proton->Write();
            _h_avgdedx_pion->Write();
            _h_avgdedx_muon->Write();
        }

        return true;
    }

}
#endif
