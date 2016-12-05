#ifndef LARLITE_MC_CCNUMU_MUONMCTRACKCONTAINEDFILTER_CXX
#define LARLITE_MC_CCNUMU_MUONMCTRACKCONTAINEDFILTER_CXX

#include "MC_CCnumu_MuonMCTrackContainedFilter.h"

namespace larlite {

    bool MC_CCnumu_MuonMCTrackContainedFilter::initialize() {

        _n_total_events = 0;
        _n_kept_events = 0;

        for (size_t i = 0; i < 10; ++i)
            _counters.push_back(0);

        return true;
    }

    bool MC_CCnumu_MuonMCTrackContainedFilter::analyze(storage_manager* storage) {

        _n_total_events++;
        _counters[0]++;

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
        _counters[1]++;

        auto const &mctruth = ev_mctruth->at(0);

        // Make sure the event is numuCC inside of the fiducial volume


        //Enforce CC interaction channel
        if ( mctruth.GetNeutrino().CCNC() != 0 ) return false;
        _counters[2]++;

        // If neutrino interacts outside of fiducial volume, skip event
        auto const nu_vtx = mctruth.GetNeutrino().Nu().Trajectory().back().Position().Vect();
        if (!_fidvol.Contain(nu_vtx)) return false;
        _counters[3]++;

        // If neutrino is not a numu, skip event
        if (mctruth.GetNeutrino().Nu().PdgCode() != 14) return false;
        _counters[4]++;


        // Now we make sure the muon mctrack from the interaction is fully contained in the fiducial volume

        // Grab the MCTracks
        auto ev_mctrack = storage->get_data<event_mctrack>("mcreco");
        if (!ev_mctrack) {
            print(larlite::msg::kERROR, __FUNCTION__, Form("Did not find specified data product, mctrack!"));
            return false;
        }
        if (!ev_mctrack->size()) {
            print(larlite::msg::kERROR, __FUNCTION__, Form("Zero mctracks in this event?! Skipping!"));
            return false;
        }


        // Find the mctrack that is a muon starting from the neutrino interaction
        // We also require it is 1m long at least... might want to remove this if you are using this
        // filter for some reason other than a multiple coloumb scattering analysis which requires 1m
        size_t n_found_muons = 0;
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
            if (!_fidvol.Contain(mct.back().Position().Vect())) continue;

            // Require the muon is at least one meter in length.
            if ((mct.back().Position().Vect() - mct.front().Position().Vect()).Mag() < 100.) continue;

            // At this point you have found a viable muon!
            n_found_muons++;

        }

        // If you didn't find any viable muons, skip the event
        if (!n_found_muons) return false;
        _counters[5]++;

        // If somehow you found more than one muon, something has gone wrong. Skip the event.
        if (n_found_muons > 1) {
            print(larlite::msg::kWARNING, __FUNCTION__, Form("More than one viable muon in the event?!"));
            return false;
        }

        _n_kept_events++;

        return true;
    }

    bool MC_CCnumu_MuonMCTrackContainedFilter::finalize() {

        std::cout << _n_total_events << " total events analyzed, "
                  << _n_kept_events << " events passed MC_CCnumu_MuonMCTrackContainedFilter." << std::endl;

        std::cout << "Printing out the _counters vector:"<<std::endl;
        for (size_t i = 0; i < _counters.size(); ++i)
            std::cout << "_counters[" << i << "] = " << _counters[i] << std::endl;

        return true;
    }

}
#endif
