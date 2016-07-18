#ifndef LARLITE_MC_1MU1PNN0ELSE_FILTER_MCTracks_MCTRACKS_CXX
#define LARLITE_MC_1MU1PNN0ELSE_FILTER_MCTracks_MCTRACKS_CXX

#include "MC_1mu1pNn0else_Filter_MCTracks.h"
#include "DataFormat/mctrack.h"
#include "DataFormat/mcshower.h"
#include "DataFormat/mctruth.h"

namespace larlite {

    bool MC_1mu1pNn0else_Filter_MCTracks::initialize() {

        total_events = 0;
        kept_events = 0;
        kaleko_ctr = 0;

        _fidvolBox = FidVolBox();
        
        return true;
    }

    bool MC_1mu1pNn0else_Filter_MCTracks::analyze(storage_manager* storage) {


        // std::cout << "Filter analyzing event index: " << kaleko_ctr << std::endl;
        kaleko_ctr++;
        //Grab the mctracks
        auto ev_mctrack = storage->get_data<event_mctrack>("mcreco");
        if (!ev_mctrack) {
            print(larlite::msg::kERROR, __FUNCTION__, Form("Did not find specified data product, mctrack!"));
            return false;
        }
        //Grab the mcshowers
        auto ev_mcshower = storage->get_data<event_mcshower>("mcreco");
        if (!ev_mcshower) {
            print(larlite::msg::kERROR, __FUNCTION__, Form("Did not find specified data product, mcshower!"));
            return false;
        }
        //Grab the MCTruth
        auto ev_mctruth = storage->get_data<event_mctruth>("generator");
        if (!ev_mctruth) {
            print(larlite::msg::kERROR, __FUNCTION__, Form("Did not find specified data product, mctruth!"));
            return false;
        }

        total_events++;

        // Require exactly one neutrino interaction
        if (ev_mctruth->size() != 1) {
            print(larlite::msg::kINFO, __FUNCTION__, Form("ev_mctruth size is not 1!"));
            return false;
        }

        auto const& nu_vtx = ev_mctruth->at(0).GetNeutrino().Nu().Trajectory().front().Position().Vect();
        // Require the neutrino interaction is inside the fiducial volume
        if (!_fidvolBox.Contain(nu_vtx))
            return false;

        // Make a small 1cm geosphere around neutrino vertex
        _nu_sphere.Center(nu_vtx);
        _nu_sphere.Radius(1.);
        // std::cout << "FILTER: here is the neutrino vertex " << _nu_sphere.Center() << std::endl;
        // Require no MCShowers have truth start point inside of this sphere
        for (auto const& mcs : *ev_mcshower)
            if ( _nu_sphere.Contain(mcs.Start().Position().Vect()) ) return false;


        // Require exactly two MCTracks originating from the neutrino vertex (one muon, one proton)
        // Store the muon and the proton mctracks for later investigation
        bool _found_muon = false;
        mctrack mu_mctrack;
        bool _found_proton = false;
        mctrack p_mctrack;

        for (auto const& mct : *ev_mctrack) {

            if (!mct.size()) continue;
            if (!_nu_sphere.Contain(mct.front().Position().Vect())) continue;

            if (mct.PdgCode() == 13 ) {
                //if you've already found a muon mctrack, return false
                if (_found_muon) return false;
                _found_muon = true;
                mu_mctrack = mct;
            }
            else if (mct.PdgCode() == 2212) {
                // if you've already found a proton mctrack, return false;
                if (_found_proton) return false;
                _found_proton = true;
                p_mctrack = mct;
            }

            else {
                // if (abs(mct.PdgCode()) == 211) {
                //     std::cout << "FILTER Found a pion. Here is the neutrino vertex: "
                //               << _nu_sphere.Center() << std::endl;
                //     std::cout << "FILTER: Here is the pion start point: ";
                //     mct.front().Position().Vect().Print();
                //     std::cout << "FILTER: Here is its deposited energy: " << mct.front().E() - mct.back().E() << std::endl;
                // }
                // print(larlite::msg::kWARNING, "MC_1mu1pNn0else_Filter_MCTracks",
                //       Form("Found a mctrack that wasn't p or mu! Its pdg is %d.", mct.PdgCode()));
                return false;
            }
        }


        // If you haven't found a proton and a muon, skip this event.
        if (!_found_proton || !_found_muon) return false;

        // Now you should have exactly 1 muon mctrack and 1 proton mctrack originating from close
        // to the true neutrino vertex, and you should have zero mcshowers originating from close
        // to the true neutrino vertex

        kept_events++;
        return true;

    }

    bool MC_1mu1pNn0else_Filter_MCTracks::finalize() {

        std::cout << "MC_1mu1pNn0else_Filter_MCTracks: Total events = " << total_events << std::endl;
        std::cout << "MC_1mu1pNn0else_Filter_MCTracks: Final kept events = " << kept_events << std::endl;
        std::cout << "NOTE! These events could have lots of protons, just some are below mctrack threshold." << std::endl;
        std::cout << " That is to say, this filter is based on reconstructable things ... events that "
                  << "LOOK like there is 1 mu and 1 p and 0 else, from the detector point of view!" << std::endl;
        std::cout << "Also I've found events where a neutrino makes a 1GeV proton which travels " << std::endl;
        std::cout << " like 25 cm, then 'protoninelastic' happens and a 700 MeV charged pion comes out. " << std::endl;
        std::cout << " Such events PASS this filter." << std::endl;

        return true;
    }

}
#endif
