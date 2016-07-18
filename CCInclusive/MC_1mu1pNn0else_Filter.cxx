#ifndef LARLITE_MC_1MU1PNN0ELSE_FILTER_CXX
#define LARLITE_MC_1MU1PNN0ELSE_FILTER_CXX

#include "MC_1mu1pNn0else_Filter.h"
#include "DataFormat/mctruth.h"

namespace larlite {

    bool MC_1mu1pNn0else_Filter::initialize() {

        total_events = 0;
        kept_events = 0;

        _fidvolBox = FidVolBox();
        return true;
    }

    bool MC_1mu1pNn0else_Filter::analyze(storage_manager* storage) {

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

        // Require the neutrino interaction is inside the fiducial volume
        if (!_fidvolBox.Contain(ev_mctruth->at(0).GetNeutrino().Nu().Trajectory().front().Position()))
            return false;

        //Enforce that there is exactly 1 electron, above 20MeV kinetic energy
        //Don't care about neutrons, protons.
        //Don't care about other particles if they are below 20MeV KE
        //I don't care about neutrons, weird quarks, the neutrino, etc.
        size_t n_muons = 0;
        size_t n_protons = 0;

        for (auto const& particle : ev_mctruth->at(0).GetParticles()) {

            // Only particles with status code 1 are relevant
            if ( particle.StatusCode() != 1 ) continue;

            // //Note: this KE is in units of GEV!
            // double KE = particle.Trajectory().at(0).E() - particle.Mass();

            // //Don't care about any particles with less than 20 MeV KE
            // if ( KE < 0.02 ) continue;

            //Count up the number of muons
            if ( particle.PdgCode() == 13 )
                n_muons++;

            //Count up the number of protons
            if ( particle.PdgCode() == 2212 )
                n_protons++;

            //Skip this event if any other particles exist... this filter
            //is to pick out 1mu1pNn events!
            // 22  => gammas
            // 211 => charged pions
            // 111 => pi0s
            // 11  => electrons
            // 321 => kaons
            if ( abs(particle.PdgCode()) == 22  ||
                    abs(particle.PdgCode()) == 211 ||
                    abs(particle.PdgCode()) == 111 ||
                    abs(particle.PdgCode()) == 11  ||
                    abs(particle.PdgCode()) == 321 ||
                    abs(particle.PdgCode()) == 2112 ||
                    //abs(particle.PdgCode()) == 2000000101 || //kbindino
                    abs(particle.PdgCode()) == 14 ||
                    abs(particle.PdgCode()) == 12 ||
                    particle.PdgCode() == -13 ||
                    //particle.PdgCode() == 1000180400 || //argon 40
                    abs(particle.PdgCode()) == 311 ||
                    abs(particle.PdgCode()) == 3222 || //strange baryons
                    abs(particle.PdgCode()) == 3122 ||
                    abs(particle.PdgCode()) == 1000020040 ) //deuteron
                return false;

        }

        //Skip this event if the number of muons is not 1 or # protons is not 1
        if ( n_muons != 1 || n_protons != 1 )
            return false;

        kept_events++;
        return true;

    }

    bool MC_1mu1pNn0else_Filter::finalize() {

        std::cout << "MC_1mu1pNn0else_Filter: Total events = " << total_events << std::endl;
        std::cout << "MC_1mu1pNn0else_Filter: Final kept events = " << kept_events << std::endl;

        return true;
    }

}
#endif
