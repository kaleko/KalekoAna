#ifndef LARLITE_NUMUCCFILTER_CXX
#define LARLITE_NUMUCCFILTER_CXX

#include "NuMuCCFilter.h"
#include "DataFormat/mctruth.h"
#include "DataFormat/mcflux.h"

namespace larlite {

    bool NuMuCCFilter::initialize() {


        total_events = 0;
        kept_events = 0;
        _fidvolBox = FidVolBox();

        return true;
    }

    bool NuMuCCFilter::analyze(storage_manager* storage) {

        bool ret = false;

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
        // AND interaction is CC
        // AND pdg is numu (14)
        // AND neutrino energy is above cutoff threshold (default threshold is 0)
        if (_fidvolBox.Contain(ev_mctruth->at(0).GetNeutrino().Nu().Trajectory().front().Position())
                && !ev_mctruth->at(0).GetNeutrino().CCNC()
                && ev_mctruth->at(0).GetNeutrino().Nu().PdgCode() == 14
                && ev_mctruth->at(0).GetNeutrino().Nu().Trajectory().front().E() > _min_nu_energy)
            ret = true;

        // If this filter toggled to only look at numuCC from kaons, read in the mcflux to determine neutrino ancestry
        // and throw out this event if neutrino doesn't have fndecay == 5 (K+ --> numu + mu+)
        if (_keep_only_numu_from_kaons) {

            auto ev_mcflux = storage->get_data<event_mcflux>("generator");
            if (!ev_mcflux) {
                print(larlite::msg::kERROR, __FUNCTION__, Form("Did not find specified data product, mcflux!"));
                return false;
            }
            // Require exactly one neutrino interaction
            if (ev_mcflux->size() != 1) {
                print(larlite::msg::kINFO, __FUNCTION__, Form("ev_mcflux size is not 1!"));
                return false;
            }

            if ( ev_mcflux->at(0).fndecay >= 10  ) ret = true;

        }

        if ( (!_flip && ret) || (_flip && !ret) ) kept_events++;

        return _flip ? !ret : ret;
    }

    bool NuMuCCFilter::finalize() {
        std::cout << "NuMuCCFilter has _keep_only_numu_from_kaons set to " << _keep_only_numu_from_kaons << "." << std::endl;
        std::cout << "NuMuCCFilter has _flip set to " << _flip << "." << std::endl;
        std::cout << "NuMuCCFilter: Total events = " << total_events << std::endl;
        std::cout << "NuMuCCFilter: Final kept events = " << kept_events << std::endl;

        return true;
    }

}
#endif
