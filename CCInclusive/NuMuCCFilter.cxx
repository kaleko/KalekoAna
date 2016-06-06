#ifndef LARLITE_NUMUCCFILTER_CXX
#define LARLITE_NUMUCCFILTER_CXX

#include "NuMuCCFilter.h"
#include "DataFormat/mctruth.h"

namespace larlite {

    bool NuMuCCFilter::initialize() {


        total_events = 0;
        kept_events = 0;

        // fidvol_dist = 5.;//10.;
        // fidvol_dist_y = 5.;//20.;

        // //Box here is TPC
        // _fidvolBox.Min( 0 + fidvol_dist,
        //                 -(::larutil::Geometry::GetME()->DetHalfHeight()) + fidvol_dist_y,
        //                 0 + fidvol_dist);

        // _fidvolBox.Max( 2 * (::larutil::Geometry::GetME()->DetHalfWidth()) - fidvol_dist,
        //                 ::larutil::Geometry::GetME()->DetHalfHeight() - fidvol_dist_y,
        //                 ::larutil::Geometry::GetME()->DetLength() - fidvol_dist);

        fidvol_dist_x = 20.;
        fidvol_dist_y = 20.;
        fidvol_dist_z = 10.;

        //Box here is TPC
        _fidvolBox.Min( 0 + fidvol_dist_x,
                        -(::larutil::Geometry::GetME()->DetHalfHeight()) + fidvol_dist_y,
                        0 + fidvol_dist_z);

        _fidvolBox.Max( 2 * (::larutil::Geometry::GetME()->DetHalfWidth()) - fidvol_dist_x,
                        ::larutil::Geometry::GetME()->DetHalfHeight() - fidvol_dist_y,
                        ::larutil::Geometry::GetME()->DetLength() - fidvol_dist_z);

        return true;
    }

    bool NuMuCCFilter::analyze(storage_manager* storage) {

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

        // Require interactino is CC and neutrino pdg is numu
        if (ev_mctruth->at(0).GetNeutrino().CCNC() || ev_mctruth->at(0).GetNeutrino().Nu().PdgCode() != 14)
            return false;

        kept_events++;

        return true;
    }

    bool NuMuCCFilter::finalize() {

        std::cout << "NuMuCCFilter: Total events = " << total_events << std::endl;
        std::cout << "NuMuCCFilter: Final kept events = " << kept_events << std::endl;

        return true;
    }

}
#endif
