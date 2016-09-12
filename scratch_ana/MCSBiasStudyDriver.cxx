#ifndef LARLITE_MCSBIASSTUDYDRIVER_CXX
#define LARLITE_MCSBIASSTUDYDRIVER_CXX

#include "MCSBiasStudyDriver.h"

namespace larlite {

    bool MCSBiasStudyDriver::initialize() {

        if (!_module)
            _module = new MCSBiasStudy();

        return true;
    }

    bool MCSBiasStudyDriver::analyze(storage_manager* storage) {

        auto ev_track = storage->get_data<event_track>("pandoraNuPMA");
        if (!ev_track) {
            print(larlite::msg::kERROR, __FUNCTION__, Form("Did not find specified data product, track!"));
            return false;
        }
        if (!ev_track->size())
            return false;

        for (auto const &track : *ev_track)
            _module->AnalyzeTrack(track);

        return true;
    }

    bool MCSBiasStudyDriver::finalize() {

        if (_fout) {
            _fout->cd();
            _module->GetTree()->Write();
        }

        return true;
    }

}
#endif
