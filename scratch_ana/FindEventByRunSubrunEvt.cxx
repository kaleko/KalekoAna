#ifndef LARLITE_FINDEVENTBYRUNSUBRUNEVT_CXX
#define LARLITE_FINDEVENTBYRUNSUBRUNEVT_CXX

#include "FindEventByRunSubrunEvt.h"

namespace larlite {

    bool FindEventByRunSubrunEvt::initialize() {

        if (!runno || !subrunno || !evtno) {
            print(larlite::msg::kERROR, __FUNCTION__, Form("Didn't set run number, subrun number, or event number!"));
            return false;
        }

        _found_event = false;
        _file_index = 0;

        return true;
    }

    bool FindEventByRunSubrunEvt::analyze(storage_manager* storage) {

        // If you've already found the event, skip it
        if ( _found_event ) return false;

        _file_index++;
        if ( storage->run_id() == runno && storage->subrun_id() == subrunno && storage->event_id() == evtno ) {
            std::cout << "EVENT FOUND!" << std::endl;
            _found_event = true;
            return true;
        }
        else return false;

    }

    bool FindEventByRunSubrunEvt::finalize() {

        std::cout << "Run " << runno
                  << ", subrun " << subrunno
                  << ", event " << evtno
                  << " found! It's located as index " << _file_index
                  << " in your input file(s)." << std::endl;
        std::cout << "If you're using the larlite EVD, the counters are off by 1, so you want to enter in " << _file_index - 1
                  << " into the \"Go to:\" box." << std::endl;

        return true;
    }

}
#endif
