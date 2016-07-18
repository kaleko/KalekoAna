#ifndef LARLITE_EVENTSELECTOR_CXX
#define LARLITE_EVENTSELECTOR_CXX

#include "EventSelector.h"
#include "DataFormat/vertex.h"
#include "DataFormat/opflash.h"

namespace larlite {

    bool EventSelector::initialize() {

        total_events = 0;
        passed_events = 0;

        _myGeoAABox = FidVolBox();


        return true;
    }

    bool EventSelector::analyze(storage_manager* storage) {

        total_events++;

        _flash_in_bgw = false;
        _is_areco_vtx_in_fidvol = false;
        _is_atrack_fromvtx = false;
        _longest_track_nearflash_z = false;
        _longest_trk_contained = false;
        _longest_trk_range_longenough = false;
        _passes_all_cuts = false;


        auto ev_vtx = storage->get_data<event_vertex>("pandoraNu");
        if (!ev_vtx) {
            print(larlite::msg::kERROR, __FUNCTION__, Form("Did not find specified data product, vertex!"));
            return false;
        }
        if (!ev_vtx->size()) {
            print(larlite::msg::kERROR, __FUNCTION__, Form("Zero reconstructed vertices in this event!"));
            return false;
        }

        auto ev_track = storage->get_data<event_track>("pandoraNu");
        if (!ev_track) {
            print(larlite::msg::kERROR, __FUNCTION__, Form("Did not find specified data product, track!"));
            return false;
        }
        if (!ev_track->size()) {
            print(larlite::msg::kERROR, __FUNCTION__, Form("Zero reconstructed tracks in this event!"));
            return false;
        }
        auto ev_opflash = storage->get_data<event_opflash>("opflashSat");
        if (!ev_opflash) {
            print(larlite::msg::kERROR, __FUNCTION__, Form("Did not find specified data product, opflash!"));
            return false;
        }
        if (!ev_opflash->size()) {
            print(larlite::msg::kERROR, __FUNCTION__, Form("opflash size is zero!"));
            return false;
        }

        // Loop over vertices.
        // For each vertex in fiducial volume, loop over reco tracks
        // If you find a reco track starting w/in 5cm from any vertex, store it
        // then if you find a longer reco track, store that one instead
        track longest_track;
        double longest_trk_len = 0;
        for (auto const& vtx : *ev_vtx) {
            ::geoalgo::Vector vertex = ::geoalgo::Vector(vtx.X(), vtx.Y(), vtx.Z());

            if (!_myGeoAABox.Contain(vertex)) continue;

            _is_areco_vtx_in_fidvol = true;

            for (auto const& trk : *ev_track) {
                if (!trk.NumberTrajectoryPoints()) continue;
                ::geoalgo::Vector trkstart = ::geoalgo::Vector(trk.Vertex());
                ::geoalgo::Vector trkend   = ::geoalgo::Vector(trk.End() );
                double dist_trk_start_vtx = vertex.Dist(trkstart);
                double dist_trk_end_vtx = vertex.Dist(trkend);

                if (dist_trk_start_vtx > 5. && dist_trk_end_vtx > 5.) continue;

                _is_atrack_fromvtx = true;
                double trk_len = trkstart.Dist(trkend);
                if (trk_len > longest_trk_len) {
                    longest_track = trk;
                    longest_trk_len = trk_len;
                }
            }
        } //done looping over vertices and finding the longest track

        if (longest_trk_len > 0) {
            _longest_trk_contained = _myGeoAABox.Contain(::geoalgo::Vector(longest_track.Vertex())) &&
                                     _myGeoAABox.Contain(::geoalgo::Vector(longest_track.End()));

            _longest_trk_range_longenough = longest_trk_len > 75.;


            opflash theflash;
            for (auto const& flash : *ev_opflash) {

                if (flash.Time() > 3.55 && flash.Time() < 5.15 && flash.TotalPE() > 50.) {
                    _flash_in_bgw = true;
                    // Keep track of the brightest flash in the BGW
                    if (flash.TotalPE() > theflash.TotalPE()) {
                        theflash = flash;
                        _longest_track_nearflash_z =
                            flashDistZ(longest_track, theflash.ZCenter()) < 80.;//10.;

                    }
                }
            }



        }

        _passes_all_cuts = _is_areco_vtx_in_fidvol &&
                           _is_atrack_fromvtx &&
                           _longest_track_nearflash_z &&
                           _longest_trk_contained &&
                           _longest_trk_range_longenough;

        if (_passes_all_cuts) passed_events++;

        return _passes_all_cuts;
    }

    bool EventSelector::finalize() {

        std::cout << "EventSelector finalize! "
                  << total_events << " analyzed in total, "
                  << passed_events << " passed the filter." << std::endl;

        return true;
    }

    double EventSelector::flashDistZ(const track &longest_track, const double flash_z) {

        // If flash is between start/end of track in Z (including flipping track direction)
        // return true

        /// I'm jsut going to copy christoph's exactly to try to duplicate
        // double FlashTrackDist(double flash, double start, double end) {
        double end = longest_track.End().Z();
        double start = longest_track.Vertex().Z();
        double flash = flash_z;
        if (end >= start) {
            if (flash < end && flash > start) return 0;
            else return TMath::Min(fabs(flash - start), fabs(flash - end));
        }
        else {
            if (flash > end && flash < start) return 0;
            else return TMath::Min(fabs(flash - start), fabs(flash - end));
        }
    }
}
#endif
