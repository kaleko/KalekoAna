#ifndef LARLITE_KALEKOTRACKSTITCHER_CXX
#define LARLITE_KALEKOTRACKSTITCHER_CXX

#include "KalekoTrackStitcher.h"

namespace larlite {

    bool KalekoTrackStitcher::initialize() {

        if ( _base_producer.empty() ||  _match_producer.empty() ) {
            print(larlite::msg::kERROR, __FUNCTION__, Form("Base Producer or Match Producer not set!"));
            return false;
        }

        if (!_debug_tree) {
            _debug_tree = new TTree("debug_tree", "debug_tree");
            _debug_tree->Branch("dotprod", &_dotprod, "dotprod/D");
            _debug_tree->Branch("infdist", &_infdist, "infdist/D");
            _debug_tree->Branch("startdist", &_startdist, "startdist/D");
            _debug_tree->Branch("enddist", &_enddist, "enddist/D");
        }

        print(larlite::msg::kDEBUG, __FUNCTION__,
              Form("Initializing with base producer: %s, match producer: %s.",
                   _base_producer.c_str(),
                   _match_producer.c_str()));

        _geoalg = geoalgo::GeoAlgo();

        return true;
    }

    bool KalekoTrackStitcher::analyze(storage_manager* storage) {

        // Read in base track producer
        auto ev_base_track = storage->get_data<event_track>(_base_producer);
        if (!ev_base_track) {
            print(larlite::msg::kERROR, __FUNCTION__,
                  Form("Did not find specified data product, track by %s!", _base_producer.c_str()));
            return false;
        }
        if (!ev_base_track->size()) {
            //print(larlite::msg::kERROR, __FUNCTION__, Form("Zero reconstructed tracks in this event!"));
            return false;
        }

        // Read in match track producer
        auto ev_match_track = storage->get_data<event_track>(_match_producer);
        if (!ev_match_track) {
            print(larlite::msg::kERROR, __FUNCTION__,
                  Form("Did not find specified data product, track by %s!", _match_producer.c_str()));
            return false;
        }
        if (!ev_match_track->size()) {
            //print(larlite::msg::kERROR, __FUNCTION__, Form("Zero reconstructed tracks in this event!"));
            return false;
        }


        // Loop over base tracks
        for (auto const& base_trk : *ev_base_track) {

            if (!trackViable(base_trk)) continue;

            // Loop over potential match tracks
            for (auto const& match_trk : *ev_match_track) {

                if (!trackViable(match_trk)) continue;

                bool matched = tracksMatched(base_trk, match_trk);
                //if (matched) std::cout << "Match found!" << std::endl;
            }
        }
        return true;
    }

    bool KalekoTrackStitcher::finalize() {

        if (_fout) {
            if (_debug_tree)
                _debug_tree->Write();
        }

        return true;
    }

    bool KalekoTrackStitcher::trackViable(const larlite::track &trk) {

        if ( trk.NumberTrajectoryPoints() < 2 ) return false;

        double _trk_len = (trk.End() - trk.Vertex()).Mag();
        if ( _trk_len < 3.) return false;

        return true;
    }

    bool KalekoTrackStitcher::tracksMatched(const larlite::track &base_trk,
                                            const larlite::track &match_candidate ) {


        double match_start_base_start_dist = (match_candidate.Vertex() - base_trk.Vertex()).Mag2();
        double match_end_base_start_dist = (match_candidate.End() - base_trk.Vertex()).Mag2();
        double match_start_base_end_dist = (match_candidate.Vertex() - base_trk.End()).Mag2();
        double match_end_base_end_dist = (match_candidate.End() - base_trk.End()).Mag2();
        bool flip_match = match_start_base_start_dist < match_end_base_start_dist ? false : true;

        TVector3 base_dir = base_trk.End() - base_trk.Vertex();
        TVector3 match_dir = match_candidate.End() - match_candidate.Vertex();
        if (flip_match) match_dir *= -1.;

        _dotprod = base_dir.Unit().Dot(match_dir.Unit());

        _startdist = flip_match ? match_end_base_start_dist : match_start_base_start_dist;
        _startdist = std::sqrt(_startdist);
        _enddist = flip_match ? match_end_base_end_dist : match_start_base_end_dist;
        _enddist = std::sqrt(_enddist);

        /// Create a GeoLine for the base track...
        ///   >> one point on the line is track start
        ///   >> one point on the line is track end
        ///   >> line extends infinitely in both directions
        ::geoalgo::Line base_line(
            ::geoalgo::Point_t(base_trk.Vertex()),
            ::geoalgo::Point_t(base_trk.End()));

        /// Ask how far the matched start point is from the line
        ::geoalgo::Point_t match_vtx = flip_match ?
                                       ::geoalgo::Point_t(match_candidate.End()) :
                                       ::geoalgo::Point_t(match_candidate.Vertex());

        _infdist = _geoalg.SqDist(base_line, match_vtx);
        _infdist = std::sqrt(_infdist);
        _debug_tree->Fill();


        return true;
    }

    larlite::track KalekoTrackStitcher::stitchTracks(const larlite::track &base_trk,
            const larlite::track &match_candidate ) {
        larlite::track result;
        return result;
    }

}
#endif
