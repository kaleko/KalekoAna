#ifndef LARLITE_KALEKOTRACKSTITCHER_CXX
#define LARLITE_KALEKOTRACKSTITCHER_CXX

#include "KalekoTrackStitcher.h"

namespace larlite {

    bool KalekoTrackStitcher::initialize() {

        if ( _base_producer.empty() ||  _match_producer.empty() || _output_producer.empty() ) {
            print(larlite::msg::kERROR, __FUNCTION__, Form("Producers not properly set!"));
            return false;
        }

        if (!_debug_tree) {
            _debug_tree = new TTree("debug_tree", "debug_tree");
            _debug_tree->Branch("dotprod", &_dotprod, "dotprod/D");
            _debug_tree->Branch("infdist", &_infdist, "infdist/D");
            _debug_tree->Branch("startdist", &_startdist, "startdist/D");
            _debug_tree->Branch("enddist", &_enddist, "enddist/D");
            _debug_tree->Branch("endprojdist", &_endprojdist, "endprojdist/D");
            _debug_tree->Branch("startprojdist", &_startprojdist, "startprojdist/D");
            _debug_tree->Branch("baselen", &_baselen, "baselen/D");
        }

        print(larlite::msg::kDEBUG, __FUNCTION__,
              Form("Initializing with base producer: %s, match producer: %s.",
                   _base_producer.c_str(),
                   _match_producer.c_str()));

        _geoalg = geoalgo::GeoAlgo();

        return true;
    }

    bool KalekoTrackStitcher::analyze(storage_manager* storage) {

        idx_map.clear();

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

        // std::cout<<"# of base tracks in event is "<<ev_base_track->size()<<std::endl;
        // std::cout<<"# of match tracks in event is "<<ev_match_track->size()<<std::endl;

        // Loop over base tracks
        // for (auto const& base_trk : *ev_base_track) {
        for (size_t base_idx = 0; base_idx < ev_base_track->size(); ++base_idx) {
            auto const& base_trk = ev_base_track->at(base_idx);

            // Check if this base track has already been matched to any candidates:
            if ( idx_map.find(base_idx) != idx_map.end() ) {
                // This base track was already matched to another track... skip it
                continue;
            }

            if (!trackViable(base_trk)) continue;

            // Loop over potential match tracks
            // for (auto const& match_trk : *ev_match_track) {
            for (size_t match_idx = 0; match_idx < ev_match_track->size(); ++match_idx) {
                auto const& match_trk = ev_match_track->at(match_idx);

                // Check if this match candidate track has already been matched to a base
                bool already_matched = false;
                std::map<size_t, size_t>::const_iterator it;
                for (it = idx_map.begin(); it != idx_map.end(); ++it) {
                    if (it->second == match_idx) already_matched = true;
                }
                if (already_matched) {
                    // This candidate was already matched to a base... skip it
                    continue;
                }

                if (!trackViable(match_trk)) continue;

                bool matched = tracksMatched(base_trk, match_trk);

                if (matched) {
                    idx_map.insert(std::pair<size_t, size_t>(base_idx, match_idx));
                }
                //if (matched) std::cout << "Match found!" << std::endl;

            } // End loop over candidate matched tracks
        } // End loop over base tracks


        // Time to make my own track producer output.
        // Start by copying the base producer tracks, then stitch on potential additional tracks
        auto ev_stitched_track = storage->get_data<event_track>(_output_producer);

        // Set the event ID's and such correctly for my new tracks
        storage->set_id(ev_base_track->run(), ev_base_track->subrun(), ev_base_track->event_id());

        // Loop through base tracks, check for a potential match
        // If no match, just copy the base track into my new track producer
        // If match, stitch together and add to new track producer
        for (size_t base_idx = 0; base_idx < ev_base_track->size(); ++base_idx) {

            // Check if match was found for this base track
            std::map<size_t, size_t>::iterator it;
            it = idx_map.find(base_idx);
            // If so:
            if (it != idx_map.end()) {
                ev_stitched_track->push_back(
                    stitchTracks(ev_base_track->at(base_idx), ev_match_track->at(it->second))
                );
            }

            // If not
            else
                ev_stitched_track->push_back(ev_base_track->at(base_idx));


        }// End loop over base tracks again
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

        /// Create a GeoLine for the base track...
        ///   >> one point on the line is track start
        ///   >> one point on the line is track end
        ///   >> line extends infinitely in both directions
        ::geoalgo::Line base_line(
            ::geoalgo::Point_t(base_trk.Vertex()),
            ::geoalgo::Point_t(base_trk.End()));

        // Decide if the match track should be flipped
        ::geoalgo::Point_t projpoint_matchstart =
            _geoalg.ClosestPt(base_line, ::geoalgo::Point_t(match_candidate.Vertex()));
        ::geoalgo::Point_t projpoint_matchend =
            _geoalg.ClosestPt(base_line, ::geoalgo::Point_t(match_candidate.End()));
        double base_end_proj_matchstart_dist = ::geoalgo::Point_t(base_trk.End()).Dist(projpoint_matchstart);
        double base_end_proj_matchend_dist = ::geoalgo::Point_t(base_trk.End()).Dist(projpoint_matchend);
        bool flip_match = base_end_proj_matchstart_dist < base_end_proj_matchend_dist ? false : true;




        // double match_start_base_start_dist = (match_candidate.Vertex() - base_trk.Vertex()).Mag2();
        // double match_end_base_start_dist = (match_candidate.End() - base_trk.Vertex()).Mag2();
        // double match_start_base_end_dist = (match_candidate.Vertex() - base_trk.End()).Mag2();
        // double match_end_base_end_dist = (match_candidate.End() - base_trk.End()).Mag2();


        TVector3 base_dir = base_trk.End() - base_trk.Vertex();
        TVector3 match_dir = match_candidate.End() - match_candidate.Vertex();
        if (flip_match) match_dir *= -1.;

        _baselen = base_dir.Mag();

        _dotprod = base_dir.Unit().Dot(match_dir.Unit());

        // _startdist = flip_match ? match_end_base_start_dist : match_start_base_start_dist;
        // _startdist = std::sqrt(_startdist);
        // _enddist = flip_match ? match_end_base_end_dist : match_start_base_end_dist;
        // _enddist = std::sqrt(_enddist);

        /// Ask how far the matched start point is from the line
        ::geoalgo::Point_t match_vtx = flip_match ?
                                       ::geoalgo::Point_t(match_candidate.End()) :
                                       ::geoalgo::Point_t(match_candidate.Vertex());

        _infdist = _geoalg.SqDist(base_line, match_vtx);
        _infdist = std::sqrt(_infdist);

        ::geoalgo::Point_t projpoint = _geoalg.ClosestPt(base_line, match_vtx);
        _endprojdist = ::geoalgo::Point_t(base_trk.End()).Dist(projpoint);
        // std::cout<<"Dist between "<<::geoalgo::Point_t(base_trk.End())
        // <<" and "<<projpoint<<" is "<<_endprojdist<<std::endl;
        // std::cout<<"infdist computed to be "<<_infdist<<", computing by hand gives "
        // << projpoint.Dist(match_vtx) << std::endl;
        // std::cout<<"reversing order of dist gives "<<match_vtx.Dist(projpoint)<<std::endl;
        // std::cout<<"computing dist myself vies "<<(match_vtx-projpoint).Length()<<std::endl;
        // std::cout<<"projpoint is "<<projpoint<<std::endl;
        // std::cout<<"line is "<<base_line<<std::endl;
        // std::cout<<"vertex is "<<match_vtx<<std::endl;
        // std::cout<<"Dist from projpoint ot the line (should be zero): "
        // << std::sqrt(_geoalg.SqDist(base_line,projpoint)) <<std::endl;
        _startprojdist = ::geoalgo::Point_t(base_trk.Vertex()).Dist(projpoint);

        _debug_tree->Fill();

        if (_infdist < 10. && fabs(_dotprod) > 0.999 &&
                _startprojdist > _baselen && _endprojdist < _startprojdist)
            return true;

        return false;
    }

    larlite::track KalekoTrackStitcher::stitchTracks(const larlite::track &base_trk,
            const larlite::track &match_candidate ) {

        larlite::track result;
        result.clear_data();

        result.set_track_id(base_trk.ID());

        size_t n_pts_base = base_trk.NumberTrajectoryPoints();
        size_t n_pts_match = match_candidate.NumberTrajectoryPoints();
        result.reserve(n_pts_base + n_pts_match);
        // Loop over base track and copy it exactly
        for (size_t i = 0; i < n_pts_base; ++i) {
            result.add_vertex(base_trk.LocationAtPoint(i));
            result.add_direction(base_trk.DirectionAtPoint(i));
        }
        // Loop over matched track and add those trajectory points on to the end
        for (size_t i = 0; i < n_pts_match; ++i) {
            result.add_vertex(match_candidate.LocationAtPoint(i));
            result.add_direction(match_candidate.DirectionAtPoint(i));
        }

        return result;
    }

}
#endif
