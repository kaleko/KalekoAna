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

        print(larlite::msg::kNORMAL, __FUNCTION__,
              Form("Initializing with base producer: %s, match producer: %s.",
                   _base_producer.c_str(),
                   _match_producer.c_str()));

        _geoalg = geoalgo::GeoAlgo();

        if (!_h_lendiff)
            _h_lendiff = new TH1F("h_lendiff", "After Stitch Len - Before Stitch Len", 1000, -10, 500);

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

        // Read in match track producer
        auto ev_match_track = storage->get_data<event_track>(_match_producer);
        if (!ev_match_track) {
            print(larlite::msg::kERROR, __FUNCTION__,
                  Form("Did not find specified data product, track by %s!", _match_producer.c_str()));
            return false;
        }

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
            for (size_t match_idx = 0; match_idx < ev_match_track->size(); ++match_idx) {
                auto const& match_trk = ev_match_track->at(match_idx);

                // Check if this match candidate track has already been matched to a base
                bool already_matched = false;
                std::map<size_t, std::pair<size_t, bool> >::const_iterator it;
                for (it = idx_map.begin(); it != idx_map.end(); ++it) {
                    if (it->second.first == match_idx) already_matched = true;
                }
                if (already_matched) {
                    // This candidate was already matched to a base... skip it
                    continue;
                }

                if (!trackViable(match_trk)) continue;

                std::pair<bool, bool> matched = tracksMatched(base_trk, match_trk);

                if (matched.first) {
                    idx_map.insert(
                        std::pair<size_t,
                        std::pair<size_t, bool> >(base_idx, std::make_pair(match_idx, matched.second))
                    );
                }
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
            std::map<size_t, std::pair<size_t, bool> >::iterator it;
            it = idx_map.find(base_idx);
            // If so:
            if (it != idx_map.end()) {
                ev_stitched_track->push_back(
                    stitchTracks(ev_base_track->at(base_idx),
                                 ev_match_track->at(it->second.first),
                                 it->second.second)
                );
                double before_len = (ev_base_track->at(base_idx).End() - ev_base_track->at(base_idx).Vertex()).Mag();
                double after_len = (ev_stitched_track->back().End() - ev_stitched_track->back().Vertex()).Mag();
                _h_lendiff->Fill(after_len - before_len);

            }

            // If not
            else
                ev_stitched_track->push_back(ev_base_track->at(base_idx));



        }// End loop over base tracks again
        return true;
    }

    bool KalekoTrackStitcher::finalize() {

        if (_fout) {
            if (_debug_tree) {
                _fout->cd();
                _debug_tree->Write();
                _h_lendiff->Write();
            }
        }

        return true;
    }

    bool KalekoTrackStitcher::trackViable(const larlite::track &trk) {

        if ( trk.NumberTrajectoryPoints() < 2 ) return false;

        double _trk_len = (trk.End() - trk.Vertex()).Mag();
        if ( _trk_len < 3.) return false;

        return true;
    }

    std::pair<bool, bool> KalekoTrackStitcher::tracksMatched(const larlite::track &base_trk,
            const larlite::track &match_candidate ) {

        /// Create a GeoLine for the base track...
        ///   >> one point on the line is track start
        ///   >> one point on the line is track end
        ///   >> line extends infinitely in both directions
        ::geoalgo::Line base_line(
            ::geoalgo::Point_t(base_trk.Vertex()),
            ::geoalgo::Point_t(base_trk.End()));

        /// Also create a GeoLineSegment for the base track...
        ///   >> one point on the line is track start
        ///   >> one point on the line is track end
        ///   >> line extends from start point to end point
        ::geoalgo::LineSegment base_line_seg(
            ::geoalgo::Point_t(base_trk.Vertex()),
            ::geoalgo::Point_t(base_trk.End()));

        // Decide if the match track should be flipped
        TVector3 base_dir = base_trk.End() - base_trk.Vertex();
        TVector3 match_dir = match_candidate.End() - match_candidate.Vertex();
        _baselen = base_dir.Mag();
        _dotprod = base_dir.Unit().Dot(match_dir.Unit());
        bool flip_match = _dotprod > 0 ? false : true;

        /// Ask how far the matched start point is from the line
        ::geoalgo::Point_t match_vtx = flip_match ?
                                       ::geoalgo::Point_t(match_candidate.End()) :
                                       ::geoalgo::Point_t(match_candidate.Vertex());

        _infdist = _geoalg.SqDist(base_line, match_vtx);
        _infdist = std::sqrt(_infdist);

        ::geoalgo::Point_t projpoint = _geoalg.ClosestPt(base_line, match_vtx);
        _endprojdist = ::geoalgo::Point_t(base_trk.End()).Dist(projpoint);
        _startprojdist = ::geoalgo::Point_t(base_trk.Vertex()).Dist(projpoint);

        _debug_tree->Fill();

        if (_infdist < 10. && fabs(_dotprod) > 0.999)
            return std::make_pair(true, flip_match);

        return std::make_pair(false, flip_match);
    }

    larlite::track KalekoTrackStitcher::stitchTracks(const larlite::track &base_trk,
            const larlite::track &match_candidate, bool flip_match ) {

        // First build a geotrajectory from the tracks, then turn it into a reco track
        // This is done because sometimes I want to add points to the *front* of the output
        // track, which isn't possible with an actual reco track object, but is easy with
        // a geotrajectory

        // First, copy the the base track into my geotrajectory
        ::geoalgo::Trajectory geotraj;
        size_t n_pts_base = base_trk.NumberTrajectoryPoints();
        for (size_t i = 0; i < n_pts_base; ++i) {
            geotraj.push_back(::geoalgo::Vector(base_trk.LocationAtPoint(i)));
        }

        /// Create a GeoLineSegment for the base track...
        ///   >> one point on the line is track start
        ///   >> one point on the line is track end
        ///   >> line extends between these two points.
        ::geoalgo::LineSegment base_line_seg(
            ::geoalgo::Point_t(base_trk.Vertex()),
            ::geoalgo::Point_t(base_trk.End()));

        // Loop over matched track, for each trajectory point make sure
        // its projection onto the base track line segment is not in the
        // already-existing range of base-track points, only add them
        // to the output track if not
        // If the point is farther than the end of the base track,
        // add it to the end. If the point is before the start of the
        // base track, add it to the start.

        size_t n_pts_match = match_candidate.NumberTrajectoryPoints();
        if (!flip_match) {
            for (size_t i = 0; i < n_pts_match; ++i) {
                ::geoalgo::Vector ipt = ::geoalgo::Vector(match_candidate.LocationAtPoint(i));
                auto pt = _geoalg.ClosestPt(base_line_seg, ipt);
                // If the point is before the front of the trajectory
                if (pt.SqDist(::geoalgo::Vector(base_trk.Vertex())) < 0.001 )
                    geotraj.insert(geotraj.begin(), ipt);

                // If the point is after the end of the trajectory
                else if (pt.SqDist(::geoalgo::Vector(base_trk.End())) < 0.001 )
                    geotraj.push_back(ipt);

                // If the point is overlapping the already existing base track
                else
                    continue;

            }
        }
        else {
            for (int i = n_pts_match - 1; i >= 0; i--) {

                ::geoalgo::Vector ipt = ::geoalgo::Vector(match_candidate.LocationAtPoint(i));
                auto pt = _geoalg.ClosestPt(base_line_seg, ipt);
                // If the point is before the front of the trajectory
                if (pt.SqDist(::geoalgo::Vector(base_trk.Vertex())) < 0.001 )
                    geotraj.insert(geotraj.begin(), ipt);

                // If the point is after the end of the trajectory
                else if (pt.SqDist(::geoalgo::Vector(base_trk.End())) < 0.001 )
                    geotraj.push_back(ipt);

                // If the point is overlapping the already existing base track
                else
                    continue;

            }
        }


        // Now convert the geotrajectory into a reco track
        larlite::track result;
        result.clear_data();

        result.set_track_id(base_trk.ID());
        result.reserve(geotraj.size());
        // Loop over base track and copy it exactly
        for (size_t i = 0; i < geotraj.size() - 1; ++i) {
            result.add_vertex(TVector3(geotraj.at(i).at(0), geotraj.at(i).at(1), geotraj.at(i).at(2)));
            result.add_direction(TVector3(geotraj.at(i + 1).at(0) - geotraj.at(i).at(0),
                                          geotraj.at(i + 1).at(1) - geotraj.at(i).at(1),
                                          geotraj.at(i + 1).at(2) - geotraj.at(i).at(2)));
        }

        return result;
    }

}
#endif
