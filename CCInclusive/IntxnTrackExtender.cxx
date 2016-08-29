#ifndef LARLITE_INTXNTRACKEXTENDER_CXX
#define LARLITE_INTXNTRACKEXTENDER_CXX

#include "IntxnTrackExtender.h"

namespace larlite {

    void IntxnTrackExtender::ExtendVertexTracks( larlite::KalekoNuItxn & itxn, const larlite::event_track * ev_track) {

        auto const &geovtx = ::geoalgo::Vector(itxn.Vertex().X(),
                                               itxn.Vertex().Y(),
                                               itxn.Vertex().Z());

        std::set<size_t> extended_trk_ids;

        // Loop over tracks associated with the vertex of the interaction
        for (size_t ivtxtrk = 0; ivtxtrk < itxn.Tracks().size(); ++ivtxtrk) {
            auto vtxtrk = itxn.Tracks().at(ivtxtrk);

            bool flip_vtx_trk = ::geoalgo::Vector(vtxtrk.Vertex()).SqDist(geovtx) <
                                ::geoalgo::Vector(vtxtrk.End()).SqDist(geovtx) ?
                                false : true;

            // Loop over potential matched tracks
            for (auto const& match_trk : *ev_track) {
                // If this track has already been stitched to a vertex track, skip it
                if (extended_trk_ids.find(match_trk.ID()) != extended_trk_ids.end()) continue;

                auto tracks_matched = tracksMatched(vtxtrk, flip_vtx_trk, match_trk);
                if (tracks_matched.first) {
                    auto new_vtxtrk = stitchTracks(vtxtrk, match_trk, flip_vtx_trk, tracks_matched.second);
                    itxn.ReplaceTrack(ivtxtrk, new_vtxtrk);
                    extended_trk_ids.insert(match_trk.ID());
                }
            }
        }
    }

    std::pair<bool, bool> IntxnTrackExtender::tracksMatched(const larlite::track & base_trk, bool flip_base_trk,
            const larlite::track & match_candidate ) {

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
        if (!flip_base_trk)
            ::geoalgo::LineSegment base_line_seg(
                ::geoalgo::Point_t(base_trk.Vertex()),
                ::geoalgo::Point_t(base_trk.End()));
        else
            ::geoalgo::LineSegment base_line_seg(
                ::geoalgo::Point_t(base_trk.End()),
                ::geoalgo::Point_t(base_trk.Vertex()));

        // Decide if the match track should be flipped
        TVector3 base_dir = base_trk.End() - base_trk.Vertex();
        TVector3 match_dir = match_candidate.End() - match_candidate.Vertex();
        _dotprod = base_dir.Unit().Dot(match_dir.Unit());
        bool flip_match = _dotprod > 0 ? false : true;

        /// Ask how far the matched start point is from the line
        ::geoalgo::Point_t match_vtx          = ::geoalgo::Point_t(match_candidate.Vertex());
        ::geoalgo::Point_t match_vtx_inverted = ::geoalgo::Point_t(match_candidate.End());

        _infdist = std::min(_geoalg.SqDist(base_line, match_vtx), _geoalg.SqDist(base_line, match_vtx_inverted));
        _infdist = std::sqrt(_infdist);

        ::geoalgo::Point_t projpoint = _geoalg.ClosestPt(base_line, match_vtx);

        if (_infdist < 30. && fabs(_dotprod) > 0.9 )
            return std::make_pair(true, flip_match);

        return std::make_pair(false, flip_match);
    }


    larlite::track IntxnTrackExtender::stitchTracks(const larlite::track &base_trk,
            const larlite::track &match_candidate, bool flip_base, bool flip_match ) {

        // First build a geotrajectory from the tracks, then turn it into a reco track
        // This is done because sometimes I want to add points to the *front* of the output
        // track, which isn't possible with an actual reco track object, but is easy with
        // a geotrajectory

        // First, copy the the base track into my geotrajectory
        ::geoalgo::Trajectory geotraj;
        size_t n_pts_base = base_trk.NumberTrajectoryPoints();
        if (!flip_base)
            for (size_t i = 0; i < n_pts_base; ++i)
                geotraj.push_back(::geoalgo::Vector(base_trk.LocationAtPoint(i)));

        else
            for (int i = n_pts_base - 1; i >= 0; i--)
                geotraj.push_back(::geoalgo::Vector(base_trk.LocationAtPoint(i)));

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
