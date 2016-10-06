#ifndef LARLITE_TRACKCHOPPER_CXX
#define LARLITE_TRACKCHOPPER_CXX

#include "TrackChopper.h"

namespace larlite {

    larlite::track TrackChopper::chopTrack(const larlite::track &trk) {

        // Loop over trajectory points and build a new track with only the points inside of the fid vol box

        // If the track is fully contained (which should never happen)
        // Just return the original track
        if (_box.Contain(trk.Vertex()) && _box.Contain(trk.End())) {
            // std::cout << " Hrmmm... why did TrackChopper get handed a fully-contained track? Returning original ... " << std::endl;
            return larlite::track(trk);
        }

        larlite::track result;
        result.set_track_id(trk.ID());

        size_t n_points = trk.NumberTrajectoryPoints();

        // If the front of the track is contained and the back isn't,
        // loop over trajectory points from front-to-back until you reach a point
        // that is not contained in the fiducial volume. Then stop and return
        // the portion of the track inside of the fiducial volume.
        if (_box.Contain(trk.Vertex()) && !_box.Contain(trk.End())) {
            for (int i = 0; i < n_points; ++i) {
                if (_box.Contain(trk.LocationAtPoint(i))) {
                    result.add_vertex(trk.LocationAtPoint(i));
                    // You also need to add directions just to get the NumberTrajectoryPoints to work
                    // because the reco track data product is ass-backwards
                    // so I'm filling it with a dummy now
                    result.add_direction(trk.LocationAtPoint(i));
                }
                else
                    return result;
            }
        }
        // If the back of the track is contained and the front isn't,
        // then the track is backwards, so we loop over in reverse order.
        // NOTE this will flip the direction of the track, so you don't need
        // to worry about flipping it later on down the line.
        else if (_box.Contain(trk.End()) && !_box.Contain(trk.Vertex())) {
            for (int i = n_points - 1; i >= 0; i--) {
                if (_box.Contain(trk.LocationAtPoint(i))) {
                    result.add_vertex(trk.LocationAtPoint(i));
                    // You also need to add directions just to get the NumberTrajectoryPoints to work
                    // because the reco track data product is ass-backwards
                    // so I'm filling it with a dummy now
                    result.add_direction(trk.LocationAtPoint(i));
                }
                else
                    return result;
            }
        }
        else {
            // std::cout << "Neither end of this track is contained. Returning an exact copy of the track."<<std::endl;
            return trk;
        }

        return result;
    }

    larlite::track TrackChopper::chopAndStraightenTrack(const larlite::track &trk) {

        bool debug = false;
        if (debug)
            std::cout << " Begin of chopAndStraightenTrack. start track has number of points = "
                      << trk.NumberTrajectoryPoints() << " and length " << (trk.Vertex() - trk.End()).Mag()
                      << ", with start Z " << trk.Vertex().Z() << " and end Z " << trk.End().Z()
                      << std::endl;


        double bad_z_min = 650.;
        double bad_z_max = 800.;

        larlite::track result;
        result.set_track_id(trk.ID());

        // Start by chopping the track. Then
        // erase any points inside of the bad z region so you get a straight line
        // There are four separate cases:
        // 1) Both start and end points of track are inside of the bad region
        // 2) The start point is outside of the bad z- region but end point is inside
        // 3) The end point is outside of the bad z- region but start point is inside
        // 4) Both the start and end points are outside of the bad region

        larlite::track chopped_trk = chopTrack(trk);

        if (debug)
            std::cout << " After initial chopping. track has number of points = "
                      << chopped_trk.NumberTrajectoryPoints() << " and length "
                      << (chopped_trk.Vertex() - chopped_trk.End()).Mag()
                      << ", with start Z " << chopped_trk.Vertex().Z() << " and end Z " << chopped_trk.End().Z() << std::endl;

        size_t n_points = chopped_trk.NumberTrajectoryPoints();

        // 1) If both start and end points are in the bad z- region, keep only the start and end point.
        if (chopped_trk.Vertex().Z() < bad_z_max && chopped_trk.Vertex().Z() > bad_z_min &&
                chopped_trk.End().Z() < bad_z_max && chopped_trk.End().Z() > bad_z_min) {

            result.add_vertex(chopped_trk.Vertex());
            result.add_direction(chopped_trk.Vertex());
            result.add_vertex(chopped_trk.End());
            result.add_direction(chopped_trk.End());
            if (debug)
                std::cout << " case 1. final track has number of points = "
                          << result.NumberTrajectoryPoints() << " and length " << (result.Vertex() - result.End()).Mag()
                          << ", with start Z " << result.Vertex().Z() << " and end Z " << result.End().Z() << std::endl;

            return result;
        }

        // 2) If the start point is outside of the bad z- region but end point is inside,
        if ( (chopped_trk.Vertex().Z() < bad_z_min || chopped_trk.Vertex().Z() > bad_z_max) &&
                (chopped_trk.End().Z() > bad_z_min && chopped_trk.End().Z() < bad_z_max) ) {

            // Loop over track start-to-end and ignore points
            // that are in the bad Z region
            // (n_points-1 is to skip the end point, which we know is in the bad region)
            for (int i = 0; i < n_points - 1; ++i) {
                auto traj_pt = chopped_trk.LocationAtPoint(i);
                if (traj_pt.Z() < bad_z_min || traj_pt.Z() > bad_z_max) {
                    result.add_vertex(traj_pt);
                    result.add_direction(traj_pt);
                }
            }
            // Don't forget to add the z- end point that is in the bad region!
            result.add_vertex(chopped_trk.End());
            result.add_direction(chopped_trk.End());

            if (debug)
                std::cout << " case 2. final track has number of points = "
                          << result.NumberTrajectoryPoints() << " and length " << (result.Vertex() - result.End()).Mag()
                          << ", with start Z " << result.Vertex().Z() << " and end Z " << result.End().Z() << std::endl;

            return result;
        }

        // 3) If the end point is outside of the bad z- region but start point is inside,
        if ( (chopped_trk.End().Z() < bad_z_min || chopped_trk.End().Z() > bad_z_max) &&
                (chopped_trk.Vertex().Z() > bad_z_min && chopped_trk.Vertex().Z() < bad_z_max) ) {

            // Add the start point that is in the bad region!
            result.add_vertex(chopped_trk.Vertex());
            result.add_direction(chopped_trk.Vertex());

            // Loop over track start-to-end and ignore points
            // that are in the bad Z region
            for (int i = 0; i < n_points - 1; ++i) {
                auto traj_pt = chopped_trk.LocationAtPoint(i);
                if (traj_pt.Z() < bad_z_min || traj_pt.Z() > bad_z_max) {
                    result.add_vertex(traj_pt);
                    result.add_direction(traj_pt);
                }
            }

            if (debug)
                std::cout << " case 3. final track has number of points = "
                          << result.NumberTrajectoryPoints() << " and length " << (result.Vertex() - result.End()).Mag()
                          << ", with start Z " << result.Vertex().Z() << " and end Z " << result.End().Z() << std::endl;

            return result;
        }

        // 4) Both the start and end points are outside of the bad region
        // Just return the chopped track since it's already good
        return chopped_trk;
    }

}
#endif
