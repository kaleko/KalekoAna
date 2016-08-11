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

}
#endif
