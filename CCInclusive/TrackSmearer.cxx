#ifndef LARLITE_TRACKSMEARER_CXX
#define LARLITE_TRACKSMEARER_CXX

#include "TrackSmearer.h"

namespace larlite {


    larlite::track TrackSmearer::SmearTrack(const larlite::track &trk) {

        double gaus_center = 0.;
        double gaus_width = 0.0150; //cm
        double small_gaus_width = 0.0075;
        double middle_gaus_width = 0.0125;

        double first_bad_z = 200.;
        double mid_z_min = 200.;
        double mid_z_max = 300.;
        double bad_z_min = 675.;
        double bad_z_max = 775.;

        larlite::track result;
        result.set_track_id(trk.ID());

        size_t n_points = trk.NumberTrajectoryPoints();

        // Loop over track start-to-end and smear points
        // that are in the bad Z region
        for (int i = 0; i < n_points - 1; ++i) {
            auto traj_pt = trk.LocationAtPoint(i);
            if ( traj_pt.Z() > bad_z_min && traj_pt.Z() < bad_z_max ) {
                double smeared_x = traj_pt.X() + _trandom.Gaus(gaus_center, middle_gaus_width);
                double smeared_y = traj_pt.Y() + _trandom.Gaus(gaus_center, middle_gaus_width);
                double smeared_z = traj_pt.Z() + _trandom.Gaus(gaus_center, middle_gaus_width);
                TVector3 smeared_traj_pt = TVector3(smeared_x, smeared_y, smeared_z);
                result.add_vertex(smeared_traj_pt);
                result.add_direction(smeared_traj_pt);
            }
            else if ( traj_pt.Z() > mid_z_min && traj_pt.Z() < mid_z_max ) {
                double smeared_x = traj_pt.X() + _trandom.Gaus(gaus_center, middle_gaus_width);
                double smeared_y = traj_pt.Y() + _trandom.Gaus(gaus_center, middle_gaus_width);
                double smeared_z = traj_pt.Z() + _trandom.Gaus(gaus_center, middle_gaus_width);
                TVector3 smeared_traj_pt = TVector3(smeared_x, smeared_y, smeared_z);
                result.add_vertex(smeared_traj_pt);
                result.add_direction(smeared_traj_pt);
            }
            else if ( traj_pt.Z() < first_bad_z ){
                double smeared_x = traj_pt.X() + _trandom.Gaus(gaus_center, gaus_width);
                double smeared_y = traj_pt.Y() + _trandom.Gaus(gaus_center, gaus_width);
                double smeared_z = traj_pt.Z() + _trandom.Gaus(gaus_center, gaus_width);
                TVector3 smeared_traj_pt = TVector3(smeared_x, smeared_y, smeared_z);
                result.add_vertex(smeared_traj_pt);
                result.add_direction(smeared_traj_pt);
            }
            else {
                double smeared_x = traj_pt.X() + _trandom.Gaus(gaus_center, small_gaus_width);
                double smeared_y = traj_pt.Y() + _trandom.Gaus(gaus_center, small_gaus_width);
                double smeared_z = traj_pt.Z() + _trandom.Gaus(gaus_center, small_gaus_width);
                TVector3 smeared_traj_pt = TVector3(smeared_x, smeared_y, smeared_z);
                result.add_vertex(smeared_traj_pt);
                result.add_direction(smeared_traj_pt);
            }
        }

        return result;
    }

}
#endif
