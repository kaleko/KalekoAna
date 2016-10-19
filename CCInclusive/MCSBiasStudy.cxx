#ifndef MCSBIASSTUDY_CXX
#define MCSBIASSTUDY_CXX

#include "MCSBiasStudy.h"
#include "DataFormat/track.h"
#include "DataFormat/mctrack.h"

namespace larlite {

    MCSBiasStudy::MCSBiasStudy() {

        _myspline = TrackMomentumSplines();
        _tmc = 0;
        _tmc = new kaleko::TrackMomentumCalculator();
        _tmc->SetStepSize(20.0);

        _chopper = TrackChopper();
        _smearer = TrackSmearer();

        _tree = new TTree("MCS_bias_tree", "MCS_bias_tree");
        _tree->Branch("length_analyzed", &_length_analyzed, "length_analyzed/D");
        _tree->Branch("full_length", &_full_length, "full_length/D");
        _tree->Branch("chopped_full_length", &_chopped_full_length, "chopped_full_length/D");
        _tree->Branch("full_range_energy", &_full_range_energy, "full_range_energy/D");
        _tree->Branch("MCS_energy", &_MCS_energy, "MCS_energy/D");
        _tree->Branch("full_MCS_energy", &_full_MCS_energy, "full_MCS_energy/D");
        _tree->Branch("full_MCS_energy_someflipped", &_full_MCS_energy_someflipped, "full_MCS_energy_someflipped/D");
        _tree->Branch("full_MCS_energy_chopped", &_full_MCS_energy_chopped, "full_MCS_energy_chopped/D");
        _tree->Branch("full_MCS_energy_smeared", &_full_MCS_energy_smeared, "full_MCS_energy_smeared/D");
        _tree->Branch("full_MCS_energy_reweighted", &_full_MCS_energy_reweighted, "full_MCS_energy_reweighted/D");
        _tree->Branch("full_MCS_energy_chopped_smeared", &_full_MCS_energy_chopped_smeared, "full_MCS_energy_chopped_smeared/D");
        _tree->Branch("track_start_x", &_track_start_x, "track_start_x/D");
        _tree->Branch("track_start_y", &_track_start_y, "track_start_y/D");
        _tree->Branch("track_start_z", &_track_start_z, "track_start_z/D");
        _tree->Branch("track_end_x", &_track_end_x, "track_end_x/D");
        _tree->Branch("track_end_y", &_track_end_y, "track_end_y/D");
        _tree->Branch("track_end_z", &_track_end_z, "track_end_z/D");
        _tree->Branch("track_dot_z", &_track_dot_z, "track_dot_z/D");
        _tree->Branch("n_traj_points", &_n_traj_points, "n_traj_points/I");
        _tree->Branch("chopped_wiggle", &_chopped_wiggle, "chopped_wiggle/D");
        _tree->Branch("chopped_std", &_chopped_std, "chopped_std/D");
        _tree->Branch("smeared_wiggle", &_smeared_wiggle, "smeared_wiggle/D");
        _tree->Branch("smeared_std", &_smeared_std, "smeared_std/D");
        _tree->Branch("chopped_smeared_wiggle", &_chopped_smeared_wiggle, "chopped_smeared_wiggle/D");
        _tree->Branch("chopped_smeared_std", &_chopped_smeared_std, "chopped_smeared_std/D");
        _tree->Branch("long_curve_dotprod", &_long_curve_dotprod, "long_curve_dotprod/D");
        _tree->Branch("full_track_tree_entry", &_full_track_tree_entry, "full_track_tree_entry/O");

        _seg_tree = new TTree("MCS_segment_tree", "MCS_segment_tree");
        _seg_tree->Branch("seg_start_x", &_seg_start_x, "seg_start_x/D");
        _seg_tree->Branch("seg_start_y", &_seg_start_y, "seg_start_y/D");
        _seg_tree->Branch("seg_start_z", &_seg_start_z, "seg_start_z/D");
        _seg_tree->Branch("seg_end_x", &_seg_end_x, "seg_end_x/D");
        _seg_tree->Branch("seg_end_y", &_seg_end_y, "seg_end_y/D");
        _seg_tree->Branch("seg_end_z", &_seg_end_z, "seg_end_z/D");
        _seg_tree->Branch("n_seg_traj_points", &_n_seg_traj_points, "n_seg_traj_points/I");
        _seg_tree->Branch("seg_avg_perp_dist", &_seg_avg_perp_dist, "seg_avg_perp_dist/D");
        _seg_tree->Branch("seg_std_perp_dist", &_seg_std_perp_dist, "seg_std_perp_dist/D");

    }


    void MCSBiasStudy::AnalyzeTrack(const larlite::mctrack &mct) {

        // _full_length = (mct.back().Position().Vect() - mct.front().Position().Vect()).Mag();
        // _full_range_energy = _myspline.GetMuMomentum(_full_length) / 1000. + 0.106;
        _full_MCS_energy_chopped = _tmc->GetMomentumMultiScatterLLHD(mct);

    }

    void MCSBiasStudy::AnalyzeTrack(const larlite::track &track) {

        bool apply_reweight = false;
        _full_MCS_energy = _tmc->GetMomentumMultiScatterLLHD(track, false, true, apply_reweight);

        // _length_analyzed = -999.;
        // _MCS_energy = -999.;
        //_full_MCS_energy = _tmc->GetMomentumMultiScatterLLHD(track, false, false, false);
        //_full_MCS_energy_reweighted = _tmc->GetMomentumMultiScatterLLHD(track, false, false, true);
        // _full_MCS_energy_smeared = -999.;
        // _full_MCS_energy_chopped_smeared = -999.;
        // _full_MCS_energy_chopped = -999.;
        // _chopped_wiggle = -999.;
        // _chopped_std = -999.;
        // _chopped_smeared_wiggle = -999.;
        // _chopped_smeared_std = -999.;
        // _smeared_wiggle = -999.;
        // _smeared_std = -999.;

        // _full_length = (track.End() - track.Vertex()).Mag();
        // _full_range_energy = _myspline.GetMuMomentum(_full_length) / 1000. + 0.106;
        // _full_MCS_energy = _tmc->GetMomentumMultiScatterLLHD(track);
        // _track_start_x = track.Vertex().X();
        // _track_start_y = track.Vertex().Y();
        // _track_start_z = track.Vertex().Z();
        // _track_end_x = track.End().X();
        // _track_end_y = track.End().Y();
        // _track_end_z = track.End().Z();
        // _full_track_tree_entry = true;
        // _n_traj_points = track.NumberTrajectoryPoints();

        // bool flip_trk = _track_start_y < _track_end_y;
        // // _full_MCS_energy_someflipped = _tmc->GetMomentumMultiScatterLLHD(track, flip_trk);

        // TVector3 zdir(0., 0., -1.);
        // _track_dot_z = !flip_trk ? (track.End() - track.Vertex()).Unit().Dot(zdir) :
        //                (track.Vertex() - track.End()).Unit().Dot(zdir);

        // // auto const chopped_trk = _chopper.chopTrack(track);
        // // _chopped_full_length = (chopped_trk.End() - chopped_trk.Vertex()).Mag();

        // // if (_chopped_full_length > 100.)
        // //     _full_MCS_energy_chopped = _tmc->GetMomentumMultiScatterLLHD(chopped_trk, false, true);

        // // auto const smeared_trk = _smearer.SmearTrack(track);
        // // _full_MCS_energy_smeared = _tmc->GetMomentumMultiScatterLLHD(smeared_trk);
        // // auto const chopped_smeared_trk = _smearer.SmearTrack(chopped_trk);
        // // _full_MCS_energy_chopped_smeared = _tmc->GetMomentumMultiScatterLLHD(chopped_smeared_trk);

        // // auto wiggle = ComputeWiggle(chopped_trk, 10.);
        // // _chopped_wiggle = wiggle.first;
        // // _chopped_std = wiggle.second;

        // // auto smeared_wiggle = ComputeWiggle(smeared_trk, 10.);
        // // _smeared_wiggle = smeared_wiggle.first;
        // // _smeared_std = smeared_wiggle.second;

        // // auto chopped_smeared_wiggle = ComputeWiggle(chopped_smeared_trk, 10.);
        // // _chopped_smeared_wiggle = chopped_smeared_wiggle.first;
        // // _chopped_smeared_std = chopped_smeared_wiggle.second;

        // _long_curve_dotprod = -999.;
        // // if (_chopped_full_length > 200.)
        // //     _long_curve_dotprod = ComputeLongCurve(chopped_trk);

        _tree->Fill();


        // _full_track_tree_entry = false;


        // // // slowly build up a copy of the initial track with increasing length
        // // // and compute MCS energy as length increases
        // // larlite::track dummy_trk;

        // // double length_increment = 10.; //cm

        // // double current_min_len = 100.;
        // // size_t n_points = track.NumberTrajectoryPoints();

        // // for (int i = 1; i < n_points; ++i) {

        // //   dummy_trk.add_vertex(track.LocationAtPoint(i - 1));
        // //   dummy_trk.add_direction(track.LocationAtPoint(i - 1));

        // //   _length_analyzed = (track.LocationAtPoint(i) - track.LocationAtPoint(0)).Mag();
        // //   if ( _length_analyzed < current_min_len ) continue;

        // //   _MCS_energy = _tmc->GetMomentumMultiScatterLLHD(dummy_trk);

        // //   current_min_len += length_increment;

        // //   _tree->Fill();

        // // }

    }

    std::pair<double, double> MCSBiasStudy::ComputeWiggle(const larlite::track &track, double seg_size) {

        _seg_start_x = -999;
        _seg_start_y = -999;
        _seg_start_z = -999;
        _seg_end_x = -999;
        _seg_end_y = -999;
        _seg_end_z = -999;
        _n_seg_traj_points = -999;
        _seg_avg_perp_dist = -999;
        _seg_std_perp_dist = -999;

        if (track.NumberTrajectoryPoints() < 2) return std::make_pair(-1., -1.);
        double seg_size_squared = seg_size * seg_size;//100.; //cm^2
        TVector3 start = track.Vertex();
        TVector3 seg_end = start;

        TVector3 current_seg_dir;
        bool found_first_segment = 0;
        size_t seg_counter = 0;
        std::vector<double> angle_deviations;

        // Store segment start indices for second loop later
        std::vector<size_t> seg_start_indices;
        seg_start_indices.clear();
        seg_start_indices.push_back(0);

        for (size_t i = 0; i < track.NumberTrajectoryPoints(); ++i) {

            TVector3 point = track.LocationAtPoint(i);

            // New segment
            if ((point - seg_end).Mag2() > seg_size_squared) {
                seg_counter++;
                seg_start_indices.push_back(i);

                if (!found_first_segment) {
                    found_first_segment = true;
                    current_seg_dir = point - seg_end;
                    continue;
                }

                double this_deviation_dot = current_seg_dir.Unit().Dot((point - seg_end).Unit());
                //acos returns NAN for numbers higher than about this
                angle_deviations.push_back(this_deviation_dot > 0.999999 ? 0 : std::acos(this_deviation_dot));
                current_seg_dir = point - seg_end;
                seg_end = point;
            } // End new segment

        } // End loop over trajectory points

        double avg_angle_deviation = 0.;
        double std_angle_deviation = 0.;
        for (size_t i = 0; i < angle_deviations.size(); ++i)
            avg_angle_deviation += angle_deviations.at(i);
        avg_angle_deviation /= angle_deviations.size();
        for (size_t i = 0; i < angle_deviations.size(); ++i)
            std_angle_deviation += std::pow((angle_deviations.at(i) - avg_angle_deviation), 2);
        std_angle_deviation /= angle_deviations.size();
        std_angle_deviation = std::sqrt(std_angle_deviation);


        /// second loop over traj points to compute segment-specific variables for the segment ttree
        for (size_t i = 0; i < seg_start_indices.size() - 1; ++i) {
            // Build a geo line segment from start point to end point
            TVector3 seg_start_point = track.LocationAtPoint(seg_start_indices[i]);
            TVector3 seg_end_point   = track.LocationAtPoint(seg_start_indices[i + 1]);
            ::geoalgo::LineSegment geo_seg =
                ::geoalgo::LineSegment(seg_start_point.X(), seg_start_point.Y(), seg_start_point.Z(),
                                       seg_end_point.X(), seg_end_point.Y(), seg_end_point.Z());

            _seg_start_x = seg_start_point.X();
            _seg_start_y = seg_start_point.Y();
            _seg_start_z = seg_start_point.Z();
            _seg_end_x = seg_end_point.X();
            _seg_end_y = seg_end_point.Y();
            _seg_end_z = seg_end_point.Z();

            // Loop over trajectory points in the segment
            std::vector<double> perp_dists;
            for (size_t j = seg_start_indices[i]; j < seg_start_indices[i + 1]; ++j) {
                // ::geoalgo::Point geo_pt = ::geoalgo::Point(track.LocationAtPoint(j));
                double perp_dist = _geoalg.SqDist(geo_seg, track.LocationAtPoint(j));
                perp_dists.push_back(perp_dist);
            }
            _seg_avg_perp_dist = 0.;
            _seg_std_perp_dist = 0.;
            for (size_t k = 0; k < perp_dists.size(); ++k)
                _seg_avg_perp_dist += perp_dists.at(k);
            _seg_avg_perp_dist /= perp_dists.size();
            for (size_t k = 0; k < perp_dists.size(); ++k)
                _seg_std_perp_dist += std::pow((perp_dists.at(k) - _seg_avg_perp_dist), 2);
            _seg_std_perp_dist /= perp_dists.size();
            _seg_std_perp_dist = std::sqrt(_seg_std_perp_dist);

            _n_seg_traj_points = seg_start_indices[i + 1] - seg_start_indices[i];
            _seg_tree->Fill();
        }


        return seg_counter > 0 ? std::make_pair( avg_angle_deviation, std_angle_deviation ) : std::make_pair(-1., -1.);
    }// End ComputeWiggle



    double MCSBiasStudy::ComputeLongCurve(const larlite::track &track) {

        if (track.NumberTrajectoryPoints() < 2) return -999.;

        TVector3 start = track.Vertex();
        TVector3 end = track.End();
        TVector3 first_50cm_dir;

        if ( (start - end).Mag() < 200. ) return -999.;

        for (size_t i = 0; i < track.NumberTrajectoryPoints(); ++i) {

            TVector3 point = track.LocationAtPoint(i);
            if ( (point - start).Mag() > 50. ) {
                first_50cm_dir = point - start;
                break;
            }
        } // End loop over trajectory points

        return (end - start).Unit().Dot(first_50cm_dir.Unit());
    } // End ComputeWiggle


}
#endif
