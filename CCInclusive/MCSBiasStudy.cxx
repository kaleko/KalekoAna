#ifndef MCSBIASSTUDY_CXX
#define MCSBIASSTUDY_CXX

#include "MCSBiasStudy.h"
#include "DataFormat/track.h"

namespace larlite {

  MCSBiasStudy::MCSBiasStudy() {

    _myspline = TrackMomentumSplines();
    _tmc = kaleko::TrackMomentumCalculator();
    _chopper = TrackChopper();

    _tree = new TTree("MCS_bias_tree", "MCS_bias_tree");
    _tree->Branch("length_analyzed", &_length_analyzed, "length_analyzed/D");
    _tree->Branch("full_length", &_full_length, "full_length/D");
    _tree->Branch("chopped_full_length", &_chopped_full_length, "chopped_full_length/D");
    _tree->Branch("full_range_energy", &_full_range_energy, "full_range_energy/D");
    _tree->Branch("MCS_energy", &_MCS_energy, "MCS_energy/D");
    _tree->Branch("full_MCS_energy", &_full_MCS_energy, "full_MCS_energy/D");
    _tree->Branch("full_MCS_energy_someflipped", &_full_MCS_energy_someflipped, "full_MCS_energy_someflipped/D");
    _tree->Branch("full_MCS_energy_chopped", &_full_MCS_energy_chopped, "full_MCS_energy_chopped/D");
    _tree->Branch("track_start_x", &_track_start_x, "track_start_x/D");
    _tree->Branch("track_start_y", &_track_start_y, "track_start_y/D");
    _tree->Branch("track_start_z", &_track_start_z, "track_start_z/D");
    _tree->Branch("track_end_x", &_track_end_x, "track_end_x/D");
    _tree->Branch("track_end_y", &_track_end_y, "track_end_y/D");
    _tree->Branch("track_end_z", &_track_end_z, "track_end_z/D");
    _tree->Branch("track_dot_z", &_track_dot_z, "track_dot_z/D");
    _tree->Branch("n_traj_points", &_n_traj_points, "n_traj_points/I");
    _tree->Branch("chopped_wiggle", &_chopped_wiggle, "chopped_wiggle/D");
    _tree->Branch("chopped_rms", &_chopped_rms, "chopped_rms/D");
    _tree->Branch("full_track_tree_entry", &_full_track_tree_entry, "full_track_tree_entry/O");

  }


  void MCSBiasStudy::AnalyzeTrack(const larlite::track &track) {

    _length_analyzed = -999.;
    _MCS_energy = -999.;

    _full_length = (track.End() - track.Vertex()).Mag();
    _full_range_energy = _myspline.GetMuMomentum(_full_length) / 1000. + 0.106;
    _full_MCS_energy = _tmc.GetMomentumMultiScatterLLHD(track);
    _track_start_x = track.Vertex().X();
    _track_start_y = track.Vertex().Y();
    _track_start_z = track.Vertex().Z();
    _track_end_x = track.End().X();
    _track_end_y = track.End().Y();
    _track_end_z = track.End().Z();
    _full_track_tree_entry = true;
    _n_traj_points = track.NumberTrajectoryPoints();

    bool flip_trk = _track_start_y < _track_end_y;
    _full_MCS_energy_someflipped = _tmc.GetMomentumMultiScatterLLHD(track, flip_trk);

    TVector3 zdir(0., 0., -1.);
    _track_dot_z = !flip_trk ? (track.End() - track.Vertex()).Unit().Dot(zdir) :
                   (track.Vertex() - track.End()).Unit().Dot(zdir);

    auto const chopped_trk = _chopper.chopTrack(track);
    _chopped_full_length = (chopped_trk.End() - chopped_trk.Vertex()).Mag();
    _full_MCS_energy_chopped = _tmc.GetMomentumMultiScatterLLHD(chopped_trk);

    auto wiggle = ComputeWiggle(chopped_trk);
    _chopped_wiggle = wiggle.first;
    _chopped_rms = wiggle.second;

    _tree->Fill();


    _full_track_tree_entry = false;


    // // slowly build up a copy of the initial track with increasing length
    // // and compute MCS energy as length increases
    // larlite::track dummy_trk;

    // double length_increment = 10.; //cm

    // double current_min_len = 100.;
    // size_t n_points = track.NumberTrajectoryPoints();

    // for (int i = 1; i < n_points; ++i) {

    //   dummy_trk.add_vertex(track.LocationAtPoint(i - 1));
    //   dummy_trk.add_direction(track.LocationAtPoint(i - 1));

    //   _length_analyzed = (track.LocationAtPoint(i) - track.LocationAtPoint(0)).Mag();
    //   if ( _length_analyzed < current_min_len ) continue;

    //   _MCS_energy = _tmc.GetMomentumMultiScatterLLHD(dummy_trk);

    //   current_min_len += length_increment;

    //   _tree->Fill();

    // }

  }

  std::pair<double,double> MCSBiasStudy::ComputeWiggle(const larlite::track &track) {

    double seg_size_squared = 100.; //cm^2
    TVector3 start = track.Vertex();
    TVector3 seg_end = start;
    
    TVector3 current_seg_dir;
    bool found_first_segment = 0;
    size_t seg_counter = 0;
    std::vector<double> angle_deviations;
    for (size_t i = 0; i < track.NumberTrajectoryPoints(); ++i) {

      TVector3 point = track.LocationAtPoint(i);

      // New segment
      if ((point - seg_end).Mag2() > seg_size_squared) {
        seg_counter++;
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
    double rms_angle_deviation = 0.;
    for (size_t i = 0; i < angle_deviations.size(); ++i){
      avg_angle_deviation += angle_deviations.at(i);
      rms_angle_deviation += (angle_deviations.at(i) * angle_deviations.at(i));
    }
    avg_angle_deviation /= angle_deviations.size();
    rms_angle_deviation /= angle_deviations.size();
    rms_angle_deviation = std::sqrt(rms_angle_deviation);

    return seg_counter > 0 ? std::make_pair( avg_angle_deviation, rms_angle_deviation ) : std::make_pair(-1.,-1.);
  } // End ComputeWiggle

}
#endif
