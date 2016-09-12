#ifndef MCSBIASSTUDY_CXX
#define MCSBIASSTUDY_CXX

#include "MCSBiasStudy.h"
#include "DataFormat/track.h"

namespace larlite {

  MCSBiasStudy::MCSBiasStudy() {

    _myspline = TrackMomentumSplines();
    _tmc = kaleko::TrackMomentumCalculator();

    _tree = new TTree("MCS_bias_tree", "MCS_bias_tree");
    _tree->Branch("length_analyzed", &_length_analyzed, "length_analyzed/D");
    _tree->Branch("full_length", &_full_length, "full_length/D");
    _tree->Branch("full_range_energy", &_full_range_energy, "full_range_energy/D");
    _tree->Branch("MCS_energy", &_MCS_energy, "MCS_energy/D");
    _tree->Branch("full_MCS_energy", &_full_MCS_energy, "full_MCS_energy/D");
    _tree->Branch("track_start_x", &_track_start_x, "track_start_x/D");
    _tree->Branch("track_start_y", &_track_start_y, "track_start_y/D");
    _tree->Branch("track_start_z", &_track_start_z, "track_start_z/D");
    _tree->Branch("track_end_x", &_track_end_x, "track_end_x/D");
    _tree->Branch("track_end_y", &_track_end_y, "track_end_y/D");
    _tree->Branch("track_end_z", &_track_end_z, "track_end_z/D");
    _tree->Branch("full_track_tree_entry", &_full_track_tree_entry, "full_track_tree_entry/O");

  }


  void MCSBiasStudy::AnalyzeTrack(const larlite::track &track) {

    _length_analyzed = -999.;
    _MCS_energy = -999.;

    _full_length = (track.End() - track.Vertex()).Mag();
    _full_range_energy = _myspline.GetMuMomentum(_full_length) / 1000. + 0.106;
    _full_MCS_energy = _MCS_energy = _tmc.GetMomentumMultiScatterLLHD(track);
    _track_start_x = track.Vertex().X();
    _track_start_y = track.Vertex().Y();
    _track_start_z = track.Vertex().Z();
    _track_end_x = track.End().X();
    _track_end_y = track.End().Y();
    _track_end_z = track.End().Z();
    _full_track_tree_entry = true;
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



}

#endif
