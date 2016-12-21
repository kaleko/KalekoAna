#ifndef MCSBIASSTUDY_CXX
#define MCSBIASSTUDY_CXX

#include "MCSBiasStudy.h"
#include "DataFormat/track.h"
#include "DataFormat/mctrack.h"

namespace larlite {

    MCSBiasStudy::MCSBiasStudy() {

        _ana_type = kANALYSIS_TYPE_MAX;

        _myspline = TrackMomentumSplines();
        _tmc = 0;
        _tmc = new kaleko::TrackMomentumCalculator();
        _tmc->SetStepSize(10);


        _tree = new TTree("MCS_bias_tree", "MCS_bias_tree");
        // _tree->Branch("length_analyzed", &_length_analyzed, "length_analyzed/D");
        _tree->Branch("full_length", &_full_length, "full_length/D");
        _tree->Branch("full_integrated_length", &_full_integrated_length, "full_integrated_length/D");
        _tree->Branch("full_range_energy", &_full_range_energy, "full_range_energy/D");
        _tree->Branch("full_range_momentum",&_full_range_momentum, "full_range_momentum/D");
        _tree->Branch("full_integrated_range_energy", &_full_integrated_range_energy, "full_integrated_range_energy/D");
        _tree->Branch("full_integrated_range_momentum",&_full_integrated_range_momentum, "full_integrated_range_momentum/D");
        _tree->Branch("full_MCS_energy", &_full_MCS_energy, "full_MCS_energy/D");
        _tree->Branch("full_MCS_momentum", &_full_MCS_momentum, "full_MCS_momentum/D");
        _tree->Branch("true_E", &_true_E, "true_E/D");
        _tree->Branch("true_momentum", &_true_momentum, "true_momentum/D");
        _tree->Branch("theta",&_theta,"theta/D");
        _tree->Branch("angle_wrt_x",&_angle_wrt_x,"angle_wrt_x/D");
        _tree->Branch("angle_wrt_y",&_angle_wrt_y,"angle_wrt_y/D");
        _tree->Branch("run", &_run, "run/I");
        _tree->Branch("subrun", &_subrun, "subrun/I");
        _tree->Branch("eventid", &_eventid, "eventid/I");
       

        xdir = TVector3(1.,0.,0.);
        ydir = TVector3(0.,1.,0.);
    }


    void MCSBiasStudy::AnalyzeTrack(const larlite::mctrack &mct, int run, int subrun, int eventid, bool flip) {

        if (_ana_type == kANALYSIS_TYPE_MAX) {
            std::cout<< "Did not set what kind of analysis you are doing in MCSBiasStudy!!! Bahhh!" << std::endl;
            return;
        }
        _run = run;
        _subrun = subrun;
        _eventid = eventid;

        _full_length = (mct.back().Position().Vect() - mct.front().Position().Vect()).Mag();
        _full_range_energy = _myspline.GetMuMomentum(_full_length) / 1000. + 0.106;
        _full_range_momentum = std::sqrt(_full_range_energy*_full_range_energy - 0.106*0.106);

        _full_integrated_length = 0;
        for (size_t i = 0; i < mct.size() - 1; ++i)
            _full_integrated_length += (mct.at(i).Position().Vect() - mct.at(i+1).Position().Vect()).Mag();
        
        _full_integrated_range_energy = _myspline.GetMuMomentum(_full_integrated_length) / 1000. + 0.106;
        _full_integrated_range_momentum = std::sqrt(_full_integrated_range_energy*_full_integrated_range_energy - 0.106*0.106);


        _full_MCS_momentum = _tmc->GetMomentumMultiScatterLLHD(mct,_run,_subrun,_eventid);
        _full_MCS_energy = std::sqrt(_full_MCS_momentum*_full_MCS_momentum + 0.106*0.106);
      
        // true E is including mass and is in GEV
        _true_E = mct.front().E() / 1000.;
        _true_momentum = mct.front().Momentum().Vect().Mag() / 1000.;
        _theta = mct.front().Momentum().Vect().Theta();

        _angle_wrt_x = mct.front().Momentum().Vect().Angle(xdir) * (180./3.14159);
        _angle_wrt_y = mct.front().Momentum().Vect().Angle(ydir) * (180./3.14159);

        _tree->Fill();

    }

    void MCSBiasStudy::AnalyzeTrack(const larlite::track &track, int run, int subrun, int eventid, bool flip) {

        if (_ana_type == kANALYSIS_TYPE_MAX) {
            std::cout<< "Did not set what kind of analysis you are doing in MCSBiasStudy!!! Bahhh!" << std::endl;
            return;
        }

        _full_length = -999999.;
        _full_range_energy = -999999.;
        _full_range_momentum = -999999.;
        _full_MCS_energy = -999999.;
        _full_MCS_momentum = -999999.;
        _true_E = -999999.;
        _true_momentum = -999999.;
        _run = run;
        _subrun = subrun;
        _eventid = eventid;

        //GetMomentumMultiScatterLLHD( const larlite::track   &trk, bool flip = false, bool debug = false, bool reweight= false );
        bool apply_reweight = false;
        bool debug_tree = true;

        _full_MCS_momentum = _tmc->GetMomentumMultiScatterLLHD(track, flip, debug_tree, apply_reweight, run, subrun, eventid);
        _full_MCS_energy = std::sqrt(_full_MCS_momentum*_full_MCS_momentum + 0.106*0.106);

        // _length_analyzed = -999.;

        _full_length = (track.End() - track.Vertex()).Mag();
        _full_range_energy = _myspline.GetMuMomentum(_full_length) / 1000. + 0.106;
        _full_range_momentum = std::sqrt(_full_range_energy*_full_range_energy - 0.106*0.106);

         _full_integrated_length = 0;
        for (size_t i = 0; i < track.NumberTrajectoryPoints() - 1; ++i)
            _full_integrated_length += (track.LocationAtPoint(i) - track.LocationAtPoint(i+1)).Mag();
        _full_integrated_range_energy = _myspline.GetMuMomentum(_full_integrated_length) / 1000. + 0.106;
        _full_integrated_range_momentum = std::sqrt(_full_integrated_range_energy*_full_integrated_range_energy - 0.106*0.106);

        _theta = (track.End() - track.Vertex()).Theta();
        _angle_wrt_x = (track.End() - track.Vertex()).Angle(xdir) * (180./3.14159);
        _angle_wrt_y = (track.End() - track.Vertex()).Angle(ydir) * (180./3.14159);

        _tree->Fill();

    }

    // std::pair<double, double> MCSBiasStudy::ComputeWiggle(const larlite::track &track, double seg_size) {

    //     _seg_start_x = -999;
    //     _seg_start_y = -999;
    //     _seg_start_z = -999;
    //     _seg_end_x = -999;
    //     _seg_end_y = -999;
    //     _seg_end_z = -999;
    //     _n_seg_traj_points = -999;
    //     _seg_avg_perp_dist = -999;
    //     _seg_std_perp_dist = -999;

    //     if (track.NumberTrajectoryPoints() < 2) return std::make_pair(-1., -1.);
    //     double seg_size_squared = seg_size * seg_size;//100.; //cm^2
    //     TVector3 start = track.Vertex();
    //     TVector3 seg_end = start;

    //     TVector3 current_seg_dir;
    //     bool found_first_segment = 0;
    //     size_t seg_counter = 0;
    //     std::vector<double> angle_deviations;

    //     // Store segment start indices for second loop later
    //     std::vector<size_t> seg_start_indices;
    //     seg_start_indices.clear();
    //     seg_start_indices.push_back(0);

    //     for (size_t i = 0; i < track.NumberTrajectoryPoints(); ++i) {

    //         TVector3 point = track.LocationAtPoint(i);

    //         // New segment
    //         if ((point - seg_end).Mag2() > seg_size_squared) {
    //             seg_counter++;
    //             seg_start_indices.push_back(i);

    //             if (!found_first_segment) {
    //                 found_first_segment = true;
    //                 current_seg_dir = point - seg_end;
    //                 continue;
    //             }

    //             double this_deviation_dot = current_seg_dir.Unit().Dot((point - seg_end).Unit());
    //             //acos returns NAN for numbers higher than about this
    //             angle_deviations.push_back(this_deviation_dot > 0.999999 ? 0 : std::acos(this_deviation_dot));
    //             current_seg_dir = point - seg_end;
    //             seg_end = point;
    //         } // End new segment

    //     } // End loop over trajectory points

    //     double avg_angle_deviation = 0.;
    //     double std_angle_deviation = 0.;
    //     for (size_t i = 0; i < angle_deviations.size(); ++i)
    //         avg_angle_deviation += angle_deviations.at(i);
    //     avg_angle_deviation /= angle_deviations.size();
    //     for (size_t i = 0; i < angle_deviations.size(); ++i)
    //         std_angle_deviation += std::pow((angle_deviations.at(i) - avg_angle_deviation), 2);
    //     std_angle_deviation /= angle_deviations.size();
    //     std_angle_deviation = std::sqrt(std_angle_deviation);


    //     /// second loop over traj points to compute segment-specific variables for the segment ttree
    //     for (size_t i = 0; i < seg_start_indices.size() - 1; ++i) {
    //         // Build a geo line segment from start point to end point
    //         TVector3 seg_start_point = track.LocationAtPoint(seg_start_indices[i]);
    //         TVector3 seg_end_point   = track.LocationAtPoint(seg_start_indices[i + 1]);
    //         ::geoalgo::LineSegment geo_seg =
    //             ::geoalgo::LineSegment(seg_start_point.X(), seg_start_point.Y(), seg_start_point.Z(),
    //                                    seg_end_point.X(), seg_end_point.Y(), seg_end_point.Z());

    //         _seg_start_x = seg_start_point.X();
    //         _seg_start_y = seg_start_point.Y();
    //         _seg_start_z = seg_start_point.Z();
    //         _seg_end_x = seg_end_point.X();
    //         _seg_end_y = seg_end_point.Y();
    //         _seg_end_z = seg_end_point.Z();

    //         // Loop over trajectory points in the segment
    //         std::vector<double> perp_dists;
    //         for (size_t j = seg_start_indices[i]; j < seg_start_indices[i + 1]; ++j) {
    //             // ::geoalgo::Point geo_pt = ::geoalgo::Point(track.LocationAtPoint(j));
    //             double perp_dist = _geoalg.SqDist(geo_seg, track.LocationAtPoint(j));
    //             perp_dists.push_back(perp_dist);
    //         }
    //         _seg_avg_perp_dist = 0.;
    //         _seg_std_perp_dist = 0.;
    //         for (size_t k = 0; k < perp_dists.size(); ++k)
    //             _seg_avg_perp_dist += perp_dists.at(k);
    //         _seg_avg_perp_dist /= perp_dists.size();
    //         for (size_t k = 0; k < perp_dists.size(); ++k)
    //             _seg_std_perp_dist += std::pow((perp_dists.at(k) - _seg_avg_perp_dist), 2);
    //         _seg_std_perp_dist /= perp_dists.size();
    //         _seg_std_perp_dist = std::sqrt(_seg_std_perp_dist);

    //         _n_seg_traj_points = seg_start_indices[i + 1] - seg_start_indices[i];
    //         _seg_tree->Fill();
    //     }


    //     return seg_counter > 0 ? std::make_pair( avg_angle_deviation, std_angle_deviation ) : std::make_pair(-1., -1.);
    // }// End ComputeWiggle



    // double MCSBiasStudy::ComputeLongCurve(const larlite::track &track) {

    //     if (track.NumberTrajectoryPoints() < 2) return -999.;

    //     TVector3 start = track.Vertex();
    //     TVector3 end = track.End();
    //     TVector3 first_50cm_dir;

    //     if ( (start - end).Mag() < 200. ) return -999.;

    //     for (size_t i = 0; i < track.NumberTrajectoryPoints(); ++i) {

    //         TVector3 point = track.LocationAtPoint(i);
    //         if ( (point - start).Mag() > 50. ) {
    //             first_50cm_dir = point - start;
    //             break;
    //         }
    //     } // End loop over trajectory points

    //     return (end - start).Unit().Dot(first_50cm_dir.Unit());
    // } // End ComputeWiggle


}
#endif
