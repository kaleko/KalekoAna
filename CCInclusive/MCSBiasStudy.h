/**
 * \file MCSBiasStudy.h
 *
 * \ingroup CCInclusive
 *
 * \brief Class def header for a class MCSBiasStudy
 *
 * @author davidkaleko
 */

/** \addtogroup CCInclusive

    @{*/
#ifndef MCSBIASSTUDY_H
#define MCSBIASSTUDY_H


#include <vector>
#include <iostream>
#include <math.h> //pow
#include "TMath.h"
#include "TVector3.h"
#include "CCInclusiveConstants.h"
#include "TrackMomentumSplines.h"
#include "TrackMomentumCalculator.h"
#include "TTree.h"
#include "TrackChopper.h"
#include "TrackSmearer.h"
#include "GeoAlgo/GeoAlgo.h"

/**
   \class MCSBiasStudy
   A class that analyzes contained tracks in data and MC and fills a ttree with MCS info for calibration purposes.
 */
namespace larlite {

    class MCSBiasStudy {

    public:

        /// Default constructor
        MCSBiasStudy();

        /// Default destructor
        virtual ~MCSBiasStudy() {};

        /// Should be handed a chopped track, and it should be in the correct orienation (no flip needed)
        /// note track chopper does the orientation automatically
        void AnalyzeTrack(const larlite::track &track);

        void AnalyzeTrack(const larlite::mctrack &mct);

        /// Getter for the ttree so the user can write it to a file, for example
        TTree* GetTree() { return _tree; }

        TTree* GetSegTree() { return _seg_tree; }

        TTree* GetTMCTree() { return _tmc->GetTree(); }
        
        /// Function to compute the average angle deviation for custom-segmented tracks
        std::pair<double, double> ComputeWiggle(const larlite::track &track, double seg_size);

        /// Takes the dot product between the first 50cm of the track direction and the start-to-end direction
        double ComputeLongCurve(const larlite::track &track);

    private:

        TrackMomentumSplines _myspline;
        kaleko::TrackMomentumCalculator *_tmc;
        TrackChopper _chopper;
        TrackSmearer _smearer;
        ::geoalgo::GeoAlgo _geoalg;

        TTree *_tree;
        double _length_analyzed;
        double _full_length;
        double _chopped_full_length;
        double _MCS_energy;
        double _full_MCS_energy;
        double _full_MCS_energy_reweighted;
        double _full_MCS_energy_someflipped;
        double _full_MCS_energy_chopped;
        double _full_MCS_energy_smeared;
        double _full_MCS_energy_chopped_smeared;
        double _full_range_energy;
        double _track_start_x;
        double _track_start_y;
        double _track_start_z;
        double _track_end_x;
        double _track_end_y;
        double _track_end_z;
        double _track_dot_z;
        int _n_traj_points;
        double _chopped_wiggle;
        double _chopped_std;
        double _chopped_smeared_wiggle;
        double _chopped_smeared_std;
        double _smeared_wiggle;
        double _smeared_std;
        double _long_curve_dotprod;
        bool _full_track_tree_entry;

        TTree *_seg_tree;
        double _seg_start_x;
        double _seg_start_y;
        double _seg_start_z;
        double _seg_end_x;
        double _seg_end_y;
        double _seg_end_z;
        int    _n_seg_traj_points;
        double _seg_avg_perp_dist;
        double _seg_std_perp_dist;

    };

}// end namespace larlite
#endif
/** @} */ // end of doxygen group

