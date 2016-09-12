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

        /// Getter for the ttree so the user can write it to a file, for example
        TTree* GetTree() { return _tree; }

    private:

        TrackMomentumSplines _myspline;
        kaleko::TrackMomentumCalculator _tmc;
        TrackChopper _chopper;

        TTree *_tree;
        double _length_analyzed;
        double _full_length;
        double _chopped_full_length;
        double _MCS_energy;
        double _full_MCS_energy;
        double _full_MCS_energy_someflipped;
        double _full_MCS_energy_chopped;
        double _full_range_energy;
        double _track_start_x;
        double _track_start_y;
        double _track_start_z;
        double _track_end_x;
        double _track_end_y;
        double _track_end_z;
        double _track_dot_z;
        int _n_traj_points;
        bool _full_track_tree_entry;

    };

}// end namespace larlite
#endif
/** @} */ // end of doxygen group

