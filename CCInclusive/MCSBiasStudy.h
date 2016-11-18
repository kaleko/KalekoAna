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
   A class that analyzes contained tracks in data and MC and fills a ttree with MCS info for calibration/analysis purposes.
 */
namespace larlite {

    class MCSBiasStudy {

    public:

        /// Default constructor
        MCSBiasStudy();

        // See .cxx file for descriptions of each of these types
        enum AnalysisType_t {
            kSingleMuonMCTrack, 
            kSingleMuonRecoTrack,
            kMCBNBSelectedRecoTrack,
            kMCBNBRecoTrack,
            kDataBNBSelectedRecoTrack,
            kANALYSIS_TYPE_MAX
        };

        /// Default destructor
        virtual ~MCSBiasStudy() {};

        void AnalyzeTrack(const larlite::track &track, int run, int subrun, int eventid, bool flip = false);

        void AnalyzeTrack(const larlite::mctrack &mct, int run, int subrun, int eventid, bool flip = false);

        /// Getter for the ttree so the user can write it to a file, for example
        TTree* GetTree() { return _tree; }

        // TTree* GetSegTree() { return _seg_tree; }

        TTree* GetTMCTree() { return _tmc->GetTree(); }

        /// Function to compute the average angle deviation for custom-segmented tracks
        // std::pair<double, double> ComputeWiggle(const larlite::track &track, double seg_size);

        /// Takes the dot product between the first 50cm of the track direction and the start-to-end direction
        // double ComputeLongCurve(const larlite::track &track);

        void SetAnalysisType(AnalysisType_t mytype) { _ana_type = mytype; }

    private:

        

        AnalysisType_t _ana_type;

        TrackMomentumSplines _myspline;
        kaleko::TrackMomentumCalculator *_tmc;
        // ::geoalgo::GeoAlgo _geoalg;

        TTree *_tree;
        // double _length_analyzed;
        double _full_length;
        double _full_MCS_energy;
        double _full_MCS_momentum;
        double _full_range_energy;
        double _full_range_momentum;
        double _true_E;
        int _run;
        int _subrun;
        int _eventid;

    };

}// end namespace larlite
#endif
/** @} */ // end of doxygen group

