/**
 * \file XiaoEventAna.h
 *
 * \ingroup CCInclusive
 *
 * \brief Ana module that uses XiaoNuFinder to find 1mu1p events then makes histograms about them
 *
 *     track: pandoraNuPMA
 *     vtx: pmtrack
 *     calorimetry: pandoraNuPMAcalo
 *
 *
 * @author davidkaleko
 */

/** \addtogroup CCInclusive

    @{*/

#ifndef LARLITE_XIAOEVENTANA_H
#define LARLITE_XIAOEVENTANA_H

#include "Analysis/ana_base.h"
#include "DataFormat/track.h"
#include "DataFormat/vertex.h"
#include "DataFormat/calorimetry.h"
#include "KalekoNuItxn.h"
#include "GeoAlgo/GeoVector.h"
#include "GeoAlgo/GeoSphere.h"
#include "GeoAlgo/GeoAABox.h"
#include "XiaoNuFinder.h"
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "CCInclusiveConstants.h"
#include "TrackMomentumSplines.h"
#include "NuEnergyCalc.h"
#include "FidVolBox.h"
#include "TrackMomentumCalculator.h"
#include "IntxnBooster.h"
#include "KalekoPIDFiller.h"
#include "TrackChopper.h"
#include "TrackSmearer.h"
#include "IntxnTrackExtender.h"
#include "MCSBiasStudy.h"
#include "FindSimilarTrack.h"


namespace larlite {
    /**
       \class XiaoEventAna
       User custom analysis class made by SHELL_USER_NAME
     */

    class XiaoEventAna : public ana_base {

    public:

        /// Default constructor
        XiaoEventAna() {
            _name = "XiaoEventAna";
            _fout = 0;
            // _hmult = 0;
            // _hdedx = 0;
            _hcorrect_ID = 0;
            _running_on_data = false;
            _tree = 0;
            _filetype = kINPUT_FILE_TYPE_MAX;
            _vtx_sphere_radius = 3.;
            _track_producer = "pandoraNuPMA";
            _vtx_producer = "pmtrack";
            _calo_producer = _track_producer + "calo";
            _extend_track_producer = "";
            _match_track_producer = "";
            _min_trk_len = 0.;
            _mcs_min_trk_len = 100.;
            _mcsbiasstudy = 0;
            _MCS_seg_size = 10.;
            _MCScalc = 0;

        }

        /// Default destructor
        virtual ~XiaoEventAna() {}

        virtual bool initialize();

        virtual bool analyze(storage_manager* storage);

        virtual bool finalize();

        void setRunningOnData(bool flag) { _running_on_data = flag; }

        void setInputType(InputFileType_t filetype) { _filetype = filetype; }

        void setVtxSphereRadius(double myradius) { _vtx_sphere_radius = myradius; }

        bool setBGWTimes();

        void setTrackProducer(std::string prod) { _track_producer = prod; }

        void setExtendTrackProducer(std::string prod) { _extend_track_producer = prod; }

        void setMatchTrackProducer(std::string prod) { _match_track_producer = prod; }

        void setVtxProducer(std::string prod) { _vtx_producer = prod; }

        void setCaloProducer(std::string prod) { _calo_producer = prod; }

        // This is minimum length to SELECT THE EVENTS
        void setMinTrkLen(double len) { _min_trk_len = len; }

        // This is minimum length the MCS scattering code uses
        void setMCSMinLen(double len) { _mcs_min_trk_len = len; }

        void setMCSSegSize(double len) { _MCS_seg_size = len; }


    protected:

        std::string _track_producer;
        std::string _extend_track_producer;
        std::string _match_track_producer;
        std::string _vtx_producer;
        std::string _calo_producer;

        // Minimum track length to consider a track anywhere in analysis. Default 0.
        double _min_trk_len;

        // Minimum track length for the MCS code to use
        double _mcs_min_trk_len;

        double _MCS_seg_size;

        XiaoNuFinder _nu_finder;
        TrackMomentumSplines _myspline;
        kaleko::TrackMomentumCalculator *_MCScalc;
        NuEnergyCalc _nu_E_calc;
        IntxnBooster _intxn_booster;
        KalekoPIDFiller _PID_filler;
        TrackChopper _chopper;
        IntxnTrackExtender _trkextender;
        MCSBiasStudy *_mcsbiasstudy;
        TrackSmearer _smearer;
        FindSimilarTrack _trackmatcher;

        double _vtx_sphere_radius;

        void resetTTreeVars();

        InputFileType_t _filetype;

        size_t total_events;
        size_t passed_events;

        double BGW_mintime;
        double BGW_maxtime;

        // // Fiducial volume box
        geoalgo::AABox _fidvolBox;

        geoalgo::AABox _tpcBox;

        // TH1F *_hmult;
        // TH2D *_hdedx;
        TH1F *_hcorrect_ID;

        bool _running_on_data;

        TTree *_tree;
        double _mu_start_dedx;
        double _mu_end_dedx;
        bool   _correct_ID;
        double _mu_phi;
        double _p_phi;
        bool _longest_trk_contained;
        bool _all_trks_contained;
        double _true_nu_E;
        int    _true_nu_pdg;
        bool   _true_nu_CCNC;
        int    _true_nu_mode;
        int    _true_lepton_pdg;
        double _true_lepton_momentum;
        int    _fndecay;
        double _fppdxdz;
        double _fppdydz;
        double _fpppz;
        double _fppenergy;
        double _kaon_prod_E;
        double _kaon_prod_theta;
        double _mu_p_dirdot;
        int _n_associated_tracks;

        double _longest_trk_len;
        double _longest_trk_len_infidvol;
        double _longest_trk_Length;
        double _longest_trk_Length_infidvol;
        double _longest_track_end_x_infidvol;
        double _longest_track_end_y_infidvol;
        double _longest_track_end_z_infidvol;

        double _longest_trk_cosy;
        double _second_longest_trk_len;
        double _longest_trk_theta;
        double _longest_trk_phi;
        double _longest_trk_MCS_mom;
        double _longest_trk_MCS_mom_chopped;
        double _longest_trk_MCS_mom_reweighted;
        double _matched_longest_trk_MCS_mom;
        double _matched_longest_trk_MCS_mom_chopped;

        double _longest_trk_MCS_mom_chopped_straightened;
        double _longest_trk_MCS_mom_smeared_MC;
        double _longest_trk_spline_mom;
        double _nu_E_estimate;
        double _corrected_nu_E_estimate;
        double _CCQE_E;
        double _longest_trk_avg_calo;
        double _second_longest_trk_avg_calo;
        double _longest_trk_wiggle;
        double _longest_trk_wiggle_5cm_seglen;
        double _longest_trk_wiggle_10cm_seglen;
        double _longest_trk_wiggle_15cm_seglen;
        double _longest_trk_wiggle_20cm_seglen;
        double _longest_trk_smeared_wiggle;

        double _true_nu_x;
        double _true_nu_y;
        double _true_nu_z;
        double _reco_nu_x;
        double _reco_nu_y;
        double _reco_nu_z;

        double _dist_reco_true_vtx;
        double _max_tracks_dotprod;
        double _longest_tracks_dotprod;
        double _longest_tracks_dotprod_trkendpoints;

        double _longest_track_end_x;
        double _longest_track_end_y;
        double _longest_track_end_z;

        double _brightest_BSW_flash_PE;
        double _BSW_flash_z_range;

        double _dist_vtx_to_TPC;

        double _longest_trk_dot_truemuondir;

        int _n_reco_nu_in_evt;
        double _E_lepton;
        double _E_hadrons;
        double _E_MCS;
        double _E_range;

        double _smallest_avg_calo;

        bool _longest_track_well_recod;
        bool _found_mu_mctrack;
        int _true_multiplicity;
        double _tot_E_mcshowers;
        double _tot_E_mctracks;

        double _longest_trk_Length_infidvol_smeared;
        double _longest_trk_len_infidvol_smeared;
    };
}
#endif

/** @} */ // end of doxygen group
