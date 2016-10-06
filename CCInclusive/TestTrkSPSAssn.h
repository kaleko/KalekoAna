/**
 * \file TestTrkSPSAssn.h
 *
 * \ingroup CCInclusive
 *
 * \brief Class def header for a class TestTrkSPSAssn
 *
 * @author davidkaleko
 */

/** \addtogroup CCInclusive

    @{*/

#ifndef LARLITE_TESTTRKSPSASSN_H
#define LARLITE_TESTTRKSPSASSN_H

#include "Analysis/ana_base.h"
#include "TTree.h"
#include "FidVolBox.h"
#include "TrackMomentumCalculator.h"

namespace larlite {

  class TestTrkSPSAssn : public ana_base {

  public:

    /// Default constructor
    TestTrkSPSAssn() { _name = "TestTrkSPSAssn"; _fout = 0; _tree = 0; _tmc = 0; _track_producer = ""; _sps1track0 = 0; _seglen = 10.;}

    /// Default destructor
    virtual ~TestTrkSPSAssn() {}

    virtual bool initialize();

    virtual bool analyze(storage_manager* storage);

    virtual bool finalize();

    void SetTrackProducer(std::string prod) { _track_producer = prod; }
    void SetSPS1Track0(int flag) { _sps1track0 = flag; }
    void SetSeglen(double seglen) { _seglen = seglen; }

  protected:

    std::string _track_producer;
    int _sps1track0;
    double _seglen;

    TTree *_tree;
    size_t _n_total_sps;
    size_t _n_asstd_sps;
    size_t _n_track_traj_pts;
    double _dist_firstsps_trackstart;
    double _dist_firstsps_trackend;
    double _track_true_E;
    double _mct_MCS_E;
    double _track_MCS_E;
    double _sps_MCS_E;
    double _trk_len;
    double _sps_trk_len;

    FidVolBox _fidvol;
    kaleko::TrackMomentumCalculator *_tmc;

  };
}
#endif

//**************************************************************************
//
// For Analysis framework documentation, read Manual.pdf here:
//
// http://microboone-docdb.fnal.gov:8080/cgi-bin/ShowDocument?docid=3183
//
//**************************************************************************

/** @} */ // end of doxygen group
