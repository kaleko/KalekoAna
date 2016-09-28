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

namespace larlite {

  class TestTrkSPSAssn : public ana_base {

  public:

    /// Default constructor
    TestTrkSPSAssn() { _name = "TestTrkSPSAssn"; _fout = 0; _tree = 0;}

    /// Default destructor
    virtual ~TestTrkSPSAssn() {}

    virtual bool initialize();

    virtual bool analyze(storage_manager* storage);

    virtual bool finalize();

  protected:

    TTree *_tree;
    size_t _n_total_sps;
    size_t _n_asstd_sps;
    size_t _n_track_traj_pts;
    double _dist_firstsps_trackstart;
    double _dist_firstsps_trackend;

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
