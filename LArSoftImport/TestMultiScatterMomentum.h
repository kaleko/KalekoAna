/**
 * \file TestMultiScatterMomentum.h
 *
 * \ingroup LArSoftImport
 *
 * \brief Class def header for a class TestMultiScatterMomentum
 *
 * @author davidkaleko
 */

/** \addtogroup LArSoftImport

    @{*/

#ifndef LARLITE_TESTMULTISCATTERMOMENTUM_H
#define LARLITE_TESTMULTISCATTERMOMENTUM_H

#include "Analysis/ana_base.h"
#include "TrackMomentumCalculator.h"
#include "LArUtil/Geometry.h"
#include "GeoAlgo/GeoAABox.h"

namespace larlite {
  /**
     \class TestMultiScatterMomentum
     User custom analysis class made by SHELL_USER_NAME
   */
  class TestMultiScatterMomentum : public ana_base {

  public:

    /// Default constructor
    TestMultiScatterMomentum() { _name = "TestMultiScatterMomentum"; _fout = 0; _ana_tree = 0; _using_mctracks = true;}

    /// Default destructor
    virtual ~TestMultiScatterMomentum() {}

    virtual bool initialize();

    virtual bool analyze(storage_manager* storage);

    virtual bool finalize();

    void SetUseMCTracks(bool flag) { _using_mctracks = flag; }

  protected:

    bool _using_mctracks;

    kaleko::TrackMomentumCalculator _tmc;

    TTree* _ana_tree;
    double _true_mom;
    double _mcs_reco_mom;
    double _mcs_reco_mom_backwards;
    double _true_length;
    double _reco_length;
    bool _mu_contained;

    geoalgo::AABox _fidvolBox;

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
