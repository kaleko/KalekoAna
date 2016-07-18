/**
 * \file TestKinematicNuReco.h
 *
 * \ingroup CCInclusive
 *
 * \brief Class def header for a class TestKinematicNuReco
 *
 * @author davidkaleko
 */

/** \addtogroup CCInclusive

    @{*/

#ifndef LARLITE_TESTKINEMATICNURECO_H
#define LARLITE_TESTKINEMATICNURECO_H

#include "Analysis/ana_base.h"
#include "TTree.h"
#include "TVector3.h"
#include "FidVolBox.h"
#include "NuEnergyCalc.h"
#include "TrackMomentumSplines.h"
#include "TrackMomentumCalculator.h"
#include "GeoAlgo/GeoAABox.h"


namespace larlite {
  /**
     \class TestKinematicNuReco
     User custom analysis class made by SHELL_USER_NAME
   */
  class TestKinematicNuReco : public ana_base {

  public:

    /// Default constructor
    TestKinematicNuReco() { _name = "TestKinematicNuReco"; _fout = 0; _tree = 0; mycalc = 0; myspline = 0; myMCScalc = 0;}

    /// Default destructor
    virtual ~TestKinematicNuReco() {}
    virtual bool initialize();

    virtual bool analyze(storage_manager* storage);
    virtual bool finalize();

  protected:
    void resetTTreeVars();
    TTree* _tree;
    double _true_nu_E;
    double _reco_nu_E_frompp;
    double _plane_Z_anglediff;
    double _true_mu_E;
    double _reco_mu_E_simple;
    double _reco_nu_E_frommuE;
    double _reco_CCQE_E;
    double _reco_mu_E_MCS;
    double _reco_mu_E_range;
    double _reco_mu_E_quadratic;
    bool _mu_contained;
    double _true_mu_len;

    NuEnergyCalc *mycalc;
    TrackMomentumSplines *myspline;
    TrackMomentumCalculator *myMCScalc;

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
