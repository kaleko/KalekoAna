/**
 * \file EfficiencyStudy.h
 *
 * \ingroup CCInclusive
 *
 * \brief Class def header for a class EfficiencyStudy
 *
 * @author davidkaleko
 */

/** \addtogroup CCInclusive

    @{*/

#ifndef LARLITE_EFFICIENCYSTUDY_H
#define LARLITE_EFFICIENCYSTUDY_H

#include "Analysis/ana_base.h"
#include "TTree.h"
#include "GeoAlgo/GeoAlgo.h"
#include "GeoAlgo/GeoVector.h"
#include "GeoAlgo/GeoAABox.h"
#include "DataFormat/track.h"
#include "TrackMomentumSplines.h"
#include "NuEnergyCalc.h"
#include "FidVolBox.h"

namespace larlite {
  /**
     \class EfficiencyStudy
     User custom analysis class made by SHELL_USER_NAME
   */
  class EfficiencyStudy : public ana_base {

  public:

    /// Default constructor
    EfficiencyStudy() { 
      _name = "EfficiencyStudy"; 
      _fout = 0; 
      _vtx_tree = 0; 
      _evt_tree = 0;
      _flash_tree = 0;
    }

    /// Default destructor
    virtual ~EfficiencyStudy() {}

    virtual bool initialize();

    virtual bool analyze(storage_manager* storage);

    virtual bool finalize();

    double flashDistZ(const track &longest_track, const double flash_z);

  protected:

    // This tree is filled once per reconstructed neutrino vertex
    TTree *_vtx_tree;
    double _dist_vtx_truth;

    // This tree is filled once per art event
    TTree *_evt_tree;
    bool _is_truth_fiducial;
    bool _is_numuCC;
    bool _flash_in_bgw;
    bool _is_areco_vtx_in_fidvol;
    bool _is_atrack_fromvtx;
    bool _longest_track_nearflash_z;
    bool _longest_trk_contained;
    bool _longest_trk_range_longenough;
    bool _longest_trk_start_near_truth_nu;
    double _dist_trk_flash;
    bool _passes_all_cuts;


    double _true_nu_E;
    double _ccqe_nu_E;
    double _true_mu_E;
    double _reco_mu_E;

    
    // This tree is to debug flash timing and is filled once per flash
    TTree *_flash_tree;
    double _flash_time; //in us

    //Geometry algo instances
    geoalgo::GeoAlgo _geoAlgo;
    // Fiducial volume box
    geoalgo::AABox _myGeoAABox;

    TrackMomentumSplines *myspline;
    NuEnergyCalc *mycalc;
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
