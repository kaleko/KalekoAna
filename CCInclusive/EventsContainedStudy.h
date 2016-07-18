/**
 * \file EventsContainedStudy.h
 *
 * \ingroup CCInclusive
 *
 * \brief Class def header for a class EventsContainedStudy
 *
 * @author davidkaleko
 */

/** \addtogroup CCInclusive

    @{*/

#ifndef LARLITE_EVENTSCONTAINEDSTUDY_H
#define LARLITE_EVENTSCONTAINEDSTUDY_H

#include "Analysis/ana_base.h"
#include "FidVolBox.h"
#include "TH1.h"
#include "GeoAlgo/GeoAABox.h"
#include "GeoAlgo/GeoSphere.h"
#include "TrackMomentumSplines.h"
#include "NuEnergyCalc.h"

namespace larlite {
  /**
     \class EventsContainedStudy
     User custom analysis class made by SHELL_USER_NAME
   */
  class EventsContainedStudy : public ana_base {

  public:

    /// Default constructor
    EventsContainedStudy() {
      _name = "EventsContainedStudy";
      _fout = 0;
      _tree = 0;
      _h_mu_E_contained = 0;
      _h_mu_E_all = 0;
      _h_mu_E_contained = 0;
      _h_mu_E_all = 0;
      myspline = 0;
      mycalc = 0;
    }

    /// Default destructor
    virtual ~EventsContainedStudy() {}

    virtual bool initialize();

    virtual bool analyze(storage_manager* storage);

    virtual bool finalize();

  protected:

    // Sphere around neutrino vertex
    geoalgo::Sphere _nu_sphere;

    TTree *_tree;
    bool _mu_contained;
    bool _p_contained;
    size_t _n_reco_tracks_fromvtx;
    size_t _n_reco_tracks_endingnearvtx;
    size_t _n_mc_tracks_fromvtx;
    double _true_mu_E;
    double _true_p_E;
    double _reco_mu_E;
    double _reco_p_E;
    double _true_nu_E;
    double _true_nu_x;
    double _true_nu_y;
    double _true_nu_z;
    double _reco_CCQE_E;
    double _true_mu_len;
    double _true_p_len;
    double _reco_mu_len;
    double _reco_p_len;
    double _true_neutrons_E;
    bool _is_reco_vtx_near_true_vtx;
    int runno;
    int subrunno;
    int evtno;
      // GNUMI decay codes (as in Ndecay)
    /* --------------------------
       Decay type Decay code
       K0L -> nue pi- e+  1
       K0L -> nuebar pi+ e- 2
       K0L -> numu pi- mu+  3
       K0L -> numubar pi+ mu- 4
       K+ -> numu mu+         5
       K+ -> nue pi0 e+         6
       K+ -> numu pi0 mu+ 7
       K- -> numubar mu-  8
       K- -> nuebar pi0 e-  9
       K- -> numubar pi0 mu-  10
       mu+ -> numubar nue e+  11
       mu- -> numu nuebar e-  12
       pi+ -> numu mu+          13
       pi- -> numubar mu- 14
       -------------------------- */
    int _fndecay;


    geoalgo::AABox _fidvolBox;

    size_t _evts_pass_filter;
    size_t _evts_fully_contained;
    size_t kaleko_ctr;

    TrackMomentumSplines *myspline;
    NuEnergyCalc *mycalc;


    TH1F *_h_mu_E_contained;
    TH1F *_h_mu_E_all;
    TH1F *_h_nu_E_contained;
    TH1F *_h_nu_E_all;
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
