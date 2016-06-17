/**
 * \file QuickMikePlot.h
 *
 * \ingroup CCInclusive
 *
 * \brief Class def header for a class QuickMikePlot
 *
 * @author davidkaleko
 */

/** \addtogroup CCInclusive

    @{*/

#ifndef LARLITE_QUICKMIKEPLOT_H
#define LARLITE_QUICKMIKEPLOT_H

#include "Analysis/ana_base.h"
#include "TTree.h"
#include "GeoAlgo/GeoAABox.h"
#include "LArUtil/Geometry.h"
#include "GeoAlgo/GeoSphere.h"

namespace larlite {
  /**
     \class QuickMikePlot
     User custom analysis class made by SHELL_USER_NAME
   */
  class QuickMikePlot : public ana_base {

  public:

    /// Default constructor
    QuickMikePlot() { _name = "QuickMikePlot"; _fout = 0; _tree = 0;}

    /// Default destructor
    virtual ~QuickMikePlot() {}
    virtual bool initialize();
    virtual bool analyze(storage_manager* storage);
    virtual bool finalize();


  protected:
    TTree *_tree;
    bool _CCNC;
    int _nupdg;
    int _fndecay;
    double _true_nu_E;
    double _true_mu_E;
    double _true_mu_len;
    bool _nu_infidvol;



    // Sphere around neutrino vertex
    geoalgo::Sphere _nu_sphere;

    // Fiducial volume box
    geoalgo::AABox _fidvolBox;
    double fidvol_dist;
    double fidvol_dist_y;


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
