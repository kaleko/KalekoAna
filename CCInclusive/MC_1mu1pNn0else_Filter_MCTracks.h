/**
 * \file MC_1mu1pNn0else_Filter_MCTracks.h
 *
 * \ingroup EventFilters
 *
 * \brief Class def header for a class MC_1mu1pNn0else_Filter_MCTracks
 *
 * @author davidkaleko
 */

/** \addtogroup EventFilters

    @{*/

#ifndef LARLITE_MC_1MU1PNN0ELSE_FILTER_MCTracks_MCTRACKS_H
#define LARLITE_MC_1MU1PNN0ELSE_FILTER_MCTracks_MCTRACKS_H

#include "Analysis/ana_base.h"
#include "GeoAlgo/GeoAABox.h"
#include "FidVolBox.h"
#include "GeoAlgo/GeoSphere.h"

namespace larlite {
  /**
     \class MC_1mu1pNn0else_Filter_MCTracks
     User custom analysis class made by SHELL_USER_NAME
   */
  class MC_1mu1pNn0else_Filter_MCTracks : public ana_base {

  public:

    /// Default constructor
    MC_1mu1pNn0else_Filter_MCTracks() { _name = "MC_1mu1pNn0else_Filter_MCTracks"; _fout = 0;}

    /// Default destructor
    virtual ~MC_1mu1pNn0else_Filter_MCTracks() {}

    virtual bool initialize();

    virtual bool analyze(storage_manager* storage);

    virtual bool finalize();

  protected:

    // Sphere around neutrino vertex
    geoalgo::Sphere _nu_sphere;
    // Fiducial volume box
    geoalgo::AABox _fidvolBox;
  
    size_t total_events;
    size_t kept_events;
    size_t kaleko_ctr;

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
