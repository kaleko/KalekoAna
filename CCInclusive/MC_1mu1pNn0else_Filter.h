/**
 * \file MC_1mu1pNn0else_Filter.h
 *
 * \ingroup EventFilters
 *
 * \brief Class def header for a class MC_1mu1pNn0else_Filter
 *
 * @author davidkaleko
 */

/** \addtogroup EventFilters

    @{*/

#ifndef LARLITE_MC_1MU1PNN0ELSE_FILTER_H
#define LARLITE_MC_1MU1PNN0ELSE_FILTER_H

#include "Analysis/ana_base.h"
#include "GeoAlgo/GeoAABox.h"
#include "FidVolBox.h"

namespace larlite {
  /**
     \class MC_1mu1pNn0else_Filter
     User custom analysis class made by SHELL_USER_NAME
   */
  class MC_1mu1pNn0else_Filter : public ana_base {

  public:

    /// Default constructor
    MC_1mu1pNn0else_Filter() { _name = "MC_1mu1pNn0else_Filter"; _fout = 0;}

    /// Default destructor
    virtual ~MC_1mu1pNn0else_Filter() {}

    /** IMPLEMENT in MC_1mu1pNn0else_Filter.cc!
        Initialization method to be called before the analysis event loop.
    */
    virtual bool initialize();

    /** IMPLEMENT in MC_1mu1pNn0else_Filter.cc!
        Analyze a data event-by-event
    */
    virtual bool analyze(storage_manager* storage);

    /** IMPLEMENT in MC_1mu1pNn0else_Filter.cc!
        Finalize method to be called after all events processed.
    */
    virtual bool finalize();

  protected:

    // Fiducial volume box
    geoalgo::AABox _fidvolBox;

    size_t total_events;
    size_t kept_events;

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
