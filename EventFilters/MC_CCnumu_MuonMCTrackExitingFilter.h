/**
 * \file MC_CCnumu_MuonMCTrackExitingFilter.h
 *
 * \ingroup EventFilters
 *
 * \brief Class def header for a class MC_CCnumu_MuonMCTrackExitingFilter
 *
 * @author davidkaleko
 */

/** \addtogroup EventFilters

    @{*/

#ifndef LARLITE_MC_CCNUMU_MUONMCTRACKEXITINGFILTER_H
#define LARLITE_MC_CCNUMU_MUONMCTRACKEXITINGFILTER_H

#include "Analysis/ana_base.h"
#include "DataFormat/mctruth.h"
#include "GeoAlgo/GeoAABox.h"
#include "LArUtil/Geometry.h"
#include "DataFormat/mctrack.h"
#include "FidVolBox.h"

namespace larlite {
  /**
     \class MC_CCnumu_MuonMCTrackExitingFilter
     User custom analysis class made by SHELL_USER_NAME
   */
  class MC_CCnumu_MuonMCTrackExitingFilter : public ana_base {

  public:

    /// Default constructor
    MC_CCnumu_MuonMCTrackExitingFilter() { _name = "MC_CCnumu_MuonMCTrackExitingFilter"; _fout = 0;}

    /// Default destructor
    virtual ~MC_CCnumu_MuonMCTrackExitingFilter() {}

    /** IMPLEMENT in MC_CCnumu_MuonMCTrackExitingFilter.cc!
        Initialization method to be called before the analysis event loop.
    */
    virtual bool initialize();

    /** IMPLEMENT in MC_CCnumu_MuonMCTrackExitingFilter.cc!
        Analyze a data event-by-event
    */
    virtual bool analyze(storage_manager* storage);

    /** IMPLEMENT in MC_CCnumu_MuonMCTrackExitingFilter.cc!
        Finalize method to be called after all events processed.
    */
    virtual bool finalize();

  protected:

    FidVolBox _fidvol;

    std::vector<size_t> _counters;

    size_t _n_total_events;
    size_t _n_kept_events;

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
