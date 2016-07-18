/**
 * \file EventSelector.h
 *
 * \ingroup CCInclusive
 *
 * \brief Class def header for a class EventSelector
 *
 * @author davidkaleko
 */

/** \addtogroup CCInclusive

    @{*/

#ifndef LARLITE_EVENTSELECTOR_H
#define LARLITE_EVENTSELECTOR_H

#include "Analysis/ana_base.h"
#include "GeoAlgo/GeoVector.h"
#include "GeoAlgo/GeoAABox.h"
#include "DataFormat/track.h"
#include "FidVolBox.h"

namespace larlite {
  /**
     \class EventSelector
     User custom analysis class made by SHELL_USER_NAME
   */
  class EventSelector : public ana_base {

  public:

    /// Default constructor
    EventSelector() { _name = "EventSelector"; _fout = 0;}

    /// Default destructor
    virtual ~EventSelector() {}

    virtual bool initialize();

    virtual bool analyze(storage_manager* storage);

    virtual bool finalize();


    double flashDistZ(const track &longest_track, const double flash_z);

  protected:

    size_t total_events;
    size_t passed_events;

    // Fiducial volume box
    geoalgo::AABox _myGeoAABox;

    bool _flash_in_bgw;
    bool _is_areco_vtx_in_fidvol;
    bool _is_atrack_fromvtx;
    bool _longest_track_nearflash_z;
    bool _longest_trk_contained;
    bool _longest_trk_range_longenough;
    bool _passes_all_cuts;

  };
}
#endif


/** @} */ // end of doxygen group
