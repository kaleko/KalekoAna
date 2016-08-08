/**
 * \file TrackChopper.h
 *
 * \ingroup CCInclusive
 *
 * \brief This class takes a reco track and chops off all parts of it that are outside of the fiducial volume
 *        by only looking at the portion of tracks inside of the fiducial volume, hopefully data and MC can
 *        agree better, since data has more edge-effects (E-field nonuniformities) than MC
 *
 * @author davidkaleko
 */

/** \addtogroup CCInclusive

    @{*/

#ifndef LARLITE_TRACKCHOPPER_H
#define LARLITE_TRACKCHOPPER_H

#include "Analysis/ana_base.h"
#include "DataFormat/track.h"
#include "FidVolBox.h"

namespace larlite {
  /**
     \class TrackChopper
     User custom analysis class made by kaleko
   */
  class TrackChopper : public larlite_base {

  public:

    /// Default constructor
    TrackChopper() {
      _box = FidVolBox();
    }

    /// Default destructor
    virtual ~TrackChopper() {}

    larlite::track chopTrack(const larlite::track &trk);

  protected:

    FidVolBox _box;

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
