/**
 * \file FindSimilarTrack.h
 *
 * \ingroup CCInclusive
 *
 * \brief This class takes a reco track by producer A, then a list of tracks by producer B,
 *        and returns the track by producer B that most closely matches the track by producer A.
 *        It does this by picking the track by producer B that has the most similar direction 
 *        and length to the track by producer A. It also requires that the track by producer B
 *        starts (or ends) within 5cm of the track by producer A.
 *
 * @author davidkaleko
 */

/** \addtogroup CCInclusive

    @{*/

#ifndef LARLITE_FINDSIMILARTRACK_H
#define LARLITE_FINDSIMILARTRACK_H

#include "Analysis/ana_base.h"
#include "DataFormat/track.h"

namespace larlite {
  /**
     \class FindSimilarTrack
     User custom analysis class made by kaleko
   */
  class FindSimilarTrack : public larlite_base {

  public:

    /// Default constructor
    FindSimilarTrack() {}

    /// Default destructor
    virtual ~FindSimilarTrack() {}

    larlite::track findSimilarTrack(const larlite::track &trk, const larlite::event_track &ev_track);


  protected:


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
