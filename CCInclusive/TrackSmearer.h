/**
 * \file TrackSmearer.h
 *
 * \ingroup CCInclusive
 *
 * \brief This class takes a reco track and smears the trajectory points by a set amount,
 *        hopefully getting MC to agree w/ data better by smearing MC
 *
 * @author davidkaleko
 */

/** \addtogroup CCInclusive

    @{*/

#ifndef LARLITE_TRACKSMEARER_H
#define LARLITE_TRACKSMEARER_H

#include "Analysis/ana_base.h"
#include "DataFormat/track.h"
#include "FidVolBox.h"
#include "TRandom3.h"

namespace larlite {
  /**
     \class TrackSmearer
     User custom analysis class made by kaleko
   */
  class TrackSmearer : public larlite_base {

  public:

    /// Default constructor
    TrackSmearer() {
      _box = FidVolBox();
    }

    /// Default destructor
    virtual ~TrackSmearer() {}

    /// Function that makes a new track that is a subset of the input tracks' trajectory points
    /// that are contained inside of the fiducial ("active") volume
    larlite::track SmearTrack(const larlite::track &trk);


  protected:

    FidVolBox _box;

    TRandom3 _trandom;
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
