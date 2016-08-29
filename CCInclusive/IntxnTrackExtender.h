/**
 * \file IntxnTrackExtender.h
 *
 * \ingroup CCInclusive
 *
 * \brief Class def header for a class IntxnTrackExtender
 *
 * @author davidkaleko
 */

/** \addtogroup CCInclusive

    @{*/

#ifndef LARLITE_INTXNTRACKEXTENDER_H
#define LARLITE_INTXNTRACKEXTENDER_H

#include <iostream>
#include "Analysis/ana_base.h"
#include "DataFormat/track.h"
#include "DataFormat/vertex.h"
#include "DataFormat/calorimetry.h"
#include "KalekoNuItxn.h"
#include "GeoAlgo/GeoVector.h"
#include "GeoAlgo/GeoSphere.h"
#include "GeoAlgo/GeoAABox.h"
#include "CCInclusiveConstants.h"
#include "KalekoPIDFiller.h"
#include "FidVolBox.h"
#include "GeoAlgo/GeoAlgo.h"

namespace larlite {
  /**
     \class IntxnTrackExtender
     User custom analysis class made by SHELL_USER_NAME
   */
  class IntxnTrackExtender : public ana_base {

  public:

    /// Default constructor
    IntxnTrackExtender() {
      _geoalg = geoalgo::GeoAlgo();
    }

    /// Default destructor
    ~IntxnTrackExtender() {}

    void ExtendVertexTracks( larlite::KalekoNuItxn & itxn, const larlite::event_track * ev_track );

    std::pair<bool, bool> tracksMatched(const larlite::track &base_trk,
                                        bool flip_base,
                                        const larlite::track &match_candidate );

    larlite::track stitchTracks(const larlite::track &base_trk,
                                const larlite::track &match_candidate,
                                bool flip_base,
                                bool flip_match );

  protected:

    geoalgo::GeoAlgo _geoalg;

    double _dotprod;
    double _infdist;
    double _startdist;
    double _enddist;
    double _endprojdist;
    double _startprojdist;
    double _baselen;

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
