/**
 * \file KalekoTrackStitcher.h
 *
 * \ingroup CCInclusive
 *
 * \brief Class def header for a class KalekoTrackStitcher
 *
 * @author davidkaleko
 */

/** \addtogroup CCInclusive

    @{*/

#ifndef LARLITE_KALEKOTRACKSTITCHER_H
#define LARLITE_KALEKOTRACKSTITCHER_H

#include "Analysis/ana_base.h"
#include "DataFormat/track.h"
#include "TTree.h"
#include "GeoAlgo/GeoAlgo.h"

namespace larlite {
  /**
     \class KalekoTrackStitcher
     User custom analysis class made by kaleko
   */
  class KalekoTrackStitcher : public ana_base {

  public:

    /// Default constructor
    KalekoTrackStitcher() {
      _name = "KalekoTrackStitcher";
      _fout = 0;
      _base_producer = "";
      _match_producer = "";
      _debug_tree = 0;
    }

    /// Default destructor
    virtual ~KalekoTrackStitcher() {}

    virtual bool initialize();

    virtual bool analyze(storage_manager* storage);

    virtual bool finalize();

    void setBaseProducer(std::string producer) { _base_producer = producer; }

    void setMatchProducer(std::string producer) { _match_producer = producer; }

  protected:

    bool tracksMatched(const larlite::track &base_trk, const larlite::track &match_candidate );

    larlite::track stitchTracks(const larlite::track &base_trk,
                                const larlite::track &match_candidate );

    bool trackViable(const larlite::track &trk);

    std::string _base_producer;
    std::string _match_producer;

    geoalgo::GeoAlgo _geoalg;

    TTree *_debug_tree;
    double _dotprod;
    double _infdist;
    double _startdist;
    double _enddist;

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
