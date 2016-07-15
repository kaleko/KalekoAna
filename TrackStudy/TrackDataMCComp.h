/**
 * \file TrackDataMCComp.h
 *
 * \ingroup TrackStudy
 *
 * \brief Class def header for a class TrackDataMCComp
 *
 * @author davidkaleko
 */

/** \addtogroup TrackStudy

    @{*/

#ifndef LARLITE_TRACKDATAMCCOMP_H
#define LARLITE_TRACKDATAMCCOMP_H

#include "Analysis/ana_base.h"
#include "TTree.h"
#include <string>

namespace larlite {
  /**
     \class TrackDataMCComp
     User custom analysis class made by SHELL_USER_NAME
   */
  class TrackDataMCComp : public ana_base {

  public:

    /// Default constructor
    TrackDataMCComp() { 
      _name = "TrackDataMCComp";
     _fout = 0;
     _track_producer = "";
     _ana_tree = 0;
     _mctrack_tree = 0;
     _track_tree = 0;
     _running_on_data = false;
   }

    /// Default destructor
    virtual ~TrackDataMCComp() {}

    virtual bool initialize();

    virtual bool analyze(storage_manager* storage);

    virtual bool finalize();

    void setTrackProducer(std::string producer) { _track_producer = producer; }

    void setRunningOnData(bool flag) { _running_on_data = flag; }

  protected:

    // Track producer name. Must be set.
    std::string _track_producer;

    // Whether to look for/analyze MC stuff
    bool _running_on_data;

    // This tree is filled once per event
    TTree *_ana_tree;
    // # of MCTracks of cosmic origin with at least 3cm in the detector.
    size_t _n_mctracks_cosmic;
    // # of MCTracks of neutrino origin with at least 3cm in the detector.
    size_t _n_mctracks_neutrino;
    // # of reconstructed tracks by _track_producer
    size_t _n_recotracks;

    // This tree is filled once per reco track that is at least 3cm long
    TTree *_track_tree;
    // Length of track
    double _trk_len;
    // Start vtx of track
    double _trk_start_x;
    double _trk_start_y;
    double _trk_start_z;
    // End vtx of track
    double _trk_end_x;
    double _trk_end_y;
    double _trk_end_z;

    // This tree is filled once per mctrack
    TTree *_mctrack_tree;
    // Length of track
    double _mctrk_len;
    // Start vtx of mctrack
    double _mctrk_start_x;
    double _mctrk_start_y;
    double _mctrk_start_z;
    // End vtx of mctrack
    double _mctrk_end_x;
    double _mctrk_end_y;
    double _mctrk_end_z;

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
