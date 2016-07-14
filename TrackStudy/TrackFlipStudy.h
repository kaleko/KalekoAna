/**
 * \file TrackFlipStudy.h
 *
 * \ingroup TrackStudy
 *
 * \brief Class def header for a class TrackFlipStudy
 *
 * @author davidkaleko
 */

/** \addtogroup TrackStudy

    @{*/

#ifndef LARLITE_TRACKFLIPSTUDY_H
#define LARLITE_TRACKFLIPSTUDY_H

#include "Analysis/ana_base.h"
#include "TTree.h"

namespace larlite {
  /**
     \class TrackFlipStudy
     User custom analysis class made by SHELL_USER_NAME
   */
  class TrackFlipStudy : public ana_base {

  public:

    /// Default constructor
    TrackFlipStudy() { _name = "TrackFlipStudy"; _fout = 0; _tree = 0;}

    /// Default destructor
    virtual ~TrackFlipStudy() {}

    /** IMPLEMENT in TrackFlipStudy.cc!
        Initialization method to be called before the analysis event loop.
    */
    virtual bool initialize();

    /** IMPLEMENT in TrackFlipStudy.cc!
        Analyze a data event-by-event
    */
    virtual bool analyze(storage_manager* storage);

    /** IMPLEMENT in TrackFlipStudy.cc!
        Finalize method to be called after all events processed.
    */
    virtual bool finalize();

  protected:

    TTree* _tree;
    float mccosyangle;
    float mclength;
    float recobackcosyangle;
    float recocosyangle;
    float recolength;
    bool flipped;
    int nrecotracks;
    int nmctracks;
    float startptdiff;
    float startanglediff;
    float mctheta;
    float mcphi;
    float mcstart_x;
    float mcstart_y;
    float mcstart_z;
    float recostart_x;
    float recostart_y;
    float recostart_z;
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
