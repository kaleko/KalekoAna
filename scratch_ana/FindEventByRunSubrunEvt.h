/**
 * \file FindEventByRunSubrunEvt.h
 *
 * \ingroup scratch_ana
 *
 * \brief Class def header for a class FindEventByRunSubrunEvt
 *
 * @author davidkaleko
 */

/** \addtogroup scratch_ana

    @{*/

#ifndef LARLITE_FINDEVENTBYRUNSUBRUNEVT_H
#define LARLITE_FINDEVENTBYRUNSUBRUNEVT_H

#include "Analysis/ana_base.h"

namespace larlite {
  /**
     \class FindEventByRunSubrunEvt
     User custom analysis class made by SHELL_USER_NAME
   */
  class FindEventByRunSubrunEvt : public ana_base {

  public:

    /// Default constructor
    FindEventByRunSubrunEvt() {
      _name = "FindEventByRunSubrunEvt";
      _fout = 0;
      runno = 0;
      subrunno = 0;
      evtno = 0;
    }

    /// Default destructor
    virtual ~FindEventByRunSubrunEvt() {}

    virtual bool initialize();

    virtual bool analyze(storage_manager* storage);

    virtual bool finalize();

    void setRunNo(unsigned int blah) { runno = blah; }
    void setSubRunNo(unsigned int blah) { subrunno = blah; }
    void setEvtNo(unsigned int blah) { evtno = blah; }
    
  protected:

    unsigned int runno;
    unsigned int subrunno;
    unsigned int evtno;

    bool _found_event;
    size_t _file_index;

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
