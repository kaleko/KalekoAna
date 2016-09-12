/**
 * \file MCSBiasStudyDriver.h
 *
 * \ingroup scratch_ana
 *
 * \brief Class def header for a class MCSBiasStudyDriver
 *
 * @author davidkaleko
 */

/** \addtogroup scratch_ana

    @{*/

#ifndef LARLITE_MCSBIASSTUDYDRIVER_H
#define LARLITE_MCSBIASSTUDYDRIVER_H

#include "Analysis/ana_base.h"
#include "DataFormat/track.h"
#include "MCSBiasStudy.h"


namespace larlite {
  /**
     \class MCSBiasStudyDriver
     User custom analysis class made by SHELL_USER_NAME
   */
  class MCSBiasStudyDriver : public ana_base {

  public:

    /// Default constructor
    MCSBiasStudyDriver() { _name = "MCSBiasStudyDriver"; _fout = 0; _module = 0;}

    /// Default destructor
    virtual ~MCSBiasStudyDriver() {}

    virtual bool initialize();

    virtual bool analyze(storage_manager* storage);

    virtual bool finalize();

  protected:

    MCSBiasStudy *_module;
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
