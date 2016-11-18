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
#include "DataFormat/mctrack.h"
#include "DataFormat/mctruth.h"
#include "MCSBiasStudy.h"
#include "FidVolBox.h"
#include "TTree.h"

namespace larlite {
  /**
     \class MCSBiasStudyDriver
     User custom analysis class made by SHELL_USER_NAME
   */
  class MCSBiasStudyDriver : public ana_base {

  public:

    /// Default constructor
    MCSBiasStudyDriver() { 
      _name = "MCSBiasStudyDriver";
     _fout = 0;
     _module = 0;
     _ana_type = MCSBiasStudy::kANALYSIS_TYPE_MAX;
     _driver_tree = 0;
     }

    /// Default destructor
    virtual ~MCSBiasStudyDriver() {}

    virtual bool initialize();

    virtual bool analyze(storage_manager* storage);

    virtual bool finalize();

    void SetAnalysisType(MCSBiasStudy::AnalysisType_t mytype) { _ana_type = mytype; }

  protected:

    MCSBiasStudy::AnalysisType_t _ana_type;

    MCSBiasStudy *_module;

    FidVolBox _fidvol;

    TTree *_driver_tree;
    int _MCT_PDG;
    int _MCT_origin; //1 neutrino, 2 cosmic
    int _run;
    int _subrun;
    int _eventid;
    size_t counter;
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
