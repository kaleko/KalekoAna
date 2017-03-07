/**
 * \file BryceCheck.h
 *
 * \ingroup scratch_ana
 * 
 * \brief Class def header for a class BryceCheck
 *
 * @author davidkaleko
 */

/** \addtogroup scratch_ana

    @{*/

#ifndef LARLITE_BRYCECHECK_H
#define LARLITE_BRYCECHECK_H

#include "Analysis/ana_base.h"
#include "FidVolBox.h"
#include "DataFormat/mctrack.h"
#include "DataFormat/mctruth.h"


namespace larlite {
  /**
     \class BryceCheck
     User custom analysis class made by SHELL_USER_NAME
   */
  class BryceCheck : public ana_base{
  
  public:

    /// Default constructor
    BryceCheck(){ _name="BryceCheck"; _fout=0;}

    /// Default destructor
    virtual ~BryceCheck(){}

    /** IMPLEMENT in BryceCheck.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();

    /** IMPLEMENT in BryceCheck.cc! 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);

    /** IMPLEMENT in BryceCheck.cc! 
        Finalize method to be called after all events processed.
    */
    virtual bool finalize();

  protected:
    
    FidVolBox _fidvolBox;
    ::geoalgo::AABox _TPCBox;

    size_t   evt_counter;
    size_t CC_counter;
    size_t in_fidvol;
    size_t is_numu;
    size_t one_muon;
    size_t the_muon_exits_TPC;
    size_t the_muon_longenough;
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
