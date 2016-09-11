/**
 * \file NeutrinoParentInvestigation.h
 *
 * \ingroup scratch_ana
 * 
 * \brief Class def header for a class NeutrinoParentInvestigation
 *
 * @author davidkaleko
 */

/** \addtogroup scratch_ana

    @{*/

#ifndef LARLITE_NEUTRINOPARENTINVESIGATION_H
#define LARLITE_NEUTRINOPARENTINVESIGATION_H

#include "Analysis/ana_base.h"
#include "TH1I.h"
#include "TTree.h"

namespace larlite {
  /**
     \class NeutrinoParentInvestigation
     User custom analysis class made by SHELL_USER_NAME
   */
  class NeutrinoParentInvestigation : public ana_base{
  
  public:

    /// Default constructor
    NeutrinoParentInvestigation(){ 
      _name="NeutrinoParentInvestigation"; 
      _fout=0;
      _nu_decay_code_hist=NULL;
      _tree = 0;
    }

    /// Default destructor
    virtual ~NeutrinoParentInvestigation(){}

    /** IMPLEMENT in NeutrinoParentInvestigation.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();

    /** IMPLEMENT in NeutrinoParentInvestigation.cc! 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);

    /** IMPLEMENT in NeutrinoParentInvestigation.cc! 
        Finalize method to be called after all events processed.
    */
    virtual bool finalize();

  protected:

    TH1I *_nu_decay_code_hist;
    TTree *_tree;

    int _pdg;
    int _fndecay;
    double _E;

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
