/**
 * \file TestMultiScatterPID.h
 *
 * \ingroup LArSoftImport
 * 
 * \brief Class def header for a class TestMultiScatterPID
 *
 * @author davidkaleko
 */

/** \addtogroup LArSoftImport

    @{*/

#ifndef LARLITE_TESTMULTISCATTERPID_H
#define LARLITE_TESTMULTISCATTERPID_H

#include "Analysis/ana_base.h"
#include "TrackMomentumCalculator.h"

namespace larlite {
  /**
     \class TestMultiScatterPID
     User custom analysis class made by davidkaleko
   */
  class TestMultiScatterPID : public ana_base{
  
  public:

    /// Default constructor
    TestMultiScatterPID(){ _name="TestMultiScatterPID"; _fout=0; _ana_tree=0;};

    /// Default destructor
    virtual ~TestMultiScatterPID(){};

    virtual bool initialize();

    virtual bool analyze(storage_manager* storage);

    virtual bool finalize();

    protected:

    kaleko::TrackMomentumCalculator _tmc;

    TTree* _ana_tree;
    double _llhd;
    double _length;
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
