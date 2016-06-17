/**
 * \file TestTrkHitAssn.h
 *
 * \ingroup CCInclusive
 * 
 * \brief Class def header for a class TestTrkHitAssn
 *
 * @author davidkaleko
 */

/** \addtogroup CCInclusive

    @{*/

#ifndef LARLITE_TESTTRKHITASSN_H
#define LARLITE_TESTTRKHITASSN_H

#include "Analysis/ana_base.h"
#include "TTree.h"

namespace larlite {

  class TestTrkHitAssn : public ana_base{
  
  public:

    /// Default constructor
    TestTrkHitAssn(){ _name="TestTrkHitAssn"; _fout=0; _tree = 0;}

    /// Default destructor
    virtual ~TestTrkHitAssn(){}

    virtual bool initialize();

    virtual bool analyze(storage_manager* storage);

    virtual bool finalize();

  protected:
    
    TTree *_tree;
    double _sum_hit_ADC;
    double _mct_depE;
    size_t _n_asstd_hits;
    size_t _n_tothits;
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
