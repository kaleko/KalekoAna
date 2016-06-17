/**
 * \file TestPPiCalo.h
 *
 * \ingroup CCInclusive
 *
 * \brief Class def header for a class TestPPiCalo
 *
 * @author davidkaleko
 */

/** \addtogroup CCInclusive

    @{*/

#ifndef LARLITE_TESTPPICALO_H
#define LARLITE_TESTPPICALO_H

#include "Analysis/ana_base.h"
#include "TH1.h"

namespace larlite {

  class TestPPiCalo : public ana_base {

  public:

    /// Default constructor
    TestPPiCalo() {
      _name = "TestPPiCalo";
      _fout = 0;
      _h_avgdedx_proton = 0;
      _h_avgdedx_pion = 0;
      _h_avgdedx_muon = 0;
    }

    /// Default destructor
    virtual ~TestPPiCalo() {}

    virtual bool initialize();

    virtual bool analyze(storage_manager* storage);

    virtual bool finalize();

  protected:

    TH1F *_h_avgdedx_proton;
    TH1F *_h_avgdedx_pion;
    TH1F *_h_avgdedx_muon;

    size_t n_proton_evts_tot;
    size_t n_pion_evts_tot;
    size_t n_muon_evts_tot;
    
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
