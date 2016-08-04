/**
 * \file FluxAna.h
 *
 * \ingroup CCInclusive
 *
 * \brief Class def header for a class FluxAna
 *
 * @author davidkaleko
 */

/** \addtogroup CCInclusive

    @{*/

#ifndef LARLITE_FLUXANA_H
#define LARLITE_FLUXANA_H

#include "Analysis/ana_base.h"
#include "TTree.h"
#include "DataFormat/mcflux.h"
#include "DataFormat/mctruth.h"
#include "FidVolBox.h"

namespace larlite {
  /**
     \class FluxAna
     User custom analysis class made by SHELL_USER_NAME
   */
  class FluxAna : public ana_base {

  public:

    /// Default constructor
    FluxAna() { _name = "FluxAna"; _fout = 0; _tree = 0;}

    /// Default destructor
    virtual ~FluxAna() {}

    virtual bool initialize();

    virtual bool analyze(storage_manager* storage);

    virtual bool finalize();

  protected:

    geoalgo::AABox _fidvolBox;

    TTree *_tree;
    double _true_nu_E;
    int    _true_nu_pdg;
    bool   _true_nu_CCNC;
    int    _true_nu_mode;
    int    _true_lepton_pdg;
    double _true_lepton_momentum;
    int    _fndecay;
    double _fppdxdz;
    double _fppdydz;
    double _fpppz;
    double _fppenergy;
    double _kaon_prod_px;
    double _kaon_prod_py;
    double _kaon_prod_pz;
    double _kaon_prod_E;
    double _kaon_prod_theta;
    bool   _nu_in_fidvol;

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
