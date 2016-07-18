/**
 * \file NuMuCCFilter.h
 *
 * \ingroup CCInclusive
 * 
 * \brief Class def header for a class NuMuCCFilter
 *
 * @author davidkaleko
 */

/** \addtogroup CCInclusive

    @{*/

#ifndef LARLITE_NUMUCCFILTER_H
#define LARLITE_NUMUCCFILTER_H

#include "Analysis/ana_base.h"
#include "GeoAlgo/GeoAABox.h"
#include "FidVolBox.h"

namespace larlite {
  /**
     \class NuMuCCFilter
     User custom analysis class made by SHELL_USER_NAME
   */
  class NuMuCCFilter : public ana_base{
  
  public:

    /// Default constructor
    NuMuCCFilter(){ 
      _name="NuMuCCFilter"; 
      _fout=0;
    _keep_only_numu_from_kaons = false;
    _min_nu_energy = 0.;
    _flip = false;
  }

    /// Default destructor
    virtual ~NuMuCCFilter(){}

    virtual bool initialize();

    virtual bool analyze(storage_manager* storage);

    virtual bool finalize();

    void SetNuMuFromKaonOnly(bool flag){ _keep_only_numu_from_kaons = flag; }
    void SetMinNuEnergy(double minE) { _min_nu_energy = minE; }
    void FlipFilter(bool flag) { _flip = flag; }
    
  protected:
    
    // Fiducial volume box
    geoalgo::AABox _fidvolBox;
 
    size_t total_events;
    size_t kept_events;

    bool _keep_only_numu_from_kaons;
    double _min_nu_energy;

    bool _flip;
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
