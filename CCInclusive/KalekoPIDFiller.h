/**
 * \file KalekoPIDFiller.h
 *
 * \ingroup CCInclusive
 * 
 * \brief Class def header for a class KalekoPIDFiller
 *
 * @author davidkaleko
 */

/** \addtogroup CCInclusive

    @{*/

#ifndef LARLITE_KALEKOPIDFILLER_H
#define LARLITE_KALEKOPIDFILLER_H

#include "Analysis/ana_base.h"
#include "CCInclusiveConstants.h"
#include "KalekoNuItxn.h"

namespace larlite {

  class KalekoPIDFiller{
  
  public:

    /// Default constructor
    KalekoPIDFiller(){}

    /// Default destructor
    virtual ~KalekoPIDFiller(){}

    // Takes in reco neutrino interaction, loops through tracks
    // and sets a PID value for each one
    // simple version just says: longest track is muon, <20cm is proton, >20cm is pion
    bool fillKalekoPIDs_simple(KalekoNuItxn &kaleko_itxn);

    // This version uses: longest track muon, then calorimetry associations
    bool fillKalekoPIDs(KalekoNuItxn &kaleko_itxn);

  protected:
    
  };
}
#endif

/** @} */ // end of doxygen group 
