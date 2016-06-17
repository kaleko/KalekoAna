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

namespace larlite {

  class KalekoPIDFiller{
  
  public:

    /// Default constructor
    KalekoPIDFiller(){}

    /// Default destructor
    virtual ~KalekoPIDFiller(){}

    // Takes in reco neutrino interaction, loops through tracks
    // and sets a PID value for each one
    bool fillKalekoPIDs(KalekoNuItxn_t &kaleko_itxn);

  protected:
    
  };
}
#endif

/** @} */ // end of doxygen group 
