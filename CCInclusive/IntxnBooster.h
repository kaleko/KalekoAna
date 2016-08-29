/**
 * \file IntxnBooster.h
 *
 * \ingroup CCInclusive
 *
 * \brief Takes in an already created neutrino interaction, builds upon it (adds more tracks, etc)
 *
 * @author davidkaleko
 */

/** \addtogroup CCInclusive

    @{*/
#ifndef LARLITE_INTXNBOOSTER_H
#define LARLITE_INTXNBOOSTER_H

#include <iostream>
#include "Analysis/ana_base.h"
#include "DataFormat/track.h"
#include "DataFormat/vertex.h"
#include "DataFormat/calorimetry.h"
#include "KalekoNuItxn.h"
#include "GeoAlgo/GeoVector.h"
#include "GeoAlgo/GeoSphere.h"
#include "GeoAlgo/GeoAABox.h"
#include "CCInclusiveConstants.h"
#include "KalekoPIDFiller.h"
#include "FidVolBox.h"

/**
   \class IntxnBooster
   User defined class IntxnBooster ... these comments are used to generate
   doxygen documentation!
 */
namespace larlite {
	class IntxnBooster : public larlite_base {

	public:

		/// Default constructor
		IntxnBooster() {}

		/// Default destructor
		~IntxnBooster() {}

		void BoostIntxn( larlite::KalekoNuItxn & itxn, const larlite::event_track * ev_track );

		void AddKalmanTracks( larlite::KalekoNuItxn & itxn, const larlite::event_track * ev_track );
	};

}// end namespace larlite
#endif
/** @} */ // end of doxygen group

