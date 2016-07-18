/**
 * \file FidVolBox.h
 *
 * \ingroup CCInclusive
 * 
 * \brief General fiducial volume box object used by lots of other things, 
 *        best to define it all in one place so as not to copy/paste a lot
 *
 * @author davidkaleko
 */

/** \addtogroup CCInclusive

    @{*/
#ifndef FIDVOLBOX_H
#define FIDVOLBOX_H

#include <iostream>
#include "GeoAlgo/GeoAABox.h"
#include "LArUtil/Geometry.h"

/**
   \class FidVolBox
   User defined class FidVolBox ... these comments are used to generate
   doxygen documentation!
 */
namespace larlite{

class FidVolBox : public ::geoalgo::AABox {

public:

  /// Default constructor
  FidVolBox();

  /// Default destructor
  ~FidVolBox(){}

};
}
#endif
/** @} */ // end of doxygen group 

