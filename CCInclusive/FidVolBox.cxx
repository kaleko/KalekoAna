#ifndef FIDVOLBOX_CXX
#define FIDVOLBOX_CXX

#include "FidVolBox.h"
namespace larlite{

	FidVolBox::FidVolBox(){

		double fidvol_dist_x = 20.;
		double fidvol_dist_y = 26.5;
		double fidvol_dist_zmin = 20.;
		double fidvol_dist_z = 36.8;//361.8;


		// double fidvol_dist_x = 40.;//20.;
		// double fidvol_dist_y = 30.;
		// double fidvol_dist_zmin = 50.;
		// double fidvol_dist_z = 436.8;//36.8;//361.8;

		
		// DetHalfHeight is 116.5
		// DetHalfWidth  is 128.175
		// DetLength     is 1036.8

		//Box here is TPC
		Min( 0 + fidvol_dist_x,
		                -(::larutil::Geometry::GetME()->DetHalfHeight()) + fidvol_dist_y,
		                0 + fidvol_dist_zmin);

		Max( 2 * (::larutil::Geometry::GetME()->DetHalfWidth()) - fidvol_dist_x,
		                ::larutil::Geometry::GetME()->DetHalfHeight() - fidvol_dist_y,
		                ::larutil::Geometry::GetME()->DetLength() - fidvol_dist_z);

	}
}
#endif
