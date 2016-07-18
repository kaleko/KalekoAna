#ifndef FIDVOLBOX_CXX
#define FIDVOLBOX_CXX

#include "FidVolBox.h"
namespace larlite{

	FidVolBox::FidVolBox(){

		double fidvol_dist_x = 10.;
		double fidvol_dist_y = 25.;
		double fidvol_dist_zmin = 10.;
		double fidvol_dist_z = 27.;

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
