/**
 * \file MCSBiasCorrector.h
 *
 * \ingroup CCInclusive
 *
 * \brief
 *
 *
 * @author davidkaleko
 */

/** \addtogroup CCInclusive

    @{*/
#ifndef MCSBIASCORRECTOR_H
#define MCSBIASCORRECTOR_H

#include <iostream>
#include "TGraph.h"
#include "TSpline.h"

/**
   \class MCSBiasCorrector
   User defined class MCSBiasCorrector ... these comments are used to generate
   doxygen documentation!
 */
namespace larlite {

	class MCSBiasCorrector {

	public:

		/// Default constructor
		MCSBiasCorrector();

		/// Default destructor
		~MCSBiasCorrector() {}

		double GetBiasFactor(const double muon_range, bool data_true_MC_false);


	private:

		TGraph *biasgraph_MC;
		TGraph *biasgraph_data;
		TSpline3 *biasspline_MC;
		TSpline3 *biasspline_data;

	};

}
#endif
/** @} */ // end of doxygen group

