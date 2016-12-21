/**
 * \file StoppingPowerSpline.h
 *
 * \ingroup CCInclusive
 *
 * \brief This is kaleko using GraphClick to import the PDG stopping power plot as an image into something we can use
 *  the image is the bottom figure in Figure 1 of ]http://pdg.lbl.gov/2016/AtomicNuclearProperties/adndt.pdf
 *  *** THIS ONLY WORKS FOR MUONS !! ****
 *
 * @author davidkaleko
 */

/** \addtogroup CCInclusive

    @{*/
#ifndef LARLITE_STOPPINGPOWERSPLINE_H
#define LARLITE_STOPPINGPOWERSPLINE_H

#include <iostream>
#include "TGraph.h"
#include "TSpline.h"

/**
   \class StoppingPowerSpline
   User defined class StoppingPowerSpline ... these comments are used to generate
   doxygen documentation!
 */
namespace larlite {

	class StoppingPowerSpline {

	public:

		/// Default constructor
		StoppingPowerSpline();

		/// Default destructor
		~StoppingPowerSpline() {}

		/// Given a muon momentum in GeV, this will tell you the dE/dx that muon is currently depositing in MeV/cm.
		/// For example, a muon with momentum around 0.3 GeV is minimally ionizing so the dEdx should be about 2.105 (assuming LAr)
		/// However a muon with momentum around 10 MeV should have a much higher dEdx because it is almost ranging out
		double GetMudEdx(const double muon_momentum);

		TGraph* GetGraph() { return stopping_power_vs_momentum; }

	private:

		TGraph *stopping_power_vs_momentum;
		TSpline3 *stopping_power_vs_momentum_spline3;
		
	};

}
#endif
/** @} */ // end of doxygen group

