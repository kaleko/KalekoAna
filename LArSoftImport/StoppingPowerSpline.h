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

		/// Given an initial muon momentum in GeV and a length, this function will tell you the remaining momentum
		/// after that length. For example a 1 GeV muon is roughly MIP so after 10cm it will have lost
		/// about 22 MeV of energy, so its momentum will be about 0.978 GeV (assuming p = E in this comment, which this
		/// function does not)
		/// This function uses the spline to determine the dE/dx at ten evenly spaced points along the length travelled,
		/// then adds up the total dE/dx lost (this is roughly integrating the stopping power spline)
		double Getp(const double initial_muon_momentum, const double length_travelled);

		double GetE(const double initial_muon_momentum, const double length_travelled){
			return p_to_E(Getp(initial_muon_momentum, length_travelled));
		}
		
	private:

		/// Some utility conversion functions to go between momentum and energy
		/// ALWAYS ASSUMING THE PARTICLE IS A MUON (USES MUON MASS)
		double p_to_E(const double p_GEV) { return std::sqrt(p_GEV * p_GEV + _mu_mass_GEV_squared); }
		double E_to_p(const double E_GEV) { return std::sqrt(E_GEV * E_GEV - _mu_mass_GEV_squared); }

		TGraph *stopping_power_vs_momentum;
		TSpline3 *stopping_power_vs_momentum_spline3;

		double _argon_density;
		double _mu_mass_GEV;
		double _mu_mass_GEV_squared;

	};

}
#endif
/** @} */ // end of doxygen group

