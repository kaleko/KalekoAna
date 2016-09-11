/**
 * \file NuEnergyCalc.h
 *
 * \ingroup CCInclusive
 *
 * \brief Class def header for a class NuEnergyCalc
 *
 * @author davidkaleko
 */

/** \addtogroup CCInclusive

    @{*/
#ifndef NUENERGYCALC_H
#define NUENERGYCALC_H


#include <vector>
#include <iostream>
#include <math.h> //pow
#include "TMath.h"
#include "TVector3.h"
#include "GeoAlgo/GeoAABox.h"
#include "FidVolBox.h"
#include "CCInclusiveConstants.h"
#include "KalekoNuItxn.h"
#include "TrackMomentumSplines.h"
#include "TrackMomentumCalculator.h"
#include "TrackChopper.h"
#include "MCSBiasCorrector.h"

/**
   \class NuEnergyCalc
   A class holding functions to calculate neutrino energy. Returns in units of GEV.
 */
namespace larlite {

    class NuEnergyCalc {

    public:

        /// Default constructor
        NuEnergyCalc() {
            _fidvolBox = FidVolBox();
            _myspline = TrackMomentumSplines();
            _tmc = new kaleko::TrackMomentumCalculator();
            _chopper = TrackChopper();
            _corrector = new MCSBiasCorrector();
            _mcs_bias_factor = 1.;
        };

        void setMCSMinLen(double len) { _tmc->SetMinLength(len); }

        /// Default destructor
        virtual ~NuEnergyCalc() {};

        /// Method using manually-input energy (IE if you smear energy first)
        /// Energy should be in MeV, direction can be (doesn't have to be) unit-normalized
        /// Energy is TOTAL ENERGY (including MASS)
        /// Return value is GEV even though input is MEV
        double ComputeECCQE(double energy, const std::vector<double> &lepton_dir, bool is_electron = true);

        double ComputeECCQE(double energy, const TVector3 &lepton_dir, bool is_electron = true);

        /// Method using 4-momentum conservation for numu interacting on neutron (with unknown magnitude
        /// of fermi momentum), with exiting proton and muon.
        /// Input are muon and proton momentum magnitude and direction
        /// The direction can be unit vector, the momentum magnitude should match
        /// E_tot^2 = sqrt( p^2 + m^2 )
        /// and should be in units of MEV
        /// Return value is GEV even though input is MEV
        double ComputeEnu1mu1p(const TVector3 &mu_dir, double mu_mom_mag_MEV, const TVector3 &p_dir, double p_mom_mag_MEV);

        /// Method leveraging Enu=Etot_mu + KE_p agrees with the 4-momentum conservation solution above
        /// It solves one equation for Enu, plugs into the other, and gets an equation for
        /// E_muon as a function of muon direction, proton direction, proton momentum magnitude
        /// Note it is the solution of a quadratic so there is a plus-or-minus in there
        /// Input proton momentum magnitude is in MEV
        /// Return value is in GEV
        double ComputeEmu1mu1pQuadratic(const TVector3 &mu_dir, const TVector3 &p_dir, double p_mom_mag_MEV);
        double ComputeEmu1mu1pQuadratic(const TVector3 &mu_dir, const TVector3 &p_dir, double p_mom_mag_MEV, double delta);
        double ComputeEmu1mu1pQuadratic(double E_tot_n_MEV, double E_tot_p_MEV, double thetamu, double thetap, double p_mom_mag_MEV);
        double ComputeEmu1mu1pQuadratic(double E_tot_n_MEV, double E_tot_p_MEV, double thetamu, double thetap, double p_mom_mag_MEV, double delta);

        /// Same as above but doesn't use Emu = pmu... and algebra is horrific so it solves it approximately, numerically
        double ComputeEmu1mu1pQuadraticNumeric(const TVector3 &mu_dir, const TVector3 &p_dir, double p_mom_mag_MEV);

        /// 1) Use pmu = Emu in equations => Solve to give Emu_1
        /// 2) Use pmu = Emu - (Emu_1 - sqrt(Emu_1^2-Mmu^2)) in equations => Solve to give Emu_2
        /// 3) Use pmu = Emu - (Emu_2 - sqrt(Emu_2^2-Mmu^2)) in equations => Solve to give Emu_3
        double ComputeEmu1mu1pQuadraticIterative(const TVector3 &mu_dir, const TVector3 &p_dir, double p_mom_mag_MEV);
        double ComputeEmu1mu1pQuadraticIterative(double E_tot_n_MEV, double E_tot_p_MEV, double thetamu, double thetap, double p_mom_mag_MEV);

        /// Loop over associated tracks in interaction and their (already determined) PID info.
        /// For contained tracks use track range splines (charged pions use muon spline, protons use proton spline)
        /// If mu track exits, use multiple coloumb scattering
        double ComputeEnuNTracksFromPID(const KalekoNuItxn itxn,
                                        double &E_lepton,
                                        double &E_hadrons,
                                        double &E_MCS,
                                        double &E_range,
                                        bool apply_correction,
                                        bool data_true_MC_false);
        // Wrapper
        double ComputeEnuNTracksFromPID(const KalekoNuItxn itxn);

        bool ViableForCorrection(const double range_energy, const double MCS_energy);

        //depricated
        void SetMCSBiasFactor(double factor) { _mcs_bias_factor = factor; }

    private:
        // Fiducial volume box
        geoalgo::AABox _fidvolBox;

        TrackMomentumSplines _myspline;
        kaleko::TrackMomentumCalculator *_tmc;
        TrackChopper _chopper;
        MCSBiasCorrector *_corrector;

        // Every time MCS is run, the final momentum (energy) is multiplied by this factor
        // this can be measured in data and MC, and set differently for each
        // Default value is 1.
        double _mcs_bias_factor;

    };

}// end namespace larlite
#endif
/** @} */ // end of doxygen group

