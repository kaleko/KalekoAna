#ifndef NUENERGYCALC_CXX
#define NUENERGYCALC_CXX

#include "NuEnergyCalc.h"
namespace larlite {


  double NuEnergyCalc::ComputeECCQE(double totalenergy, const std::vector<double> &lepton_dir, bool is_electron) {

    if ( lepton_dir.size() != 3 ) {
      std::cerr << "From ComputeECCQE: input direction vector doesn't have size 3! Quitting..." << std::endl;
      return -99.;
    }

    double M_n = 939.565;    // MeV/c2
    double M_p = 938.272;    // MeV/c2
    double leptonmass = 0.511;      // MeV/c2
    if (!is_electron) leptonmass = 105.6583;

    double bindingE = 30.0;  // MeV

    double l_energy = totalenergy;
    double l_mom = pow(pow(l_energy, 2) - pow(leptonmass, 2), 0.5);

    // Only truth info goes into theta calculation
    double l_theta =
      TMath::ACos(lepton_dir.at(2) /
                  pow(
                    (
                      pow(lepton_dir.at(0), 2) +
                      pow(lepton_dir.at(1), 2) +
                      pow(lepton_dir.at(2), 2)
                    ), 0.5
                  )
                 );

    double nu_energy_num = pow(M_p, 2) - pow(M_n - bindingE, 2)
                           - pow(leptonmass, 2) + 2.0 * (M_n - bindingE) * l_energy;
    double nu_energy_den = 2.0 * (M_n - bindingE - l_energy + l_mom * TMath::Cos(l_theta));

    // For a result in GEV, divide by 1000.
    return ( (nu_energy_num / nu_energy_den) ) / 1000.;

  }

  double NuEnergyCalc::ComputeECCQE(double totalenergy, const TVector3 &lepton_dir, bool is_electron) {

    std::vector<double> temp;
    temp.clear();
    temp.push_back(lepton_dir.X());
    temp.push_back(lepton_dir.Y());
    temp.push_back(lepton_dir.Z());

    return ComputeECCQE(totalenergy, temp, is_electron);

  }

  double NuEnergyCalc::ComputeEnu1mu1p(const TVector3 &mu_dir,
                                       double mu_mom_mag_MEV,
                                       const TVector3 &p_dir,
                                       double p_mom_mag_MEV) {

    // This formula actually uses M_n == M_p
    double M_mu_MEV2 = 11163.676;  // 105.6583 MeV/c2 SQUARED
    double M_p_MEV2 = 880354.345; //  938.272 MeV/c2 SQUARED

    double E_tot_mu = pow(pow(mu_mom_mag_MEV, 2) + M_mu_MEV2, 0.5);
    double E_tot_p  = pow(pow(p_mom_mag_MEV, 2)  + M_p_MEV2,  0.5);

    double thetamu  = mu_dir.Theta();
    double thetap   = p_dir.Theta();
    double thetamup = thetamu + thetap;//mu_dir.Angle(p_dir);
    // std::cout << "meanwhile in compute Enu1mu1p, thetamu, thetap, thetamup are "
    //           << thetamu << ","
    //           << thetap << ","
    //           << thetamup << std::endl;


    double numerator = 0.5 * M_mu_MEV2 + E_tot_mu * E_tot_p - p_mom_mag_MEV * mu_mom_mag_MEV * std::cos(thetamup);
    // std::cout << "numerator is " << numerator << std::endl;
    double denominator = E_tot_p - p_mom_mag_MEV * std::cos(thetap) + E_tot_mu - mu_mom_mag_MEV * std::cos(thetamu);
    // std::cout << "denominator is " << denominator << std::endl;
    // std::cout << "E_tot_p is " << E_tot_p
    //           << ", p_mom_mag_MEV is " << p_mom_mag_MEV
    //           << ", cos thetap is " << std::cos(thetap)
    //           << ", E_tot_mu is " << E_tot_mu
    //           << ", mu_mom_mag_MEV is " << mu_mom_mag_MEV << std::endl;

    if (!denominator) {
      std::cerr << "Wow what are the odds. Denominator in ComputeEnu1mu1p is 0. ERROR!" << std::endl;
      return -1.;
    }

    return (numerator / denominator) / 1000.;
  }

  double NuEnergyCalc::ComputeEmu1mu1pQuadratic(const TVector3 &mu_dir,
      const TVector3 &p_dir,
      double p_mom_mag_MEV)
  {
    // This formula actually uses Emu = pmu which might not be necessary
    // but makes the algebra easier
    // Update: the algebra is incomprehensible (like 50 pages of mathematica output)
    // so I'm writing a new function to solve it numerically

    // double M_mu_MEV2 = 11163.676;  // 105.6583 MeV/c2 SQUARED
    double M_p_MEV2 = 880354.345; //  938.272 MeV/c2 SQUARED
    double M_n_MEV  = 939.57; // MeV/c2
    // double E_b_MEV = 29.5;

    double E_tot_p_MEV = pow(pow(p_mom_mag_MEV, 2)  + M_p_MEV2,  0.5);
    double E_tot_n_MEV = M_n_MEV + 25.;

    double thetamu  = mu_dir.Theta();
    double thetap   = p_dir.Theta();
    // double thetamup = thetamu + thetap;//mu_dir.Angle(p_dir);

    // double Ap = E_tot_p - p_mom_mag_MEV * std::cos(thetap) + E_b;
    // double Bp = E_tot_p - p_mom_mag_MEV * std::cos(thetamup) + E_b;
    // double Cp = E_tot_p - ( M_n_MEV + 23.5 ) + E_b;

    // double quad_a = (1 - std::cos(thetamu));
    // double quad_b = Ap + Cp * (1 - std::cos(thetamu)) - Bp;
    // double quad_c = Cp * Ap - 0.5 * M_mu_MEV2 - E_tot_p * E_b;

    // double descriminant = std::sqrt(quad_b * quad_b - 4 * quad_a * quad_c);
    // double res1 = (-quad_b + descriminant) / (2 * quad_a);
    // double res2 = (-quad_b - descriminant) / (2 * quad_a);

    // res1 /= 1000.;
    // res2 /= 1000.;

    // std::cout << "res1 is " << res1 << ", res2 is " << res2 << std::endl;
    // return res2;

    return ComputeEmu1mu1pQuadratic(E_tot_n_MEV, E_tot_p_MEV, thetamu, thetap, p_mom_mag_MEV);
    // double result = (1 / (2 * (2 - 2 * std::cos(thetamu))));
    // result *= (2 * E_tot_n_MEV - 2 * E_tot_p_MEV - 2 * E_tot_n_MEV * std::cos(thetamu) + 2 * E_tot_p_MEV * std::cos(thetamu)
    //            + 2 * p_mom_mag_MEV * std::cos(thetap) - 2 * p_mom_mag_MEV * std::cos(thetamup)
    //            - std::sqrt(-4 * (2 - 2 * std::cos(thetamu)) * (-M_mu_MEV2 - 2 * E_tot_n_MEV * E_tot_p_MEV + 2 * E_tot_p_MEV * E_tot_p_MEV
    //                        + 2 * E_tot_n_MEV * p_mom_mag_MEV * std::cos(thetap) - 2 * E_tot_p_MEV * p_mom_mag_MEV * std::cos(thetap))
    //                        + pow((-2 * E_tot_n_MEV + 2 * E_tot_p_MEV + 2 * E_tot_n_MEV * std::cos(thetamu) - 2 * E_tot_p_MEV * std::cos(thetamu)
    //                               - 2 * p_mom_mag_MEV * std::cos(thetap) + 2 * p_mom_mag_MEV * std::cos(thetamup)), 2)));
    // return result / 1000.;

  }
  double NuEnergyCalc::ComputeEmu1mu1pQuadratic(const TVector3 &mu_dir,
      const TVector3 &p_dir,
      double p_mom_mag_MEV,
      double delta)
  {

    double M_p_MEV2 = 880354.345; //  938.272 MeV/c2 SQUARED
    double M_n_MEV  = 939.57; // MeV/c2

    double E_tot_p_MEV = pow(pow(p_mom_mag_MEV, 2)  + M_p_MEV2,  0.5);
    double E_tot_n_MEV = M_n_MEV + 25.;

    double thetamu  = mu_dir.Theta();
    double thetap   = p_dir.Theta();

    return ComputeEmu1mu1pQuadratic(E_tot_n_MEV, E_tot_p_MEV, thetamu, thetap, p_mom_mag_MEV, delta);

  }

  double NuEnergyCalc::ComputeEmu1mu1pQuadratic(double E_tot_n_MEV,
      double E_tot_p_MEV,
      double thetamu,
      double thetap,
      double p_mom_mag_MEV)
  {
    double M_mu_MEV2 = 11163.676;  // 105.6583 MeV/c2 SQUARED
    double thetamup = thetamu + thetap;
    double result = (1 / (2 * (2 - 2 * std::cos(thetamu))));
    result *= (2 * E_tot_n_MEV - 2 * E_tot_p_MEV - 2 * E_tot_n_MEV * std::cos(thetamu) + 2 * E_tot_p_MEV * std::cos(thetamu)
               + 2 * p_mom_mag_MEV * std::cos(thetap) - 2 * p_mom_mag_MEV * std::cos(thetamup)
               - std::sqrt(-4 * (2 - 2 * std::cos(thetamu)) * (-M_mu_MEV2 - 2 * E_tot_n_MEV * E_tot_p_MEV + 2 * E_tot_p_MEV * E_tot_p_MEV
                           + 2 * E_tot_n_MEV * p_mom_mag_MEV * std::cos(thetap) - 2 * E_tot_p_MEV * p_mom_mag_MEV * std::cos(thetap))
                           + pow((-2 * E_tot_n_MEV + 2 * E_tot_p_MEV + 2 * E_tot_n_MEV * std::cos(thetamu) - 2 * E_tot_p_MEV * std::cos(thetamu)
                                  - 2 * p_mom_mag_MEV * std::cos(thetap) + 2 * p_mom_mag_MEV * std::cos(thetamup)), 2)));
    return result / 1000.;

  }

  double NuEnergyCalc::ComputeEmu1mu1pQuadratic(double E_tot_n_MEV,
      double E_tot_p_MEV,
      double thetamu,
      double thetap,
      double p_mom_mag_MEV,
      double delta)
  {
    double M_mu_MEV2 = 11163.676;  // 105.6583 MeV/c2 SQUARED
    double thetamup = thetamu + thetap;
    double result = (1 / (2 * (2 - 2 * std::cos(thetamu))));

    result *= (2 * E_tot_n_MEV - 2 * E_tot_p_MEV - 2 * delta * std::cos(thetamu) - 2 * E_tot_n_MEV * std::cos(thetamu) + 2 * E_tot_p_MEV * std::cos(thetamu)
               + 2 * p_mom_mag_MEV * std::cos(thetap) - 2 * p_mom_mag_MEV * std::cos(thetamup)
               - std::sqrt(
                 pow((-2 * E_tot_n_MEV + 2 * E_tot_p_MEV + 2 * delta * std::cos(thetamu) + 2 * E_tot_n_MEV * std::cos(thetamu)
                      - 2 * E_tot_p_MEV * std::cos(thetamu) - 2 * p_mom_mag_MEV * std::cos(thetap)
                      + 2 * p_mom_mag_MEV * std::cos(thetamup)), 2)

                 - 4 * (2 - 2 * std::cos(thetamu)) * (-M_mu_MEV2 - 2 * E_tot_n_MEV * E_tot_p_MEV + 2 * E_tot_p_MEV * E_tot_p_MEV
                     - 2 * delta * E_tot_n_MEV * std::cos(thetamu) + 2 * delta * E_tot_p_MEV * std::cos(thetamu)

                     + 2 * E_tot_n_MEV * p_mom_mag_MEV * std::cos(thetap) - 2 * E_tot_p_MEV * p_mom_mag_MEV * std::cos(thetap)
                     - 2 * delta * p_mom_mag_MEV * std::cos(thetamup))));
    return result / 1000.;

  }

  double NuEnergyCalc::ComputeEmu1mu1pQuadraticNumeric(const TVector3 &mu_dir,
      const TVector3 &p_dir,
      double p_mom_mag_MEV)
  {
    // This formula does not use uses Emu = pmu
    // but then the analytic solution is 50 pages of mathematica output
    // so I'm rough-solving it numerically here

    double M_mu_MEV2 = 11163.676;  // 105.6583 MeV/c2 SQUARED
    double M_p_MEV2 = 880354.345; //  938.272 MeV/c2 SQUARED
    double M_n_MEV  = 939.57; // MeV/c2

    double E_tot_p  = pow(pow(p_mom_mag_MEV, 2)  + M_p_MEV2,  0.5);
    double E_tot_n = M_n_MEV + 23.5;


    double thetamu  = mu_dir.Theta();
    double thetap   = p_dir.Theta();
    double thetamup = thetamu + thetap; //mu_dir.Angle(p_dir);

    double valmin = 99999999999.;
    double imin = -1.;
    for (int i = 0; i < 1500; i += 1) {
      double E_tot_mu = (double)i;
      double mu_mom_mag_MEV = pow(E_tot_mu * E_tot_mu - M_mu_MEV2, 0.5);
      double val = (E_tot_mu + E_tot_p - E_tot_n) *
                   ( 2 * E_tot_p - 2 * p_mom_mag_MEV * std::cos(thetap)
                     + 2 * E_tot_mu - 2 * mu_mom_mag_MEV * std::cos(thetamu) )
                   - M_mu_MEV2 - 2 * E_tot_mu * E_tot_p +
                   2 * mu_mom_mag_MEV * p_mom_mag_MEV * std::cos(thetamup);
      val = std::abs(val);
      // std::cout<<i<<","<<val<<std::endl;
      if (val < valmin) {
        valmin = val;
        imin = E_tot_mu;
      }
    }
    //val should equal 0 if the correct E_tot_mu guess was made
    std::cout << "best fit numerically computed energy is " << imin << std::endl;
    // std::cout<<"valmin is "<<valmin<<", mu total energy computed is "<<imin/1000.<<std::endl;
    return imin / 1000.;
  }

  double NuEnergyCalc::ComputeEmu1mu1pQuadraticIterative(const TVector3 &mu_dir,
      const TVector3 &p_dir, double p_mom_mag_MEV) {

    double thetamu = mu_dir.Theta();
    double thetap = p_dir.Theta();
    double M_p_MEV2 = 880354.345; //  938.272 MeV/c2 SQUARED
    double E_tot_p_MEV = std::sqrt(p_mom_mag_MEV * p_mom_mag_MEV + M_p_MEV2);
    double E_tot_n_MEV = 938. + 25.;

    return ComputeEmu1mu1pQuadraticIterative(E_tot_n_MEV, E_tot_p_MEV, thetamu, thetap, p_mom_mag_MEV);
    // double M_mu_MEV2 = 11163.676;
    // double M_mu_MEV  = 105.65837;

    // double muEguess = ComputeEmu1mu1pQuadratic(mu_dir, p_dir, p_mom_mag_MEV) * 1000.;

    // for (size_t i = 0; i < 10; ++i) {
    //   double delta = muEguess - std::max(M_mu_MEV, std::sqrt(muEguess * muEguess - M_mu_MEV2));
    //   muEguess = ComputeEmu1mu1pQuadratic(mu_dir, p_dir, p_mom_mag_MEV, delta) * 1000.;
    // }

    // std::cout << "Self-contained iterative approach thinks that muE guess is " << muEguess << std::endl;
    // return muEguess;
  }

  double NuEnergyCalc::ComputeEmu1mu1pQuadraticIterative(double E_tot_n_MEV,
      double E_tot_p_MEV,
      double thetamu,
      double thetap,
      double p_mom_mag_MEV) {

    double M_mu_MEV2 = 11163.676;
    double M_mu_MEV  = 105.65837;

    double muEguess = ComputeEmu1mu1pQuadratic(E_tot_n_MEV, E_tot_p_MEV, thetamu, thetap, p_mom_mag_MEV) * 1000.;

    for (size_t i = 0; i < 10; ++i) {
      double delta = muEguess - std::max(M_mu_MEV, std::sqrt(muEguess * muEguess - M_mu_MEV2));
      muEguess = ComputeEmu1mu1pQuadratic(E_tot_n_MEV, E_tot_p_MEV, thetamu, thetap, p_mom_mag_MEV, delta) * 1000.;
    }

    // std::cout << "Self-contained iterative approach thinks that muE guess is " << muEguess << std::endl;
    return muEguess;
  }

  double NuEnergyCalc::ComputeEnuNTracksFromPID(const KalekoNuItxn itxn){
    double dummy1 = 0.;
    double dummy2 = 0.;
    return ComputeEnuNTracksFromPID(itxn,dummy1,dummy2);
  }
  
  double NuEnergyCalc::ComputeEnuNTracksFromPID(const KalekoNuItxn itxn, double &E_lepton, double &E_hadrons) {

    bool debug = false;
    double tot_nu_energy = 0.;
    E_lepton = 0.;
    E_hadrons = 0.;
    // Loop over associated tracks.
    // They are not necessarily pointing in the correct direction ...
    // Assume the start/end point closest to the interaction vertex is the start of the track
    if (debug) std::cout << ":::NuEnergyCalc::: looping over " << itxn.Tracks().size() << " tracks ... " << std::endl;
    for (size_t itrk = 0; itrk < itxn.Tracks().size(); ++itrk) {
      auto const track = itxn.Tracks().at(itrk);

      //let's try using start-to-end distance instead
      // double itrklen = track.Length();
      double itrklen = (track.End() - track.Vertex()).Mag();
      bool trkcontained = _fidvolBox.Contain(track.Vertex()) && _fidvolBox.Contain(track.End());


      auto const &geovtx = ::geoalgo::Vector(itxn.Vertex().X(), itxn.Vertex().Y(), itxn.Vertex().Z());
      bool flip_trk = ::geoalgo::Vector(track.Vertex()).SqDist(geovtx) <
                      ::geoalgo::Vector(track.End()).SqDist(geovtx) ?
                      false : true;


      if (itxn.PIDs().at(itrk) == kKalekoMuon) {
        if (debug) std::cout << " :::NuEnergyCalc::: this track is thought to be a muon " << std::endl;
        if (trkcontained) {
          tot_nu_energy += _myspline.GetMuMomentum(itrklen) / 1000. + 0.106;
          E_lepton = _myspline.GetMuMomentum(itrklen) / 1000. + 0.106;
          if (debug) std::cout << " :::NuEnergyCalc::: this track contained. added spline energy of "
                                 << _myspline.GetMuMomentum(itrklen) / 1000. + 0.106 << std::endl;
        }
        else {
          if (debug) std::cout << " :::NuEnergyCalc::: this track NOT contained." << std::endl;
          double mcs_energy = _tmc.GetMomentumMultiScatterLLHD(track,flip_trk) + 0.106;

          //New addition: if range energy is more than MCS energy, then always use range energy!
          // there's no way range energy is going to overestimate.
          double spline_energy = _myspline.GetMuMomentum(itrklen) / 1000. + 0.106;
          if (spline_energy > mcs_energy) {
            tot_nu_energy += spline_energy;
            E_lepton = spline_energy;
          }
          else if (mcs_energy > 0) {
            tot_nu_energy += mcs_energy;
            E_lepton = mcs_energy;
            if (debug) std::cout << " :::NuEnergyCalc::: MCS worked fine (length = " << itrklen
                                   << ". adding in energy of " << mcs_energy + 0.106 << std::endl;
          }
          else {
            std::cout << "SHOULD NEVER GET HERE" << std::endl;
            // if (debug) std::cout << " :::NuEnergyCalc::: MCS failed. using tracklength contained. adding "
            //                        << _myspline.GetMuMomentum(itrklen) / 1000. + 0.106 << std::endl;
          }
        }
      }
      else if (itxn.PIDs().at(itrk) == kKalekoChargedPion) {
        if (debug) std::cout << " :::NuEnergyCalc::: this track is thought to be a pion " << std::endl;
        if (trkcontained) {
          tot_nu_energy += _myspline.GetMuMomentum(itrklen) / 1000. + 0.140;
          E_hadrons += _myspline.GetMuMomentum(itrklen) / 1000. + 0.140;
          if (debug) std::cout << " :::NuEnergyCalc::: this track contained. added spline energy of "
                                 << _myspline.GetMuMomentum(itrklen) / 1000. + 0.140 << std::endl;
        }
        else {
          if (debug) std::cout << " :::NuEnergyCalc::: this track NOT contained." << std::endl;
          double mcs_energy = _tmc.GetMomentumMultiScatterLLHD(track,flip_trk) + 0.140;

          //New addition: if range energy is more than MCS energy, then always use range energy!
          // there's no way range energy is going to overestimate.
          double spline_energy = _myspline.GetMuMomentum(itrklen) / 1000. + 0.140;
          if (spline_energy > mcs_energy) {
            tot_nu_energy += spline_energy;
            E_hadrons += spline_energy;
          }
          else if (mcs_energy > 0) {
            tot_nu_energy += mcs_energy;
            E_hadrons += mcs_energy;
            if (debug) std::cout << " :::NuEnergyCalc::: MCS worked fine. adding in energy of " << mcs_energy + 0.140 << std::endl;
          }
          else {
            std::cout << "SHOULD NEVER GET HERE (pion loop)" << std::endl;
            // if (debug) std::cout << " :::NuEnergyCalc::: MCS failed. using tracklength contained. adding "
            //                        << _myspline.GetMuMomentum(itrklen) / 1000. + 0.146 << std::endl;
          }
        }
      }
      else if (itxn.PIDs().at(itrk) == kKalekoProton) {
        if (debug) std::cout << " :::NuEnergyCalc::: this track is thought to be a proton " << std::endl;
        tot_nu_energy += _myspline.GetPMomentum(itrklen) / 1000.;
        E_hadrons += _myspline.GetPMomentum(itrklen) / 1000.;
        if (debug) std::cout << " :::NuEnergyCalc::: adding in energy from spline: "
                               << _myspline.GetPMomentum(itrklen) / 1000. << std::endl;
      }
      else {
        std::cerr << "UHHH NuEnergyCalc::ComputeEnuNTracksFromPID was handed tracks that didn't have PIDs associated." << std::endl;
        return -1.;
      } // end if no PID
    } // end loop over tracks

    return tot_nu_energy;
  } // end computeenuNtracksfromPID
}

#endif
