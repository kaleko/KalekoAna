#ifndef LARLITE_TESTCCQECALC_CXX
#define LARLITE_TESTCCQECALC_CXX

#include "TestCCQECalc.h"
#include "DataFormat/mctruth.h"
#include "DataFormat/mctrack.h"

namespace larlite {

  bool TestCCQECalc::initialize() {

    _h = new TH2D("h", "CCQE Energy vs True Energy", 100, 0, 3, 100, 0, 3);

    mycalc = new NuEnergyCalc();

    _myGeoAABox = FidVolBox();

    return true;
  }

  bool TestCCQECalc::analyze(storage_manager* storage) {

    auto ev_mctruth = storage->get_data<event_mctruth>("generator");
    if (!ev_mctruth) {
      print(larlite::msg::kERROR, __FUNCTION__, Form("Did not find specified data product, mctruth!"));
      return false;
    }
    if (ev_mctruth->size() != 1) {
      print(larlite::msg::kERROR, __FUNCTION__, Form("MCTruth size is not equal to 1... it equals %lu!", ev_mctruth->size()));
      return false;
    }


    auto ev_mctrack = storage->get_data<event_mctrack>("mcreco");
    if (!ev_mctrack) {
      print(larlite::msg::kERROR, __FUNCTION__, Form("Did not find specified data product, mctrack!"));
      return false;
    }


    const mcnu neutrino = ev_mctruth->at(0).GetNeutrino();
    bool _is_numuCC = !neutrino.CCNC() && neutrino.Nu().PdgCode() == 14 && neutrino.Lepton().PdgCode() == 13;
    if (!_is_numuCC) return false;

    size_t n_muon_mctracks = 0;
    mctrack mutrack;
    for (auto const& mct : *ev_mctrack) {
      if (mct.PdgCode() == 13 && mct.size() >= 2) {
        mutrack = mct;
        n_muon_mctracks++;
      }
    }

    if (mutrack.size() < 2) return false;

    if ( !_myGeoAABox.Contain(::geoalgo::Vector(mutrack.front().Position().Vect())) ||
         !_myGeoAABox.Contain(::geoalgo::Vector(mutrack.back().Position().Vect()))) return false;
    if (n_muon_mctracks > 1) return false;

    std::vector<double> mu_dir;
    mu_dir.clear();
    mu_dir.push_back(mutrack.front().Momentum().X());
    mu_dir.push_back(mutrack.front().Momentum().Y());
    mu_dir.push_back(mutrack.front().Momentum().Z());

    double muon_mass_GEV = 0.10566;
    double mu_depE = mutrack.front().E() - mutrack.back().E();
    double eccqe = mycalc->ComputeECCQE(mu_depE + muon_mass_GEV * 1000., mu_dir, false);
    double nu_E_true = neutrino.Nu().Trajectory().front().E();

    _h->Fill(nu_E_true, eccqe);

    return true;
  }

  bool TestCCQECalc::finalize() {

    if (_fout) {
      _h->Write();
    }

    return true;
  }

}
#endif
