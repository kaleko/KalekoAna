#ifndef LARLITE_NEUTRINOPARENTINVESIGATION_CXX
#define LARLITE_NEUTRINOPARENTINVESIGATION_CXX

#include "NeutrinoParentInvestigation.h"
#include "DataFormat/mcflux.h"
#include "DataFormat/mctruth.h"

namespace larlite {

  bool NeutrinoParentInvestigation::initialize() {

    if (!_nu_decay_code_hist)
      _nu_decay_code_hist = new TH1I("nu_decay_modes", "Neutrino Decay Modes", 14, 0.5, 14.5);

    if (!_tree) {
      _tree = new TTree("tree", "tree");
      _tree->Branch("pdg", &_pdg, "pdg/I");
      _tree->Branch("fndecay", &_fndecay, "fndecay/I");
      _tree->Branch("E",&_E,"E/D");
    }

    return true;
  }

  bool NeutrinoParentInvestigation::analyze(storage_manager* storage) {


    auto ev_mcflux = storage->get_data<event_mcflux>("generator");

    if (!ev_mcflux) {
      print(larlite::msg::kERROR, __FUNCTION__, Form("Did not find specified data product, MCFlux!"));
      return false;
    }

    auto my_flux = ev_mcflux->at(0);

    _nu_decay_code_hist->Fill(my_flux.fndecay);
    _fndecay = my_flux.fndecay;

    auto ev_mctruth = storage->get_data<event_mctruth>("generator");

    if (!ev_mctruth) {
      print(larlite::msg::kERROR, __FUNCTION__, Form("Did not find specified data product, MCTruth!"));
      return false;
    }
    if (ev_mctruth->size() != 1) {
      return false;
    }

    _pdg = ev_mctruth->at(0).GetNeutrino().Nu().PdgCode();
    _E = ev_mctruth->at(0).GetNeutrino().Nu().Trajectory().front().E();

    _tree->Fill();

    return true;
  }

  bool NeutrinoParentInvestigation::finalize() {

    if (_fout) {
      _fout->cd();
      if (_nu_decay_code_hist)
        _nu_decay_code_hist->Write();
      if (_tree)
        _tree->Write();
    }

    return true;
  }

}
#endif
