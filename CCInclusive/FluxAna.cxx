#ifndef LARLITE_FLUXANA_CXX
#define LARLITE_FLUXANA_CXX

#include "FluxAna.h"

namespace larlite {

    bool FluxAna::initialize() {

        if (!_tree) {
            _tree = new TTree("fluxtree", "fluxtree");

            _tree->Branch("true_nu_E", &_true_nu_E, "true_nu_E/D");
            _tree->Branch("true_nu_pdg", &_true_nu_pdg, "true_nu_pdg/I");
            _tree->Branch("true_nu_CCNC", &_true_nu_CCNC, "true_nu_CCNC/O");
            _tree->Branch("true_nu_mode", &_true_nu_mode, "true_nu_mode/I");
            _tree->Branch("true_lepton_pdg", &_true_lepton_pdg, "true_lepton_pdg/I");
            _tree->Branch("true_lepton_momentum", &_true_lepton_momentum, "true_lepton_momentum/D");
            _tree->Branch("fndecay", &_fndecay, "fndecay/I");
            _tree->Branch("fppdxdz", &_fppdxdz, "fppdxdz/D");
            _tree->Branch("fppdydz", &_fppdydz, "fppdydz/D");
            _tree->Branch("fpppz", &_fpppz, "fpppz/D");
            _tree->Branch("fppenergy", &_fppenergy, "fppenergy/D");
            _tree->Branch("kaon_prod_px", &_kaon_prod_px, "kaon_prod_px/D");
            _tree->Branch("kaon_prod_py", &_kaon_prod_py, "kaon_prod_py/D");
            _tree->Branch("kaon_prod_pz", &_kaon_prod_pz, "kaon_prod_pz/D");
            _tree->Branch("kaon_prod_E", &_kaon_prod_E, "kaon_prod_E/D");
            _tree->Branch("kaon_prod_theta", &_kaon_prod_theta, "kaon_prod_theta/D");
            _tree->Branch("nu_in_fidvol", &_nu_in_fidvol, "nu_in_fidvol/O");
        }

        _fidvolBox = FidVolBox();

        return true;
    }

    bool FluxAna::analyze(storage_manager* storage) {

        _true_nu_E = -999.;
        _true_nu_pdg = -999;
        _true_nu_CCNC = false;
        _true_nu_mode = -999;
        _true_lepton_pdg = -999;
        _true_lepton_momentum = -999.;
        _fndecay = -999;
        _fppdxdz = -999.;
        _fppdydz = -999.;
        _fpppz = -999.;
        _fppenergy = -999.;
        _kaon_prod_px = -999.;
        _kaon_prod_py = -999.;
        _kaon_prod_pz = -999.;
        _kaon_prod_E = -999.;
        _kaon_prod_theta = -999.;
        _nu_in_fidvol = false;


        double k_plus_mass = 0.493667; // GEV

        auto ev_mctruth = storage->get_data<event_mctruth>("generator");
        if (!ev_mctruth) {
            print(larlite::msg::kERROR, __FUNCTION__, Form("Did not find specified data product, mctruth!"));
            return false;
        }
        if (ev_mctruth->size() != 1) {
            print(larlite::msg::kERROR, __FUNCTION__, Form("MCTruth size doesn't equal one!"));
            return false;
        }
        auto ev_mcflux = storage->get_data<event_mcflux>("generator");
        if (!ev_mcflux) {
            print(larlite::msg::kERROR, __FUNCTION__, Form("Did not find specified data product, mcflux!"));
            return false;
        }
        // Require exactly one neutrino interaction
        if (ev_mcflux->size() != 1) {
            print(larlite::msg::kINFO, __FUNCTION__, Form("ev_mcflux size is not 1!"));
            return false;
        }
        _fndecay = ev_mcflux->at(0).fndecay;
        _fppdxdz = ev_mcflux->at(0).fppdxdz;
        _fppdydz = ev_mcflux->at(0).fppdydz;
        _fpppz = ev_mcflux->at(0).fpppz;
        _fppenergy = ev_mcflux->at(0).fppenergy;

        if ( _fndecay <= 10 ) {
            _kaon_prod_px = _fpppz * _fppdxdz;
            _kaon_prod_py = _fpppz * _fppdydz;
            _kaon_prod_pz = _fpppz;
            _kaon_prod_E = std::sqrt( _kaon_prod_px * _kaon_prod_px +
                                      _kaon_prod_py * _kaon_prod_py +
                                      _kaon_prod_pz * _kaon_prod_pz +
                                      k_plus_mass * k_plus_mass);
            _kaon_prod_theta = ::geoalgo::Vector(
                                   _kaon_prod_px, _kaon_prod_py, _kaon_prod_pz
                               ).Theta() * (180. / 3.14159);

        }

        auto const mcnu = ev_mctruth->at(0).GetNeutrino();
        _true_nu_E =            mcnu.Nu().Trajectory().front().E();
        _true_nu_pdg =          mcnu.Nu().PdgCode();
        _true_nu_CCNC =         mcnu.CCNC();
        _true_nu_mode =         mcnu.Mode();
        _true_lepton_pdg =      mcnu.Lepton().PdgCode();
        _true_lepton_momentum = mcnu.Lepton().Trajectory().front().Momentum().Vect().Mag();
        _nu_in_fidvol = _fidvolBox.Contain(::geoalgo::Vector(mcnu.Nu().Trajectory().front().Position()));



        _tree->Fill();

        return true;
    }

    bool FluxAna::finalize() {

        if (_fout) { _fout->cd(); _tree->Write(); }

        else
            print(larlite::msg::kERROR, __FUNCTION__, "Did not find an output file pointer!!! File not opened?");


        return true;
    }

}
#endif
