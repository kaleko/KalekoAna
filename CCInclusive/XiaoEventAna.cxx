#ifndef LARLITE_XIAOEVENTANA_CXX
#define LARLITE_XIAOEVENTANA_CXX

#include "XiaoEventAna.h"
#include "DataFormat/opflash.h"
#include "DataFormat/mctruth.h"
#include "DataFormat/mcflux.h"

namespace larlite {

    bool XiaoEventAna::initialize() {

        _nu_finder = XiaoNuFinder();
        _myspline = TrackMomentumSplines();
        _nu_E_calc = NuEnergyCalc();

        if (_filetype == kINPUT_FILE_TYPE_MAX) {
            print(larlite::msg::kERROR, __FUNCTION__, Form("DID NOT SET INPUT FILE TYPE!"));
            return false;
        }
        _nu_finder.setInputType(_filetype);

        // myspline = new TrackMomentumSplines();

        total_events = 0;
        passed_events = 0;

        fidvol_dist_x = 20.;
        fidvol_dist_y = 20.;
        fidvol_dist_z = 10.;

        //Box here is TPC
        _fidvolBox.Min( 0 + fidvol_dist_x,
                        -(::larutil::Geometry::GetME()->DetHalfHeight()) + fidvol_dist_y,
                        0 + fidvol_dist_z);

        _fidvolBox.Max( 2 * (::larutil::Geometry::GetME()->DetHalfWidth()) - fidvol_dist_x,
                        ::larutil::Geometry::GetME()->DetHalfHeight() - fidvol_dist_y,
                        ::larutil::Geometry::GetME()->DetLength() - fidvol_dist_z);

        // _hmult = new TH1F("hmult", "Track Multiplicity", 10, -0.5, 9.5);
        // _hdedx = new TH2D("hdedx", "End dEdx vs Start dEdx;End dEdx;Start dEdx", 50, 0, 20, 50, 0, 20);
        _hcorrect_ID = new TH1F("hcorrectID", "Was Neutrino Vtx Correctly Identified?", 2, -0.5, 1.5);


        if (!_tree) {
            _tree = new TTree("tree", "tree");
            _tree->Branch("true_nu_pdg", &_true_nu_pdg, "true_nu_pdg/I");
            _tree->Branch("true_nu_E", &_true_nu_E, "true_nu_E/D");
            _tree->Branch("true_nu_CCNC", &_true_nu_CCNC, "true_nu_CCNC/O");
            _tree->Branch("true_nu_mode", &_true_nu_mode, "true_nu_mode/I");
            _tree->Branch("mu_contained", &_mu_contained, "mu_contained/O");
            _tree->Branch("p_phi", &_p_phi, "p_phi/D");
            _tree->Branch("mu_phi", &_mu_phi, "mu_phi/D");
            _tree->Branch("correct_ID", &_correct_ID, "correct_ID/O");
            _tree->Branch("mu_end_dedx", &_mu_end_dedx, "mu_end_dedx/D");
            _tree->Branch("mu_start_dedx", &_mu_start_dedx, "mu_start_dedx/D");
            _tree->Branch("fndecay", &_fndecay, "fndecay/I");
            _tree->Branch("mu_p_dirdot", &_mu_p_dirdot, "mu_p_dirdot/D");
            _tree->Branch("true_lepton_pdg", &_true_lepton_pdg, "true_lepton_pdg/I");
            _tree->Branch("true_lepton_momentum", &_true_lepton_momentum, "true_lepton_momentum/D");
            _tree->Branch("n_associated_tracks", &_n_associated_tracks, "n_associated_tracks/I");
            _tree->Branch("longest_trk_len", &_longest_trk_len, "longest_trk_len/D");
            _tree->Branch("longest_trk_theta", &_longest_trk_theta, "longest_trk_theta/D");
            _tree->Branch("longest_trk_MCS_mom", &_longest_trk_MCS_mom, "longest_trk_MCS_mom/D");
            _tree->Branch("nu_E_estimate", &_nu_E_estimate, "nu_E_estimate/D");
            _tree->Branch("true_nu_x", &_true_nu_x, "true_nu_x/D");
            _tree->Branch("true_nu_y", &_true_nu_y, "true_nu_y/D");
            _tree->Branch("true_nu_z", &_true_nu_z, "true_nu_z/D");
        }

        return true;
    }

    void XiaoEventAna::resetTTreeVars() {
        _mu_start_dedx = -999.;
        _mu_end_dedx = -999.;
        _correct_ID = false;
        _mu_phi = -999.;
        _p_phi = -999.;
        _mu_contained = false;
        _true_nu_E = -999.;
        _true_nu_pdg = -999;
        _true_nu_CCNC = false;
        _true_nu_mode = -999;
        _fndecay = 0;
        _mu_p_dirdot = 999.;
        _true_lepton_pdg = -999;
        _true_lepton_momentum = -999.;
        _n_associated_tracks = 0;
        _longest_trk_len = -999.;
        _longest_trk_theta = -999.;
        _longest_trk_MCS_mom = -999.;
        _nu_E_estimate = -999.;
        _true_nu_x = -999.;
                _true_nu_y = -999.;
                        _true_nu_z = -999.;
    }

    bool XiaoEventAna::analyze(storage_manager* storage) {

        total_events++;

        resetTTreeVars();

        auto ev_vtx = storage->get_data<event_vertex>("pmtrack");
        if (!ev_vtx) {
            print(larlite::msg::kERROR, __FUNCTION__, Form("Did not find specified data product, vertex!"));
            return false;
        }
        if (!ev_vtx->size()) {
            print(larlite::msg::kERROR, __FUNCTION__, Form("Zero reconstructed vertices in this event!"));
            return false;
        }

        auto ev_track = storage->get_data<event_track>("pandoraNuPMA");
        if (!ev_track) {
            print(larlite::msg::kERROR, __FUNCTION__, Form("Did not find specified data product, track!"));
            return false;
        }
        if (!ev_track->size()) {
            print(larlite::msg::kERROR, __FUNCTION__, Form("Zero reconstructed tracks in this event!"));
            return false;
        }
        auto ev_opflash = storage->get_data<event_opflash>("opflashSat");
        if (!ev_opflash) {
            print(larlite::msg::kERROR, __FUNCTION__, Form("Did not find specified data product, opflash!"));
            return false;
        }
        if (!ev_opflash->size()) {
            print(larlite::msg::kERROR, __FUNCTION__, Form("opflash size is zero!"));
            return false;
        }
        event_calorimetry* ev_calo = nullptr;
        auto const& ass_calo_v = storage->find_one_ass(ev_track->id(),
                                 ev_calo,
                                 Form("%scalo", ev_track->name().c_str()));

        if (!ev_calo) {
            std::cerr << "no event_calorimetry!" << std::endl;
            return false;
        }

        if (ev_calo->empty()) {
            std::cout << "ev_calo empty" << std::endl;
            return false;
        }


        // Try to find a neutrino vertex in this event... return a reco vertex,
        // and a std::vector of tracks that are associated with that vertex
        // std::pair<larlite::vertex, std::vector<larlite::track> > reco_neutrino;
        KalekoNuItxn_t reco_neutrino;
        try {
            reco_neutrino = _nu_finder.findNeutrino(ev_track, ev_calo, ass_calo_v, ev_vtx, ev_opflash);
        }
        catch (...) {
            ///////// TEMP
            // auto ev_mctruth = storage->get_data<event_mctruth>("generator");
            // if (!ev_mctruth) {
            //     print(larlite::msg::kERROR, __FUNCTION__, Form("Did not find specified data product, mctruth!"));
            //     return false;
            // }
            // if (ev_mctruth->size() != 1) {
            //     print(larlite::msg::kERROR, __FUNCTION__, Form("MCTruth size doesn't equal one!"));
            //     return false;
            // }
            // auto const& nuvtx = ev_mctruth->at(0).GetNeutrino().Nu().Trajectory().front().Position();
            // std::cout << " no neutrino reconstructed! tchain index is " << storage->get_index() << std::endl;
            // std::cout<<"    with true nu vtx ("<<nuvtx.X()<<","<<nuvtx.Y()<<","<<nuvtx.Z()<<")"<<std::endl;
            // std::cout<<"    and true nu energy of "<<ev_mctruth->at(0).GetNeutrino().Nu().Trajectory().front().E()<<std::endl;
            // /////////// END TEMP
            return false;
        }


        _n_associated_tracks = (int)reco_neutrino.second.size();
        larlite::mcnu mcnu;
        // If we found a vertex and we are running over MC, let's check if it is accurate
        if (!_running_on_data) {
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
            // std::cout << "The reconstructed vertex is at : " << thevertexsphere.Center() << std::endl;
            // std::cout << "The true vertex is at : "
            //           <<::geoalgo::Vector(ev_mctruth->at(0).GetNeutrino().Nu().Trajectory().front().Position()) << std::endl;
            // auto const& mcnu = ev_mctruth->at(0).GetNeutrino();
            mcnu = ev_mctruth->at(0).GetNeutrino();
            _true_nu_E =            mcnu.Nu().Trajectory().front().E();
            _true_nu_pdg =          mcnu.Nu().PdgCode();
            _true_nu_CCNC =         mcnu.CCNC();
            _true_nu_mode =         mcnu.Mode();
            _true_lepton_pdg =      mcnu.Lepton().PdgCode();
            _true_lepton_momentum = mcnu.Lepton().Trajectory().front().Momentum().Vect().Mag();
            _true_nu_x            = mcnu.Nu().Trajectory().front().Position().X();
            _true_nu_y            = mcnu.Nu().Trajectory().front().Position().Y();
            _true_nu_z            = mcnu.Nu().Trajectory().front().Position().Z();

            ::geoalgo::Sphere thevertexsphere(reco_neutrino.first.X(),
                                              reco_neutrino.first.Y(),
                                              reco_neutrino.first.Z(),
                                              3.0);
            _correct_ID = thevertexsphere.Contain(ev_mctruth->at(0).GetNeutrino().Nu().Trajectory().front().Position());
            _hcorrect_ID->Fill(_correct_ID);
        }

        // Quick loop over associated track lengths to find the longest one:
        for (auto const& asstd_trk_pair : reco_neutrino.second) {
            auto const &asstd_trk = asstd_trk_pair.second;
            if ( asstd_trk.Length() > _longest_trk_len ) {
                _longest_trk_len = asstd_trk.Length();
                _longest_trk_theta = asstd_trk.Theta();
                _longest_trk_MCS_mom = _myspline.GetMuMomentum(asstd_trk.Length());
            }
        }


        // Some ttree entries are only for events with ==2 tracks associated with the vertex:
        if (_n_associated_tracks == 2) {
            //     print(larlite::msg::kWARNING,
            //           __FUNCTION__,
            //           Form("XiaoNuFinder found a neutrino event with %i` != 2 tracks associated with it!",
            //                _n_associated_tracks));

            //     _tree->Fill();
            //     return false;
            // }

            // Now that we found a neutrino interaction with ==2 tracks, let's pick which track is the muon and which is the proton
            // then store some ttree variables about each
            auto const &mutrack = reco_neutrino.second.at(0).second.Length() > reco_neutrino.second.at(1).second.Length() ?
                                  reco_neutrino.second.at(0).second          : reco_neutrino.second.at(1).second;
            auto const &ptrack  = reco_neutrino.second.at(0).second.Length() > reco_neutrino.second.at(1).second.Length() ?
                                  reco_neutrino.second.at(1).second          : reco_neutrino.second.at(0).second;

            _mu_contained = _fidvolBox.Contain(::geoalgo::Vector(mutrack.Vertex())) &&
                            _fidvolBox.Contain(::geoalgo::Vector(mutrack.End()));
            _mu_phi = ::geoalgo::Vector(mutrack.VertexDirection()).Phi();
            _p_phi  = ::geoalgo::Vector( ptrack.VertexDirection()).Phi();

            // There is now an additional cut requiring the dot product between the two track directions
            // is less than 0.95, to reduce broken tracks being recod as 2track events (cosmic background)
            // Make a unit TVector3 for each of the two tracks, ensuring each are pointing away from the vertex
            auto const &geovtx = ::geoalgo::Vector(reco_neutrino.first.X(), reco_neutrino.first.Y(), reco_neutrino.first.Z());
            auto mudir = ::geoalgo::Vector(mutrack.Vertex()).SqDist(geovtx) <
                         ::geoalgo::Vector(mutrack.End()).SqDist(geovtx) ?
                         ::geoalgo::Vector(mutrack.VertexDirection()) :
                         ::geoalgo::Vector(mutrack.EndDirection()) * -1.;
            auto pdir = ::geoalgo::Vector(ptrack.Vertex()).SqDist(geovtx) <
                        ::geoalgo::Vector(ptrack.End()).SqDist(geovtx) ?
                        ::geoalgo::Vector(ptrack.VertexDirection()) :
                        ::geoalgo::Vector(ptrack.EndDirection()) * -1.;
            mudir.Normalize();
            pdir.Normalize();
            _mu_p_dirdot = mudir.Dot(pdir);

        }

        _nu_E_estimate = _nu_E_calc.ComputeEnuNTracksFromPID(reco_neutrino);
        _tree->Fill();
        passed_events++;


        // if (_nu_E_estimate > 4 && _fndecay > 11 ) {
        //     std::cout << "-- found event with true Nu E < 2.4 (" << _true_nu_E << ") "
        //               << "and estimated nu E > 2.4 (" << _nu_E_estimate << ")." << std::endl;
        //     std::cout << "  - the true vertex is at " << Form("(%0.2f,%0.2f,%0.2f)", mcnu.Nu().Trajectory().front().Position().X()
        //               , mcnu.Nu().Trajectory().front().Position().Y()
        //               , mcnu.Nu().Trajectory().front().Position().Z()) << std::endl;
        //     std::cout << "  - the reco vertex is at " << Form("(%0.2f,%0.2f,%0.2f)", reco_neutrino.first.X()
        //               , reco_neutrino.first.Y()
        //               , reco_neutrino.first.Z()) << std::endl;
        //     std::cout << "  - and the ttree index is " << storage->get_index() << std::endl;
        // }

        return true;
    }


    bool XiaoEventAna::finalize() {

        std::cout << "XiaoEventAna finalize! "
                  << total_events << " analyzed in total, "
                  << passed_events << " passed the filter (uhhh meaning == 2 tracks, right now.. check ttree entries)." << std::endl;


        if (_fout) {
            _fout->cd();
            // _hmult->Write();
            // _hdedx->Write();
            _hcorrect_ID->Write();
            _tree->Write();
        }

        _nu_finder.printNumbers();

        return true;
    }

}

#endif
