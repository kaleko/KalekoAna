#ifndef LARLITE_TESTKINEMATICNURECO_CXX
#define LARLITE_TESTKINEMATICNURECO_CXX

#include "TestKinematicNuReco.h"
#include "DataFormat/mctrack.h"
#include "DataFormat/mctruth.h"

namespace larlite {

    bool TestKinematicNuReco::initialize() {

        _fidvolBox = FidVolBox();

        if (!_tree) {
            _tree = new TTree("tree", "tree");
            _tree->Branch("true_nu_E", &_true_nu_E, "true_nu_E/D");
            // _tree->Branch("reco_nu_E_frompp", &_reco_nu_E_frompp, "reco_nu_E_frompp/D");
            _tree->Branch("plane_Z_anglediff", &_plane_Z_anglediff, "plane_Z_anglediff/D");
            _tree->Branch("true_mu_E", &_true_mu_E, "true_mu_E/D");
            _tree->Branch("reco_mu_E_simple", &_reco_mu_E_simple, "reco_mu_E_simple/D");
            // _tree->Branch("reco_nu_E_frommuE", &_reco_nu_E_frommuE, "reco_nu_E_frommuE/D");
            // _tree->Branch("reco_CCQE_E", &_reco_CCQE_E, "reco_CCQE_E/D");
            _tree->Branch("reco_mu_E_quadratic", &_reco_mu_E_quadratic, "reco_mu_E_quadratic/D");
            _tree->Branch("reco_mu_E_range", &_reco_mu_E_range, "reco_mu_E_range/D");
            _tree->Branch("reco_mu_E_MCS", &_reco_mu_E_MCS, "reco_mu_E_MCS/D");
            _tree->Branch("mu_contained", &_mu_contained, "mu_contained/O");
            _tree->Branch("true_mu_len", &_true_mu_len, "true_mu_len/D");
        }

        mycalc = new NuEnergyCalc();
        myspline = new TrackMomentumSplines();
        myMCScalc = new TrackMomentumCalculator();

        return true;
    }

    bool TestKinematicNuReco::analyze(storage_manager* storage) {

        resetTTreeVars();

        double mumass_MEV = 105.6583715;

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
        if (ev_mctrack->size() != 2) {
            print(larlite::msg::kWARNING, __FUNCTION__, Form("MCTrack size is not equal to 2... it equals %lu!", ev_mctrack->size()));
            return false;
        }

        // Store the muon and the proton mctracks for later investigation
        mctrack mu_mctrack;
        mctrack p_mctrack;
        for (auto const& mct : *ev_mctrack) {

            if (!mct.size()) continue;

            if (mct.PdgCode() == 13)
                mu_mctrack = mct;
            else if (mct.PdgCode() == 2212)
                p_mctrack = mct;
            else {
                print(larlite::msg::kWARNING, __FUNCTION__, Form("Found a mctrack that wasn't p or mu!"));
                return false;
            }
        }

        TVector3 muvec = mu_mctrack.front().Momentum().Vect();
        TVector3 pvec = p_mctrack.front().Momentum().Vect();

        if (!muvec.Mag2() || !pvec.Mag2()) {
            print(larlite::msg::kWARNING, __FUNCTION__, Form("Failed to find both p and mu track!"));
            return false;
        }

        _true_nu_E = ev_mctruth->at(0).GetNeutrino().Nu().Trajectory().front().E();

        _true_mu_E = mu_mctrack.front().E();

        double prot_depE = p_mctrack.front().E() - p_mctrack.back().E();


        TVector3 cross = muvec.Cross(pvec);
        TVector3 zdir(0, 0, 1.);

        _plane_Z_anglediff = cross.Angle(zdir) * (180. / 3.14159265);

        double thetamu = muvec.Angle(zdir);
        double thetap  = pvec.Angle(zdir);

        _reco_nu_E_frompp = ( pvec * (sin(thetap) / tan(thetamu) + cos(thetap)) ).Mag() / 1000.;

        double reco_mu_p = pvec.Mag() * sin(thetap) / sin(thetamu);
        _reco_mu_E_simple = pow(pow(reco_mu_p, 2) + pow(mumass_MEV / 1000., 2), 0.5);

        _reco_nu_E_frommuE = (_reco_mu_E_simple + prot_depE) / 1000.;

        _reco_CCQE_E = mycalc->ComputeECCQE((_reco_mu_E_simple + mumass_MEV), muvec, false);

        // multiple scattering muon energy:
        _reco_mu_E_MCS = myMCScalc->GetMomentumMultiScatterLLHD(mu_mctrack) * 1000.;
        double tracklen = (mu_mctrack.front().Position().Vect() - mu_mctrack.back().Position().Vect()).Mag();
        _reco_mu_E_range = pow(pow(myspline->GetMuMomentum(tracklen), 2) + mumass_MEV * mumass_MEV, 0.5);

        // std::cout<<"true mu E is "<<_true_mu_E<<std::endl;

        // double avg_quad_mu_E = 0.;
        // std::cout<<"proton E is "<<p_mctrack.front().E()<<", thetamu is "<<thetamu<<", thetap is "<<thetap
        // <<", p mag is "<<p_mctrack.front().Momentum().Vect().Mag()<<std::endl;
        // for(size_t i = 0; i <= 10; ++i){
        //     double n_KE = 10. + i*4.;
        //     double computed_E = mycalc->ComputeEmu1mu1pQuadraticIterative(938. + n_KE, p_mctrack.front().E(),
        //                        thetamu, thetap,
        //                        p_mctrack.front().Momentum().Vect().Mag());
        //     avg_quad_mu_E += computed_E;
        // std::cout<<" n_KE = "<<n_KE<<", computed E = "<<computed_E<<std::endl;
        // }
        // avg_quad_mu_E /= 11.;
        _reco_mu_E_quadratic = mycalc->ComputeEmu1mu1pQuadraticIterative(938. + 25., p_mctrack.front().E(),
                               thetamu, thetap,
                               p_mctrack.front().Momentum().Vect().Mag());


        _mu_contained = _fidvolBox.Contain(mu_mctrack.front().Position().Vect()) &&
                        _fidvolBox.Contain(mu_mctrack.back().Position().Vect());

        _true_mu_len = (mu_mctrack.front().Position().Vect() - mu_mctrack.back().Position().Vect()).Mag();

        _tree->Fill();

        return true;
    }

    bool TestKinematicNuReco::finalize() {

        if (_fout) {
            _fout->cd();
            _tree->Write();
        }

        return true;
    }

    void TestKinematicNuReco::resetTTreeVars() {
        _true_nu_E = -999.;
        _reco_nu_E_frompp = -999.;
        _plane_Z_anglediff = -999.;
        _true_mu_E = -999.;
        _reco_mu_E_simple = -999.;
        _reco_nu_E_frommuE = -999.;
        _reco_CCQE_E = -999.;
        _reco_mu_E_MCS = -999.;
        _reco_mu_E_range = -999.;
        _reco_mu_E_quadratic = -999.;
        _mu_contained = false;
        _true_mu_len = -999.;
    }

}
#endif
