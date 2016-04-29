#ifndef LARLITE_TESTKINEMATICNURECO_CXX
#define LARLITE_TESTKINEMATICNURECO_CXX

#include "TestKinematicNuReco.h"
#include "DataFormat/mctrack.h"
#include "DataFormat/mctruth.h"

namespace larlite {

    bool TestKinematicNuReco::initialize() {

        if (!_tree) {
            _tree = new TTree("tree", "tree");
            _tree->Branch("true_nu_E", &_true_nu_E, "true_nu_E/D");
            _tree->Branch("reco_nu_E", &_reco_nu_E, "reco_nu_E/D");
            _tree->Branch("plane_Z_anglediff", &_plane_Z_anglediff, "plane_Z_anglediff/D");
        }

        return true;
    }

    bool TestKinematicNuReco::analyze(storage_manager* storage) {

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

        _true_nu_E = ev_mctruth->at(0).GetNeutrino().Nu().Trajectory().front().E();

        TVector3 muvec;
        TVector3 pvec;

        for (auto const& mct : *ev_mctrack) {
            
            if(!mct.size()) continue;

            if (mct.PdgCode() == 13)
                muvec = mct.front().Momentum().Vect();
            else if (mct.PdgCode() == 2212)
                pvec = mct.front().Momentum().Vect();
            else {
                print(larlite::msg::kWARNING, __FUNCTION__, Form("Found a mctrack that wasn't p or mu!"));
                return false;
            }
        }

        if (!muvec.Mag2() || !pvec.Mag2()) {
            print(larlite::msg::kWARNING, __FUNCTION__, Form("Failed to find both p and mu track!"));
            return false;
        }

        TVector3 cross = muvec.Cross(pvec);
        TVector3 zdir(0,0,1.);

        _plane_Z_anglediff = cross.Angle(zdir) * (180./3.14159265);

        double thetamu = muvec.Angle(zdir);
        double thetap  = pvec.Angle(zdir);

        _reco_nu_E = ( pvec*(sin(thetap)/tan(thetamu) + cos(thetap)) ).Mag() / 1000.;

        _tree->Fill();
        return true;
    }

    bool TestKinematicNuReco::finalize() {

        if(_fout){
            _fout->cd();
            _tree->Write();
        }

        return true;
    }

    void TestKinematicNuReco::resetTTreeVars() {
        _true_nu_E = -999.;
        _reco_nu_E = -999.;
        _plane_Z_anglediff = -999.;
    }

}
#endif
