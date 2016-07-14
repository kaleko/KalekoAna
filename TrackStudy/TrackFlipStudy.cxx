#ifndef LARLITE_TRACKFLIPSTUDY_CXX
#define LARLITE_TRACKFLIPSTUDY_CXX

#include "TrackFlipStudy.h"
#include "DataFormat/mctrack.h"
#include "DataFormat/track.h"

namespace larlite {

    bool TrackFlipStudy::initialize() {

        _tree = new TTree("tree", "tree");
        _tree->Branch("mccosyangle", &mccosyangle, "mccosyangle/F");
        _tree->Branch("recobackcosyangle", &recobackcosyangle, "recobackcosyangle/F");
        _tree->Branch("recocosyangle", &recocosyangle, "recocosyangle/F");
        _tree->Branch("flipped", &flipped, "flipped/O");
        _tree->Branch("nrecotracks", &nrecotracks, "nrecotracks/I");
        _tree->Branch("nmctracks", &nmctracks, "nmctracks/I");
        _tree->Branch("mclength", &mclength, "mclength/F");
        _tree->Branch("mctheta", &mctheta, "mctheta/F");
        _tree->Branch("mcphi", &mcphi, "mcphi/F");
        _tree->Branch("recolength", &recolength, "recolength/F");
        _tree->Branch("startptdiff", &startptdiff, "startptdiff/F");
        _tree->Branch("startanglediff", &startanglediff, "startanglediff/F");

        _tree->Branch("mcstart_x", &mcstart_x, "mcstart_x/F");
        _tree->Branch("mcstart_y", &mcstart_y, "mcstart_y/F");
        _tree->Branch("mcstart_z", &mcstart_z, "mcstart_z/F");
        _tree->Branch("recostart_x", &recostart_x, "recostart_x/F");
        _tree->Branch("recostart_y", &recostart_y, "recostart_y/F");
        _tree->Branch("recostart_z", &recostart_z, "recostart_z/F");

        return true;
    }

    bool TrackFlipStudy::analyze(storage_manager* storage) {

        mccosyangle = -999.;
        recobackcosyangle = -999.;
        recocosyangle = -999.;
        flipped = false;
        nrecotracks = -1;
        nmctracks = -1;
        mclength = -1.;
        recolength = -1.;
        startptdiff = -1.;
        startanglediff = -1.;
        mcstart_x = -1.;
        mcstart_y = -1.;
        mcstart_z = -1.;
        mctheta = -999.;
        mcphi = -999.;
        recostart_x = -1.;
        recostart_y = -1.;
        recostart_z = -1.;

        ///Read in mctracks
        auto ev_mct = storage->get_data<event_mctrack>("mcreco");
        if (!ev_mct) {
            print(larlite::msg::kERROR, __FUNCTION__, Form("Did not find specified data product, mctrack!"));
            return false;
        }

        ///Read in reco tracks
        auto ev_track = storage->get_data<event_track>("pandoraNuPMA");
        if (!ev_track) {
            print(larlite::msg::kERROR, __FUNCTION__, Form("Did not find specified data product, track!"));
            return false;
        }

        nrecotracks = ev_track->size();
        nmctracks = ev_mct->size();

        if (ev_mct->size() != 1) {
            print(larlite::msg::kERROR, __FUNCTION__, Form("# of mctracks in event != 1"));
            _tree->Fill();
            return false;
        }

        if (!ev_track->size()) {
            print(larlite::msg::kERROR, __FUNCTION__, Form("track exists but has zero size!"));
            _tree->Fill();
            return false;
        }

        auto &mct = ev_mct->at(0);
        mccosyangle = mct.Start().Momentum().Vect().Y() / mct.Start().Momentum().Vect().Mag();
        mclength = -1.*(mct.back().Position() - mct.front().Position()).Mag();
        mcstart_x = mct.Start().Position().X();
        mcstart_y = mct.Start().Position().Y();
        mcstart_z = mct.Start().Position().Z();
        mctheta = (mct.back().Position().Vect() - mct.front().Position().Vect()).Theta();
        mcphi =  (mct.back().Position().Vect() - mct.front().Position().Vect()).Phi();


        auto &track = ev_track->at(0);
        recocosyangle = track.DirectionAtPoint(0).Y() / track.DirectionAtPoint(0).Mag();
        recobackcosyangle = ( -1. * track.EndDirection().Y() ) / track.EndDirection().Mag();
        recolength = (track.LocationAtPoint(0) - track.End()).Mag();
        recostart_x = track.LocationAtPoint(0).X();
        recostart_y = track.LocationAtPoint(0).Y();
        recostart_z = track.LocationAtPoint(0).Z();

        flipped = (mct.Start().Position().Vect() - track.LocationAtPoint(0)).Mag() <
                  (mct.End().Position().Vect() - track.LocationAtPoint(0)).Mag() ? false : true;

        startptdiff = (mct.Start().Position().Vect() - track.LocationAtPoint(0)).Mag();
        startanglediff = acos((mct.Start().Momentum().Vect().Unit().Dot(track.DirectionAtPoint(0).Unit())));

        _tree->Fill();

        return true;
    }

    bool TrackFlipStudy::finalize() {

        if (_fout) { _fout->cd(); _tree->Write(); }

        return true;
    }

}
#endif
