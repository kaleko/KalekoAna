#ifndef LARLITE_EFFICIENCYSTUDY_CXX
#define LARLITE_EFFICIENCYSTUDY_CXX

#include "EfficiencyStudy.h"
#include "DataFormat/vertex.h"
#include "DataFormat/mctruth.h"
#include "DataFormat/opflash.h"
#include "DataFormat/mctrajectory.h"


namespace larlite {

    bool EfficiencyStudy::initialize() {

        _myGeoAABox = FidVolBox();

        if (!_vtx_tree) {
            _vtx_tree = new TTree("vtxtree", "vtxtree");
            _vtx_tree->Branch("dist_vtx_truth", &_dist_vtx_truth, "dist_vtx_truth/D");
        }

        if (!_evt_tree) {
            _evt_tree = new TTree("evttree", "evttree");
            _evt_tree->Branch("is_truth_fiducial", &_is_truth_fiducial, "is_truth_fiducial/O");
            _evt_tree->Branch("is_numuCC", &_is_numuCC, "is_numuCC/O");
            _evt_tree->Branch("flash_in_bgw", &_flash_in_bgw, "flash_in_bgw/O");
            _evt_tree->Branch("is_areco_vtx_in_fidvol", &_is_areco_vtx_in_fidvol, "_is_areco_vtx_in_fidvol/O");
            _evt_tree->Branch("is_atrack_fromvtx", &_is_atrack_fromvtx, "_is_atrack_fromvtx/O");
            _evt_tree->Branch("longest_track_nearflash_z", &_longest_track_nearflash_z, "_longest_track_nearflash_z/O");
            _evt_tree->Branch("longest_trk_contained", &_longest_trk_contained, "_longest_trk_contained/O");
            _evt_tree->Branch("longest_trk_range_longenough", &_longest_trk_range_longenough, "_longest_trk_range_longenough/O");
            _evt_tree->Branch("longest_trk_start_near_truth_nu", &_longest_trk_start_near_truth_nu, "_longest_trk_start_near_truth_nu/O");
            _evt_tree->Branch("dist_trk_flash", &_dist_trk_flash, "dist_trk_flash/D");
            _evt_tree->Branch("true_nu_E", &_true_nu_E, "true_nu_E/D");
            _evt_tree->Branch("ccqe_nu_E", &_ccqe_nu_E, "ccqe_nu_E/D");
            _evt_tree->Branch("true_mu_E", &_true_mu_E, "true_mu_E/D");
            _evt_tree->Branch("reco_mu_E", &_reco_mu_E, "reco_mu_E/D");
            _evt_tree->Branch("passes_all_cuts", &_passes_all_cuts, "_passes_all_cuts/O");

        }

        if (!_flash_tree) {
            _flash_tree = new TTree("flashtree", "flashtree");
            _flash_tree->Branch("flash_time", &_flash_time, "flash_time/D");
        }


        myspline = new TrackMomentumSplines();
        mycalc = new NuEnergyCalc();

        return true;
    }

    bool EfficiencyStudy::analyze(storage_manager* storage) {

        _dist_vtx_truth = -999.;
        _is_truth_fiducial = false;
        _is_numuCC = false;
        _flash_in_bgw = false;
        _flash_time = -999.;
        _is_areco_vtx_in_fidvol = false;
        _is_atrack_fromvtx = false;
        _longest_track_nearflash_z = false;
        _longest_trk_contained = false;
        _longest_trk_range_longenough = false;
        _longest_trk_start_near_truth_nu = false;
        _dist_trk_flash = 99999.;
        _true_nu_E = -9999.;
        _ccqe_nu_E = -9999.;
        _passes_all_cuts = false;
        _true_mu_E = -9999.;
        _reco_mu_E = -9999.;

        auto ev_vtx = storage->get_data<event_vertex>("pandoraNu");
        if (!ev_vtx) {
            print(larlite::msg::kERROR, __FUNCTION__, Form("Did not find specified data product, vertex!"));
            return false;
        }
        if (!ev_vtx->size()) {
            print(larlite::msg::kERROR, __FUNCTION__, Form("Zero reconstructed vertices in this event!"));
            return false;
        }

        auto ev_track = storage->get_data<event_track>("pandoraNu");
        if (!ev_track) {
            print(larlite::msg::kERROR, __FUNCTION__, Form("Did not find specified data product, track!"));
            return false;
        }
        if (!ev_track->size()) {
            print(larlite::msg::kERROR, __FUNCTION__, Form("Zero reconstructed tracks in this event!"));
            return false;
        }

        auto ev_mctruth = storage->get_data<event_mctruth>("generator");
        if (!ev_mctruth) {
            print(larlite::msg::kERROR, __FUNCTION__, Form("Did not find specified data product, mctruth!"));
            return false;
        }
        if (ev_mctruth->size() != 1) {
            print(larlite::msg::kERROR, __FUNCTION__, Form("MCTruth size is not equal to 1... it equals %lu!", ev_mctruth->size()));
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

        double muon_mass_GEV = 0.10566;

        const mcnu neutrino = ev_mctruth->at(0).GetNeutrino();
        _is_numuCC = !neutrino.CCNC() && neutrino.Nu().PdgCode() == 14 && neutrino.Lepton().PdgCode() == 13;
        _true_nu_E = neutrino.Nu().Trajectory().front().E();
        _true_mu_E = neutrino.Lepton().Trajectory().front().E() - muon_mass_GEV; //muon mass subtracted

        const ::geoalgo::Vector nustart = ::geoalgo::Vector(neutrino.Nu().Trajectory().at(0).X(),
                                          neutrino.Nu().Trajectory().at(0).Y(),
                                          neutrino.Nu().Trajectory().at(0).Z());

        // 10cm fiducial volume
        _is_truth_fiducial = _myGeoAABox.Contain(nustart);

        // Loop over vertices.
        // For each vertex in fiducial volume, loop over reco tracks
        // If you find a reco track starting w/in 5cm from any vertex, store it
        // then if you find a longer reco track, store that one instead
        track longest_track;
        double longest_trk_len = 0;
        for (auto const& vtx : *ev_vtx) {
            ::geoalgo::Vector vertex = ::geoalgo::Vector(vtx.X(), vtx.Y(), vtx.Z());
            _dist_vtx_truth = vertex.Dist(nustart);
            _vtx_tree->Fill();

            if (!_myGeoAABox.Contain(vertex)) continue;

            _is_areco_vtx_in_fidvol = true;

            for (auto const& trk : *ev_track) {
                if (!trk.NumberTrajectoryPoints()) continue;
                ::geoalgo::Vector trkstart = ::geoalgo::Vector(trk.Vertex());
                ::geoalgo::Vector trkend   = ::geoalgo::Vector(trk.End() );
                double dist_trk_start_vtx = vertex.Dist(trkstart);
                double dist_trk_end_vtx = vertex.Dist(trkend);

                if (dist_trk_start_vtx > 5. && dist_trk_end_vtx > 5.) continue;

                _is_atrack_fromvtx = true;
                double trk_len = trkstart.Dist(trkend);
                if (trk_len > longest_trk_len) {
                    longest_track = trk;
                    longest_trk_len = trk_len;
                }
            }
        } //done looping over vertices and finding the longest track

        if (longest_trk_len > 0) {
            _longest_trk_contained = _myGeoAABox.Contain(::geoalgo::Vector(longest_track.Vertex())) &&
                                     _myGeoAABox.Contain(::geoalgo::Vector(longest_track.End()));
            _reco_mu_E = myspline->GetMuMomentum(longest_trk_len) / 1000.;
            std::vector<double> track_dir;
            track_dir.clear();
            for (size_t x = 0; x < 3; ++x) track_dir.push_back(longest_track.VertexDirection()[x]);
            _ccqe_nu_E = mycalc->ComputeECCQE((_reco_mu_E + muon_mass_GEV) * 1000., track_dir, false);
            _longest_trk_range_longenough = longest_trk_len > 75.;

            _longest_trk_start_near_truth_nu = nustart.Dist(::geoalgo::Vector(longest_track.Vertex())) < 5.;


            opflash theflash;
            for (auto const& flash : *ev_opflash) {
                _flash_time = flash.Time();
                _flash_tree->Fill();

                if (flash.Time() > 3.55 && flash.Time() < 5.15 && flash.TotalPE() > 50.) {
                    _flash_in_bgw = true;
                    // Keep track of the brightest flash in the BGW
                    if (flash.TotalPE() > theflash.TotalPE()) {
                        theflash = flash;
                        _dist_trk_flash = flashDistZ(longest_track, theflash.ZCenter());
                        _longest_track_nearflash_z =  _dist_trk_flash < 80.;//10.;

                    }
                }
            }



        }

        _passes_all_cuts = _is_truth_fiducial &&
                           _is_areco_vtx_in_fidvol &&
                           _is_atrack_fromvtx &&
                           _longest_track_nearflash_z &&
                           _longest_trk_contained &&
                           _longest_trk_range_longenough;

        _evt_tree->Fill();


        return true;
    }

    bool EfficiencyStudy::finalize() {

        if (_fout) {
            _fout->cd();
            _vtx_tree->Write();
            _evt_tree->Write();
            _flash_tree->Write();
        }

        return true;
    }

    double EfficiencyStudy::flashDistZ(const track &longest_track, const double flash_z) {

        // If flash is between start/end of track in Z (including flipping track direction)
        // return true

        /// I'm jsut going to copy christoph's exactly to try to duplicate
        // double FlashTrackDist(double flash, double start, double end) {
        double end = longest_track.End().Z();
        double start = longest_track.Vertex().Z();
        double flash = flash_z;
        if (end >= start) {
            if (flash < end && flash > start) return 0;
            else return TMath::Min(fabs(flash - start), fabs(flash - end));
        }
        else {
            if (flash > end && flash < start) return 0;
            else return TMath::Min(fabs(flash - start), fabs(flash - end));
        }
    }
    // if ( (flash_z > longest_track.Vertex().Z() && flash_z < longest_track.End().Z()) ||
    //         (flash_z > longest_track.End().Z() && flash_z < longest_track.Vertex().Z()) )
    //     return 0.;


    // else {
    //     double mymin = std::min(
    //                        std::abs(longest_track.Vertex().Z() - flash_z),
    //                        std::abs(longest_track.End().Z() - flash_z)
    //                    );
    //     return mymin;
    // }
// }

}
#endif
