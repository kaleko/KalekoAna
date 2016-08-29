#ifndef LARLITE_XIAOEVENTANA_CXX
#define LARLITE_XIAOEVENTANA_CXX

#include "XiaoEventAna.h"
#include "DataFormat/opflash.h"
#include "DataFormat/mctruth.h"
#include "DataFormat/mcflux.h"
#include "DataFormat/mctrack.h"
#include "DataFormat/mcshower.h"

namespace larlite {

    bool XiaoEventAna::initialize() {

        _nu_finder = XiaoNuFinder();
        _nu_finder.setMinTrkLen(_min_trk_len);
        _myspline = TrackMomentumSplines();
        _MCScalc = TrackMomentumCalculator();
        _MCScalc.SetMinLength(_mcs_min_trk_len);
        _nu_E_calc = NuEnergyCalc();
        _nu_E_calc.setMCSMinLen(_mcs_min_trk_len);
        _intxn_booster = IntxnBooster();
        _PID_filler = KalekoPIDFiller();
        _chopper = TrackChopper();
        _trkextender = IntxnTrackExtender();


        if (_filetype == kINPUT_FILE_TYPE_MAX) {
            print(larlite::msg::kERROR, __FUNCTION__, Form("DID NOT SET INPUT FILE TYPE!"));
            return false;
        }
        _nu_finder.setInputType(_filetype);
        _nu_finder.setVtxSphereRadius(_vtx_sphere_radius);
        setBGWTimes();
        // myspline = new TrackMomentumSplines();

        total_events = 0;
        passed_events = 0;

        _fidvolBox = FidVolBox();

        //Box here is TPC
        _tpcBox.Min( 1,
                     -(::larutil::Geometry::GetME()->DetHalfHeight()) + 1,
                     1);

        _tpcBox.Max( 2 * (::larutil::Geometry::GetME()->DetHalfWidth()) - 1,
                     ::larutil::Geometry::GetME()->DetHalfHeight() - 1,
                     ::larutil::Geometry::GetME()->DetLength() - 1);

        // _hmult = new TH1F("hmult", "Track Multiplicity", 10, -0.5, 9.5);
        // _hdedx = new TH2D("hdedx", "End dEdx vs Start dEdx;End dEdx;Start dEdx", 50, 0, 20, 50, 0, 20);
        _hcorrect_ID = new TH1F("hcorrectID", "Was Neutrino Vtx Correctly Identified?", 2, -0.5, 1.5);


        if (!_tree) {
            _tree = new TTree("tree", "tree");
            _tree->Branch("true_nu_pdg", &_true_nu_pdg, "true_nu_pdg/I");
            _tree->Branch("true_nu_E", &_true_nu_E, "true_nu_E/D");
            _tree->Branch("true_nu_CCNC", &_true_nu_CCNC, "true_nu_CCNC/O");
            _tree->Branch("true_nu_mode", &_true_nu_mode, "true_nu_mode/I");
            _tree->Branch("true_multiplicity", &_true_multiplicity, "true_multiplicity/I");
            _tree->Branch("longest_trk_contained", &_longest_trk_contained, "longest_trk_contained/O");
            _tree->Branch("all_trks_contained", &_all_trks_contained, "all_trks_contained/O");
            _tree->Branch("p_phi", &_p_phi, "p_phi/D");
            _tree->Branch("mu_phi", &_mu_phi, "mu_phi/D");
            _tree->Branch("correct_ID", &_correct_ID, "correct_ID/O");
            _tree->Branch("mu_end_dedx", &_mu_end_dedx, "mu_end_dedx/D");
            _tree->Branch("mu_start_dedx", &_mu_start_dedx, "mu_start_dedx/D");
            _tree->Branch("fndecay", &_fndecay, "fndecay/I");
            _tree->Branch("fppdxdz", &_fppdxdz, "fppdxdz/D");
            _tree->Branch("fppdydz", &_fppdydz, "fppdydz/D");
            _tree->Branch("fpppz", &_fpppz, "fpppz/D");
            _tree->Branch("fppenergy", &_fppenergy, "fppenergy/D");
            _tree->Branch("kaon_prod_E", &_kaon_prod_E, "kaon_prod_E/D");
            _tree->Branch("kaon_prod_theta", &_kaon_prod_theta, "kaon_prod_theta/D");
            _tree->Branch("mu_p_dirdot", &_mu_p_dirdot, "mu_p_dirdot/D");
            _tree->Branch("true_lepton_pdg", &_true_lepton_pdg, "true_lepton_pdg/I");
            _tree->Branch("true_lepton_momentum", &_true_lepton_momentum, "true_lepton_momentum/D");
            _tree->Branch("n_associated_tracks", &_n_associated_tracks, "n_associated_tracks/I");
            _tree->Branch("longest_trk_len", &_longest_trk_len, "longest_trk_len/D");
            _tree->Branch("longest_trk_len_infidvol", &_longest_trk_len_infidvol, "longest_trk_len_infidvol/D");
            _tree->Branch("longest_trk_Length", &_longest_trk_Length, "longest_trk_Length/D");
            _tree->Branch("longest_trk_Length_infidvol", &_longest_trk_Length_infidvol, "longest_trk_Length_infidvol/D");
            _tree->Branch("longest_track_end_x_infidvol", &_longest_track_end_x_infidvol, "longest_track_end_x_infidvol/D");
            _tree->Branch("longest_track_end_y_infidvol", &_longest_track_end_y_infidvol, "longest_track_end_y_infidvol/D");
            _tree->Branch("longest_track_end_z_infidvol", &_longest_track_end_z_infidvol, "longest_track_end_z_infidvol/D");
            _tree->Branch("longest_trk_cosy", &_longest_trk_cosy, "longest_trk_cosy/D");
            _tree->Branch("second_longest_trk_len", &_second_longest_trk_len, "second_longest_trk_len/D");
            _tree->Branch("longest_trk_theta", &_longest_trk_theta, "longest_trk_theta/D");
            _tree->Branch("longest_trk_phi", &_longest_trk_phi, "longest_trk_phi/D");
            _tree->Branch("longest_trk_MCS_mom", &_longest_trk_MCS_mom, "longest_trk_MCS_mom/D");
            _tree->Branch("longest_trk_spline_mom", &_longest_trk_spline_mom, "longest_trk_spline_mom/D");
            _tree->Branch("nu_E_estimate", &_nu_E_estimate, "nu_E_estimate/D");
            _tree->Branch("CCQE_E", &_CCQE_E, "CCQE_E/D");
            _tree->Branch("longest_trk_avg_calo", &_longest_trk_avg_calo, "longest_trk_avg_calo/D");
            _tree->Branch("second_longest_trk_avg_calo", &_second_longest_trk_avg_calo, "second_longest_trk_avg_calo/D");
            _tree->Branch("true_nu_x", &_true_nu_x, "true_nu_x/D");
            _tree->Branch("true_nu_y", &_true_nu_y, "true_nu_y/D");
            _tree->Branch("true_nu_z", &_true_nu_z, "true_nu_z/D");
            _tree->Branch("reco_nu_x", &_reco_nu_x, "reco_nu_x/D");
            _tree->Branch("reco_nu_y", &_reco_nu_y, "reco_nu_y/D");
            _tree->Branch("reco_nu_z", &_reco_nu_z, "reco_nu_z/D");
            _tree->Branch("dist_reco_true_vtx", &_dist_reco_true_vtx, "dist_reco_true_vtx/D");
            _tree->Branch("max_tracks_dotprod", &_max_tracks_dotprod, "max_tracks_dotprod/D");
            _tree->Branch("longest_tracks_dotprod", &_longest_tracks_dotprod, "longest_tracks_dotprod/D");
            _tree->Branch("longest_tracks_dotprod_trkendpoints", &_longest_tracks_dotprod_trkendpoints, "longest_tracks_dotprod_trkendpoints/D");
            _tree->Branch("longest_track_end_x", &_longest_track_end_x, "longest_track_end_x/D");
            _tree->Branch("longest_track_end_y", &_longest_track_end_y, "longest_track_end_y/D");
            _tree->Branch("longest_track_end_z", &_longest_track_end_z, "longest_track_end_z/D");
            _tree->Branch("brightest_BSW_flash_PE", &_brightest_BSW_flash_PE, "brightest_BSW_flash_PE/D");
            _tree->Branch("BSW_flash_z_range", &_BSW_flash_z_range, "BSW_flash_z_range/D");
            _tree->Branch("longest_trk_dot_truemuondir", &_longest_trk_dot_truemuondir, "longest_trk_dot_truemuondir/D");
            _tree->Branch("n_reco_nu_in_evt", &_n_reco_nu_in_evt, "n_reco_nu_in_evt/I");
            _tree->Branch("E_lepton", &_E_lepton, "E_lepton/D");
            _tree->Branch("E_hadrons", &_E_hadrons, "E_hadrons/D");
            _tree->Branch("E_MCS", &_E_MCS, "E_MCS/D");
            _tree->Branch("E_range", &_E_range, "E_range/D");
            _tree->Branch("smallest_avg_calo", &_smallest_avg_calo, "smallest_avg_calo/D");
            _tree->Branch("longest_track_well_recod", &_longest_track_well_recod, "longest_track_well_recod/O");
            _tree->Branch("found_mu_mctrack", &_found_mu_mctrack, "found_mu_mctrack/O");
            _tree->Branch("tot_E_mcshowers", &_tot_E_mcshowers, "tot_E_mcshowers/D");
            _tree->Branch("tot_E_mctracks", &_tot_E_mctracks, "tot_E_mctracks/D");
        }


        return true;
    }


    bool XiaoEventAna::setBGWTimes() {

        if (_filetype == kINPUT_FILE_TYPE_MAX) {
            print(larlite::msg::kERROR, __FUNCTION__, Form("DID NOT SET INPUT FILE TYPE!"));
            return false;
        }
        else if ( _filetype == kOnBeam ) {
            BGW_mintime = 3.3;
            BGW_maxtime = 4.9;
        }
        else if ( _filetype == kOffBeam ) {
            BGW_mintime = 3.65;
            BGW_maxtime = 5.25;
        }
        else if ( _filetype == kCorsikaInTime ) {
            BGW_mintime = 3.2;
            BGW_maxtime = 4.8;
        }
        else if ( _filetype == kBNBOnly ) {
            BGW_mintime = 3.55;
            BGW_maxtime = 5.15;
        }
        else if ( _filetype == kBNBCosmic ) {
            BGW_mintime = 3.55;
            BGW_maxtime = 5.15;
        }
        return true;
    }

    void XiaoEventAna::resetTTreeVars() {
        _mu_start_dedx = -999.;
        _mu_end_dedx = -999.;
        _correct_ID = false;
        _mu_phi = -999.;
        _p_phi = -999.;
        _longest_trk_contained = false;
        _all_trks_contained = false;
        _true_nu_E = -999.;
        _true_nu_pdg = -999;
        _true_nu_CCNC = false;
        _true_nu_mode = -999;
        _fndecay = 0;
        _kaon_prod_E = -999.;
        _kaon_prod_theta = -999.;
        _mu_p_dirdot = 999.;
        _true_lepton_pdg = -999;
        _true_lepton_momentum = -999.;
        _n_associated_tracks = 0;
        _longest_trk_len = -999.;
        _longest_trk_len_infidvol = -999.;
        _longest_trk_Length = -999.;
        _longest_trk_Length_infidvol = -999.;
        _second_longest_trk_len = -999.;
        _longest_trk_cosy = -999;
        _longest_trk_theta = -999.;
        _longest_trk_phi = -999.;
        _longest_trk_MCS_mom = -999.;
        _longest_trk_spline_mom = -999.;
        _longest_trk_avg_calo = -999.;
        _second_longest_trk_avg_calo = -999.;
        _nu_E_estimate = -999.;
        _CCQE_E = -999.;
        _true_nu_x = -999.;
        _true_nu_y = -999.;
        _true_nu_z = -999.;
        _reco_nu_x = -999.;
        _reco_nu_y = -999.;
        _reco_nu_z = -999.;
        _dist_reco_true_vtx = -999.;
        _max_tracks_dotprod = -999.;
        _longest_tracks_dotprod = -999.;
        _longest_tracks_dotprod_trkendpoints = -999.;
        _longest_track_end_x = -999.;
        _longest_track_end_y = -999.;
        _longest_track_end_z = -999.;
        _longest_track_end_x_infidvol = -999.;
        _longest_track_end_y_infidvol = -999.;
        _longest_track_end_z_infidvol = -999.;
        _brightest_BSW_flash_PE = -999.;
        _BSW_flash_z_range = -999.;
        _fppdxdz = -999.;
        _fppdydz = -999.;
        _fpppz = -999.;
        _fppenergy = -999.;
        _longest_trk_dot_truemuondir = -999.;
        _n_reco_nu_in_evt = 0;
        _E_lepton = -999.;
        _E_hadrons = -999.;
        _E_MCS = -999.;
        _E_range = -999.;
        _longest_track_well_recod = false;
        _found_mu_mctrack = false;
        _true_multiplicity = 0;
    }

    bool XiaoEventAna::analyze(storage_manager* storage) {


        // if (!(storage->run_id() == 5411 && storage->subrun_id() == 114 && storage->event_id() == 5718)) return false;
        total_events++;
        _n_reco_nu_in_evt = 0;
        resetTTreeVars();

        auto ev_vtx = storage->get_data<event_vertex>(_vtx_producer);
        if (!ev_vtx) {
            print(larlite::msg::kERROR, __FUNCTION__, Form("Did not find specified data product, vertex!"));
            return false;
        }
        if (!ev_vtx->size()) {
            //print(larlite::msg::kERROR, __FUNCTION__, Form("Zero reconstructed vertices in this event!"));
            return false;
        }

//KalekopandoraNuPMAPlustrackalmanhit
        auto ev_track = storage->get_data<event_track>(_track_producer);
        if (!ev_track) {
            print(larlite::msg::kERROR, __FUNCTION__, Form("Did not find specified data product, track!"));
            return false;
        }
        if (!ev_track->size()) {
            //print(larlite::msg::kERROR, __FUNCTION__, Form("Zero reconstructed tracks in this event!"));
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
        //        auto const& ass_calo_v = storage->find_one_ass(ev_track->id(),
        //                                 ev_calo,
        //                                 Form("%scalo", ev_track->name().c_str()));
        std::string track_calo_prod = _calo_producer.substr(0, _calo_producer.find("calo"));
        larlite::product_id myid(data::kTrack, track_calo_prod);
        //  std::cout<<"looking for id with track producer "<< _track_producer<<" and calo producer "<<_calo_producer<<std::endl;
        auto const& ass_calo_v = storage->find_one_ass(myid,
                                 ev_calo,
                                 _calo_producer);

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
        // Update: code allows for multiple neutrinos in each event, so I return a vector of them
        std::vector<larlite::KalekoNuItxn> reco_neutrinos;
        try {
            reco_neutrinos = _nu_finder.findNeutrino(ev_track, ev_calo, ass_calo_v, ev_vtx, ev_opflash);
        }
        catch (...) {
            return false;
        }

        if (!reco_neutrinos.size()) {
            print(larlite::msg::kERROR, __FUNCTION__, Form("Nu Finder returned vector of length 0. Shouldn't happen!"));
            return false;
        }

        _n_reco_nu_in_evt = reco_neutrinos.size();

        // Loop over the reconstructed neutrinos and fill TTree one per neutrino
        for (auto & reco_neutrino : reco_neutrinos) {


            // Let's "boost" the interaction by potentially adding more tracks:
            // _intxn_booster.BoostIntxn(reco_neutrino, ev_track);

            // Let's "extend" the tracks from the vertex by stitching them to new producers
            if (!_extend_track_producer.empty()) {
                auto ev_extend_track = storage->get_data<event_track>(_extend_track_producer);
                if (!ev_extend_track) {
                    print(larlite::msg::kERROR, __FUNCTION__,
                          Form("Did not find specified data product, track for extended track!"));
                    return false;
                }
                if (!ev_extend_track->size())
                    return false;
                _trkextender.ExtendVertexTracks( reco_neutrino, ev_extend_track);

            }
            // Fill the PIDs for the tracks in this interaction
            if ( !_PID_filler.fillKalekoPIDs(reco_neutrino) ) {
                print(larlite::msg::kERROR, __FUNCTION__, Form("Failed filling PIDs for some reason!"));
                throw std::exception();
            }

            // Store TTree variables that have to do with the flashes in the BSW
            _brightest_BSW_flash_PE = -999.;
            _BSW_flash_z_range = -999.;
            for (auto const& flash : *ev_opflash) {
                if (flash.Time() > BGW_mintime && flash.Time() < BGW_maxtime) {
                    if ( flash.TotalPE() > _brightest_BSW_flash_PE ) _brightest_BSW_flash_PE = flash.TotalPE();
                    if ( flash.ZWidth() > _BSW_flash_z_range ) _BSW_flash_z_range = flash.ZWidth();
                }
            }

            _n_associated_tracks = (int)reco_neutrino.Tracks().size();

            _reco_nu_x = reco_neutrino.Vertex().X();
            _reco_nu_y = reco_neutrino.Vertex().Y();
            _reco_nu_z = reco_neutrino.Vertex().Z();
            auto const &geovtx = ::geoalgo::Vector(_reco_nu_x, _reco_nu_y, _reco_nu_z);
            auto longest_trackdir = ::geoalgo::Vector(0., 0., 0.);
            auto longest_trackdir_endpoints = ::geoalgo::Vector(0., 0., 0.);
            _all_trks_contained = true;
            // Quick loop over associated track lengths to find the longest one:
            larlite::track longest_track;
            _longest_trk_len = -999.;
            _longest_trk_len_infidvol = -999.;
            _longest_trk_Length = -999.;
            _longest_trk_Length_infidvol = -999.;
            _smallest_avg_calo = 999.;
            size_t longest_trk_id = 0;
            for (size_t david = 0; david < reco_neutrino.Tracks().size(); david++) {
                auto const asstd_trk = reco_neutrino.Tracks().at(david);
                bool contained = _fidvolBox.Contain(::geoalgo::Vector(asstd_trk.Vertex())) &&
                                 _fidvolBox.Contain(::geoalgo::Vector(asstd_trk.End()));
                if (!contained) _all_trks_contained = false;

                // Compute smallest calo of all vertex-associated tracks
                // to maybe help cut out poorly reconstructed events.
                // Seems some noise sources in data (tracks going thru cathode)
                // might make tracks with poorly reconstructed calos.
                auto const asstd_calo = reco_neutrino.Calos().at(david);
                double dummy_calo = 0;
                for (size_t j = 0; j < asstd_calo.dEdx().size(); ++j)
                    dummy_calo += asstd_calo.dEdx().at(j);
                dummy_calo /= asstd_calo.dEdx().size();
                if (dummy_calo < _smallest_avg_calo) _smallest_avg_calo = dummy_calo;


                double asstd_trk_len = (asstd_trk.End() - asstd_trk.Vertex()).Mag();
                if ( asstd_trk_len > _longest_trk_len ) {
                    longest_track = asstd_trk;
                    auto const chopped_trk = _chopper.chopTrack(asstd_trk);

                    bool flip_longest_trk = ::geoalgo::Vector(asstd_trk.Vertex()).SqDist(geovtx) <
                                            ::geoalgo::Vector(asstd_trk.End()).SqDist(geovtx) ?
                                            false : true;

                    _longest_track_end_x_infidvol = flip_longest_trk ? chopped_trk.Vertex().X() : chopped_trk.End().X();
                    _longest_track_end_y_infidvol = flip_longest_trk ? chopped_trk.Vertex().Y() : chopped_trk.End().Y();
                    _longest_track_end_z_infidvol = flip_longest_trk ? chopped_trk.Vertex().Z() : chopped_trk.End().Z();

                    _longest_trk_len = (asstd_trk.End() - asstd_trk.Vertex()).Mag();
                    _longest_trk_Length = asstd_trk.Length();
                    if (!contained) {
                        _longest_trk_len_infidvol = (chopped_trk.End() - chopped_trk.Vertex()).Mag();
                        _longest_trk_Length_infidvol = chopped_trk.Length();
                    }
                    else {
                        _longest_trk_len_infidvol = _longest_trk_len;
                        _longest_trk_Length_infidvol = _longest_trk_Length;
                    }
                    longest_trk_id = asstd_trk.ID();
                    _longest_trk_theta = (asstd_trk.End() - asstd_trk.Vertex()).Theta();
                    _longest_trk_phi = (asstd_trk.End() - asstd_trk.Vertex()).Phi();
                    _longest_trk_spline_mom = _myspline.GetMuMomentum(asstd_trk_len);

                    longest_trackdir = ::geoalgo::Vector(asstd_trk.Vertex()).SqDist(geovtx) <
                                       ::geoalgo::Vector(asstd_trk.End()).SqDist(geovtx) ?
                                       ::geoalgo::Vector(asstd_trk.VertexDirection()) :
                                       ::geoalgo::Vector(asstd_trk.EndDirection()) * -1.;
                    longest_trackdir_endpoints = ::geoalgo::Vector(asstd_trk.Vertex()).SqDist(geovtx) <
                                                 ::geoalgo::Vector(asstd_trk.End()).SqDist(geovtx) ?
                                                 ::geoalgo::Vector(asstd_trk.End() - asstd_trk.Vertex()) :
                                                 ::geoalgo::Vector(asstd_trk.Vertex() - asstd_trk.End());
                    _longest_trk_contained = contained;


                    _longest_trk_MCS_mom = _MCScalc.GetMomentumMultiScatterLLHD(asstd_trk, flip_longest_trk);

                    _longest_track_end_x = flip_longest_trk ? asstd_trk.Vertex().X() : asstd_trk.End().X();
                    _longest_track_end_y = flip_longest_trk ? asstd_trk.Vertex().Y() : asstd_trk.End().Y();
                    _longest_track_end_z = flip_longest_trk ? asstd_trk.Vertex().Z() : asstd_trk.End().Z();


                    // // Choose the calo object for this track by the one
                    // // with the most number of hits in the dEdx vector
                    // // (this is how analysis tree does it)
                    // int long_track_idx = -1;
                    // for (size_t i = 0; i < ev_track->size(); ++i) {
                    //     auto const& trk = ev_track->at(i);
                    //     if (trk.ID() == asstd_trk.ID()) {
                    //         long_track_idx = i;
                    //         break;
                    //     }
                    // }
                    // if (long_track_idx == -1) {
                    //     std::cout << "ERROR ERROR ERROR DIDNT FIND LONGEST TRACK INDEX BY ID" << std::endl;
                    //     return false;
                    // }

                    auto const asstd_calo = reco_neutrino.Calos().at(david);
                    _longest_trk_avg_calo = 0;
                    for (size_t j = 0; j < asstd_calo.dEdx().size(); ++j)
                        _longest_trk_avg_calo += asstd_calo.dEdx().at(j);
                    _longest_trk_avg_calo /= asstd_calo.dEdx().size();

                    // size_t tmp_nhits = 0;

                    // for (size_t i = 0; i < 3; ++i)
                    //     if (ev_calo->at(ass_calo_v[long_track_idx][i]).dEdx().size() > tmp_nhits) {
                    //         auto const& thecalo = ev_calo->at(ass_calo_v[long_track_idx][i]);
                    //         tmp_nhits = thecalo.dEdx().size();
                    //         _longest_trk_avg_calo = 0;
                    //         for (size_t j = 0; j < tmp_nhits; ++j)
                    //             _longest_trk_avg_calo += thecalo.dEdx().at(j);
                    //         _longest_trk_avg_calo /= tmp_nhits;
                    //     }
                }
            }

            if (longest_trackdir_endpoints.Length() < 0.00001) {
                std::cout << "length problem. " << longest_trackdir_endpoints.Length() << std::endl;
            }

            auto second_longest_trackdir = ::geoalgo::Vector(0., 0., 0.);
            auto second_longest_trackdir_endpoints = ::geoalgo::Vector(0., 0., 0.);
            // Another quick loop to find the second longest one

            for (size_t david = 0; david < reco_neutrino.Tracks().size(); david++) {
                auto const asstd_trk = reco_neutrino.Tracks().at(david);
                // for (auto const& asstd_trk_pair : reco_neutrino.second) {
                // auto const &asstd_trk = asstd_trk_pair.second;
                // Skip the already-found longest track
                if (asstd_trk.ID() == longest_trk_id) continue;
                double asstd_trk_len = (asstd_trk.End() - asstd_trk.Vertex()).Mag();
                if ( asstd_trk_len > _second_longest_trk_len ) {
                    _second_longest_trk_len = asstd_trk_len;
                    second_longest_trackdir = ::geoalgo::Vector(asstd_trk.Vertex()).SqDist(geovtx) <
                                              ::geoalgo::Vector(asstd_trk.End()).SqDist(geovtx) ?
                                              ::geoalgo::Vector(asstd_trk.VertexDirection()) :
                                              ::geoalgo::Vector(asstd_trk.EndDirection()) * -1.;
                    second_longest_trackdir_endpoints = ::geoalgo::Vector(asstd_trk.Vertex()).SqDist(geovtx) <
                                                        ::geoalgo::Vector(asstd_trk.End()).SqDist(geovtx) ?
                                                        ::geoalgo::Vector(asstd_trk.End() - asstd_trk.Vertex()) :
                                                        ::geoalgo::Vector(asstd_trk.Vertex() - asstd_trk.End());

                    auto const asstd_calo = reco_neutrino.Calos().at(david);
                    _second_longest_trk_avg_calo = 0;
                    for (size_t j = 0; j < asstd_calo.dEdx().size(); ++j)
                        _second_longest_trk_avg_calo += asstd_calo.dEdx().at(j);
                    _second_longest_trk_avg_calo /= asstd_calo.dEdx().size();
                }
            }

            // Find the dot product of directions between the longest two tracks
            longest_trackdir_endpoints.Normalize();
            _longest_trk_cosy = longest_trackdir_endpoints.at(1);
            second_longest_trackdir_endpoints.Normalize();
            _longest_tracks_dotprod = longest_trackdir.Dot(second_longest_trackdir);
            _longest_tracks_dotprod_trkendpoints = longest_trackdir_endpoints.Dot(second_longest_trackdir_endpoints);

            // Loop over all track combinations  and compute the highest absolute value of dot product
            // of directions to try to select broken track MIDs

            for (auto const& asstd_trk1 : reco_neutrino.Tracks()) {
                // auto const &asstd_trk1 = asstd_trk_pair1.second;

                auto track1dir = ::geoalgo::Vector(asstd_trk1.Vertex()).SqDist(geovtx) <
                                 ::geoalgo::Vector(asstd_trk1.End()).SqDist(geovtx) ?
                                 ::geoalgo::Vector(asstd_trk1.VertexDirection()) :
                                 ::geoalgo::Vector(asstd_trk1.EndDirection()) * -1.;

                for (auto const& asstd_trk2 : reco_neutrino.Tracks()) {
                    // auto const &asstd_trk2 = asstd_trk_pair2.second;
                    // Don't compare a track to itself.. doing this by track ID
                    if (asstd_trk1.ID() == asstd_trk2.ID()) continue;

                    auto track2dir = ::geoalgo::Vector(asstd_trk2.Vertex()).SqDist(geovtx) <
                                     ::geoalgo::Vector(asstd_trk2.End()).SqDist(geovtx) ?
                                     ::geoalgo::Vector(asstd_trk2.VertexDirection()) :
                                     ::geoalgo::Vector(asstd_trk2.EndDirection()) * -1.;

                    if (fabs(track1dir.Dot(track2dir)) > _max_tracks_dotprod)
                        _max_tracks_dotprod = fabs(track1dir.Dot(track2dir));
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

                auto const mutrack = reco_neutrino.Tracks().at(0).Length() > reco_neutrino.Tracks().at(1).Length() ?
                                     reco_neutrino.Tracks().at(0)          : reco_neutrino.Tracks().at(1);
                auto const ptrack  = reco_neutrino.Tracks().at(0).Length() > reco_neutrino.Tracks().at(1).Length() ?
                                     reco_neutrino.Tracks().at(1)          : reco_neutrino.Tracks().at(0);

                _mu_phi = ::geoalgo::Vector(mutrack.VertexDirection()).Phi();
                _p_phi  = ::geoalgo::Vector( ptrack.VertexDirection()).Phi();

                // There is now an additional cut requiring the dot product between the two track directions
                // is less than 0.95, to reduce broken tracks being recod as 2track events (cosmic background)
                // Make a unit TVector3 for each of the two tracks, ensuring each are pointing away from the vertex
                //            auto const &geovtx = ::geoalgo::Vector(reco_neutrino.first.X(), reco_neutrino.first.Y(), reco_neutrino.first.Z());
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

            // This now fills E lepton and E hadrons
            _nu_E_estimate = _nu_E_calc.ComputeEnuNTracksFromPID(reco_neutrino, _E_lepton, _E_hadrons, _E_MCS, _E_range);
            _CCQE_E = _nu_E_calc.ComputeECCQE(_E_lepton * 1000., longest_trackdir_endpoints, false);

            larlite::mcnu mcnu;
            larlite::mctrack mu_mctrack;
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

                auto ev_mctrack = storage->get_data<event_mctrack>("mcreco");
                if (!ev_mctrack || !ev_mctrack->size()) {
                    print(larlite::msg::kERROR, __FUNCTION__, Form("Did not find specified data product, mctrack (or size == 0)!"));
                    return false;
                }
                auto ev_mcshower = storage->get_data<event_mcshower>("mcreco");
                if (!ev_mcshower || !ev_mcshower->size()) {
                    print(larlite::msg::kERROR, __FUNCTION__, Form("Did not find specified data product, mcshower (or size == 0)!"));
                    return false;
                }

                _fndecay = ev_mcflux->at(0).fndecay;
                _fppdxdz = ev_mcflux->at(0).fppdxdz;
                _fppdydz = ev_mcflux->at(0).fppdydz;
                _fpppz = ev_mcflux->at(0).fpppz;
                _fppenergy = ev_mcflux->at(0).fppenergy;

                double k_plus_mass = 0.493667; // GEV

                if ( _fndecay <= 10 ) {
                    double _kaon_prod_px = _fpppz * _fppdxdz;
                    double _kaon_prod_py = _fpppz * _fppdydz;
                    double _kaon_prod_pz = _fpppz;
                    _kaon_prod_E = std::sqrt( _kaon_prod_px * _kaon_prod_px +
                                              _kaon_prod_py * _kaon_prod_py +
                                              _kaon_prod_pz * _kaon_prod_pz +
                                              k_plus_mass * k_plus_mass);
                    _kaon_prod_theta = ::geoalgo::Vector(
                                           _kaon_prod_px, _kaon_prod_py, _kaon_prod_pz
                                       ).Theta() * (180. / 3.14159);

                }


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
                TVector3 reco_vtx_tvec = TVector3(reco_neutrino.Vertex().X(), reco_neutrino.Vertex().Y(), reco_neutrino.Vertex().Z());
                _dist_reco_true_vtx   = (mcnu.Nu().Trajectory().front().Position().Vect() - reco_vtx_tvec).Mag();
                ::geoalgo::Sphere thevertexsphere(reco_neutrino.Vertex().X(),
                                                  reco_neutrino.Vertex().Y(),
                                                  reco_neutrino.Vertex().Z(),
                                                  5.0);
                _correct_ID = thevertexsphere.Contain(ev_mctruth->at(0).GetNeutrino().Nu().Trajectory().front().Position());
                _hcorrect_ID->Fill(_correct_ID);

                ::geoalgo::Vector true_mu_dir = ::geoalgo::Vector(mcnu.Lepton().Trajectory().front().Momentum().Vect());
                true_mu_dir.Normalize();
                _longest_trk_dot_truemuondir = longest_trackdir.Dot(true_mu_dir);

                // Find the muon mctrack if there is one
                _found_mu_mctrack = false;
                _true_multiplicity = 0;
                _tot_E_mctracks = 0.;
                for (auto const& mct : *ev_mctrack) {
                    if (mct.Origin() != 1 || !mct.size()) continue; //throw out cosmic mctracks

                    double dist_to_true_vtx = (mct.front().Position().Vect() - mcnu.Nu().Trajectory().front().Position().Vect()).Mag();
                    double mct_len = (mct.back().Position().Vect() - mct.front().Position().Vect()).Mag();
                    //true multiplicity is # of mctracks coming from vertex that are longer than 1cm
                    if ( mct_len > 1. and dist_to_true_vtx < 0.01 ) {
                        if (abs(mct.PdgCode()) == 13)
                            _tot_E_mctracks += (mct.front().E());
                        if (abs(mct.PdgCode()) == 211)
                            _tot_E_mctracks += (mct.front().E());
                        if (abs(mct.PdgCode()) == 2212)
                            _tot_E_mctracks += (mct.front().E() - 938.);

                        _true_multiplicity++;
                    }

                    if (abs(mct.PdgCode()) == 13 && dist_to_true_vtx < 0.01 ) {
                        mu_mctrack = mct;
                        _found_mu_mctrack = true;
                    }

                }
                _tot_E_mctracks /= 1000.;

                // count up energy lost to showers
                _tot_E_mcshowers = 0;
                for (auto const& mcs : *ev_mcshower) {
                    if (mcs.Origin() != 1) continue;
                    _tot_E_mcshowers += mcs.DetProfile().E();
                }
                _tot_E_mcshowers /= 1000.;
                // std::cout << "# reco neutrinos in event is " << reco_neutrinos.size() << std::endl;
                // std::cout << "found mu mctrack is " << _found_mu_mctrack << std::endl;
                // std::cout << "size of it is " << mu_mctrack.size() << std::endl;


                // longest track is well recod if the front of reco track matches front of mctrack AND
                // back of reco track matches back of mctrack OR
                // the same but if hte reco track is flipped direction
                if (_found_mu_mctrack)
                    _longest_track_well_recod =
                        ( (mu_mctrack.front().Position().Vect() - longest_track.Vertex()).Mag() < 3. &&
                          (mu_mctrack.back().Position().Vect() -     longest_track.End()).Mag() < 3. ) ||
                        ( (mu_mctrack.front().Position().Vect() -    longest_track.End()).Mag() < 3. &&
                          (mu_mctrack.back().Position().Vect() -  longest_track.Vertex()).Mag() < 3. ) ?
                        true : false;
            }



            // Only apply second longest trk len cut and collinearity cut for ==2 associated tracks
            // if (_n_associated_tracks != 2) {
            //     _second_longest_trk_len = 999.;
            //     _max_tracks_dotprod = 999.;
            //     _longest_tracks_dotprod = 999.;
            //     _longest_tracks_dotprod_trkendpoints = 999.;
            // }

            _tree->Fill();
            passed_events++;
        } // Done loop over all neutrinos in this event
        // std::cout << storage->run_id() << " " << storage->subrun_id() << " " << storage->event_id() << std::endl;
        // if(_longest_trk_contained && _longest_trk_len_infidvol > 100. && _longest_trk_spline_mom < 500 && _longest_trk_MCS_mom > 0.75){//
        //     std::cout << storage->run_id() << " " << storage->subrun_id() << " " << storage->event_id() << std::endl;
        // }

        // std::cout << storage->run_id() << " " << storage->subrun_id() << " " << storage->event_id() << std::endl;
        // if (_n_reco_nu_in_evt == 1 && _all_trks_contained &&  !_true_nu_CCNC && _correct_ID && _longest_trk_len > 100.) {
        //     if (_true_nu_E / _nu_E_estimate > 2) {
        //         std::cout << "Found weird fully contained event." << std::endl;
        //         std::cout << "  - run " << storage->run_id() << ", subrun " << storage->subrun_id() << ", event " << storage->event_id() << std::endl;
        //         std::cout << "  - the (correctly ID'd) reco vertex is at " << Form("(%0.2f,%0.2f,%0.2f)", reco_neutrinos.front().Vertex().X()
        //                   , reco_neutrinos.front().Vertex().Y()
        //                   , reco_neutrinos.front().Vertex().Z()) << std::endl;
        //         std::cout << "  - true neutrino energy is " << _true_nu_E << std::endl;
        //         std::cout << "  - reco neutrino energy is " << _nu_E_estimate << std::endl;
        //         std::cout << "  - and the ttree index is " << storage->get_index() << std::endl;
        //         std::cout << "  - and n associated tracks is " << _n_associated_tracks << std::endl;
        //     }
        // }

        // // if ( _second_longest_trk_len > 16. && fabs(_longest_tracks_dotprod_trkendpoints) > 0.5) {
        // // // if (storage->run_id() == 5411 && storage->subrun_id() == 114 && storage->event_id() == 5718) {
        // if (_correct_ID) {
        //     if (_true_nu_E / _nu_E_estimate > 2 || _true_nu_E / _nu_E_estimate < 0.5 ) {
        //         if (_true_nu_E / _nu_E_estimate > 2 ) std::cout << "TRUE ENERGY MORE THAN TWICE RECO NU ENERGY" << std::endl;
        //         if (_true_nu_E / _nu_E_estimate < 0.5 ) std::cout << "TRUE ENERGY LESS THAN HALF RECO NU ENERGY" << std::endl;
        //         std::cout << storage->run_id() << " " << storage->subrun_id() << " " << storage->event_id() << " (index " << storage->get_index() << ")" << std::endl;
        //         std::cout << "  INFO:" << std::endl;
        //         std::cout << "  - the reco vertex is at " << Form("(%0.2f,%0.2f,%0.2f)", reco_neutrinos.front().Vertex().X()
        //                   , reco_neutrinos.front().Vertex().Y()
        //                   , reco_neutrinos.front().Vertex().Z()) << std::endl;
        //         std::cout << "  - # neutrinos in this event is " << reco_neutrinos.size() << std::endl;
        //         std::cout << "  - (only showing info for the first one)" << std::endl;
        //         std::cout << "  - true neutrino energy is " << _true_nu_E << std::endl;
        //         std::cout << "  - reco neutrino energy is " << _nu_E_estimate << std::endl;
        //         std::cout << "  - and n associated tracks is " << _n_associated_tracks << std::endl;
        //         std::cout << " longest trk len is " << _longest_trk_len << std::endl;
        //         std::cout << " longest trk len contained in fidvol is " << _longest_trk_len_infidvol << std::endl;
        //         std::cout << " longest trk mcs energy is " << _longest_trk_MCS_mom << std::endl;
        //         std::cout << " longest trk dotprod endpoints is " << _longest_tracks_dotprod_trkendpoints << std::endl;
        //         std::cout << " second longest trk len is " << _second_longest_trk_len << std::endl;
        //     }
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
