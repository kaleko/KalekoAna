#ifndef LARLITE_EVENTSCONTAINEDSTUDY_CXX
#define LARLITE_EVENTSCONTAINEDSTUDY_CXX

#include "EventsContainedStudy.h"
#include "DataFormat/mctrack.h"
#include "DataFormat/mctruth.h"
#include "DataFormat/track.h"
#include "DataFormat/vertex.h"
#include "DataFormat/mcflux.h"

namespace larlite {

	bool EventsContainedStudy::initialize() {

		_fidvolBox = FidVolBox();

		_evts_pass_filter = 0;
		_evts_fully_contained = 0;


		_h_mu_E_contained = new TH1F("_h_mu_E_contained", "Muon True Start Energy: Contained Events [GeV]", 100, 0, 3.);
		_h_mu_E_all = new TH1F("_h_mu_E_all", "Muon True Start Energy: All Events [GeV]", 100, 0, 3.);
		_h_nu_E_contained = new TH1F("_h_nu_E_contained", "Neutrino True Energy: Contained Events [GeV]", 100, 0, 3.);
		_h_nu_E_all = new TH1F("_h_nu_E_all", "Neutrino True Energy: All Events [GeV]", 100, 0, 3.);

		_tree = new TTree("tree", "tree");
		_tree->Branch("mu_contained", &_mu_contained, "mu_contained/O");
		_tree->Branch("p_contained", &_p_contained, "p_contained/O");
		_tree->Branch("n_reco_tracks_fromvtx", &_n_reco_tracks_fromvtx, "n_reco_tracks_fromvtx/I");
		_tree->Branch("n_reco_tracks_endingnearvtx", &_n_reco_tracks_endingnearvtx, "n_reco_tracks_endingnearvtx/I");
		_tree->Branch("n_mc_tracks_fromvtx", &_n_mc_tracks_fromvtx, "n_mc_tracks_fromvtx/I");
		_tree->Branch("true_mu_E", &_true_mu_E, "true_mu_E/D");
		_tree->Branch("true_p_E", &_true_p_E, "true_p_E/D");
		_tree->Branch("reco_mu_E", &_reco_mu_E, "reco_mu_E/D");
		_tree->Branch("reco_p_E", &_reco_p_E, "reco_p_E/D");
		_tree->Branch("true_nu_E", &_true_nu_E, "true_nu_E/D");
		_tree->Branch("true_nu_x", &_true_nu_x, "true_nu_x/D");
		_tree->Branch("true_nu_y", &_true_nu_y, "true_nu_y/D");
		_tree->Branch("true_nu_z", &_true_nu_z, "true_nu_z/D");
		_tree->Branch("reco_CCQE_E", &_reco_CCQE_E, "reco_CCQE_E/D");
		_tree->Branch("true_mu_len", &_true_mu_len, "true_mu_len/D");
		_tree->Branch("true_p_len", &_true_p_len, "true_p_len/D");
		_tree->Branch("reco_mu_len", &_reco_mu_len, "reco_mu_len/D");
		_tree->Branch("reco_p_len", &_reco_p_len, "reco_p_len/D");
		_tree->Branch("runno", &runno, "runno/I");
		_tree->Branch("subrunno", &subrunno, "subrunno/I");
		_tree->Branch("evtno", &evtno, "evtno/I");
		_tree->Branch("true_neutrons_E", &_true_neutrons_E, "true_neutrons_E/D");
		_tree->Branch("is_reco_vtx_near_true_vtx", &_is_reco_vtx_near_true_vtx, "is_reco_vtx_near_true_vtx/O");
		_tree->Branch("fndecay", &_fndecay, "fndecay/I");

		myspline = new TrackMomentumSplines();
		mycalc = new NuEnergyCalc();

		kaleko_ctr = 0;
		return true;
	}

	bool EventsContainedStudy::analyze(storage_manager* storage) {

		runno = storage->run_id();
		subrunno = storage->subrun_id();
		evtno = storage->event_id();

		kaleko_ctr++;

		_evts_pass_filter++;

		_mu_contained = false;
		_p_contained = false;
		_n_mc_tracks_fromvtx = 0;
		_n_reco_tracks_fromvtx = 0;
		_n_reco_tracks_endingnearvtx = 0;
		_true_mu_E = -999.;
		_true_p_E = -999.;
		_reco_mu_E = -999.;
		_reco_p_E = -999.;
		_true_nu_E = -999.;
		_true_nu_x = -999.;
		_true_nu_y = -999.;
		_true_nu_z = -999.;
		_reco_CCQE_E = -999.;
		_true_mu_len = -999.;
		_true_p_len = -999.;
		_reco_mu_len = -999.;
		_reco_p_len = -999.;
		_true_neutrons_E = 0.;
		_is_reco_vtx_near_true_vtx = false;
		_fndecay = 0;

		auto ev_mctrack = storage->get_data<event_mctrack>("mcreco");
		if (!ev_mctrack) {
			print(larlite::msg::kERROR, __FUNCTION__, Form("Did not find specified data product, mctrack!"));
			return false;
		}
		// if (ev_mctrack->size() != 2) {
		// 	print(larlite::msg::kWARNING, __FUNCTION__, Form("MCTrack size is not equal to 2... it equals %lu!", ev_mctrack->size()));
		// 	return false;
		// }
//Grab the MCTruth
		auto ev_mctruth = storage->get_data<event_mctruth>("generator");
		if (!ev_mctruth) {
			print(larlite::msg::kERROR, __FUNCTION__, Form("Did not find specified data product, mctruth!"));
			return false;
		}
		// Require exactly one neutrino interaction
		if (ev_mctruth->size() != 1) {
			print(larlite::msg::kINFO, __FUNCTION__, Form("ev_mctruth size is not 1!"));
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

		auto ev_track = storage->get_data<event_track>("pandoraNuPMA");
		if (!ev_track) {
			print(larlite::msg::kERROR, __FUNCTION__, Form("Did not find specified data product, track!"));
			return false;
		}

		auto ev_vtx = storage->get_data<event_vertex>("pmtrack");//pandoraNu
		if (!ev_vtx) {
			print(larlite::msg::kERROR, __FUNCTION__, Form("Did not find specified data product, vertex!"));
			return false;
		}

		// Sum up the true energy from neutrons in the true final state
		for (auto const& mypart : ev_mctruth->at(0).GetParticles())
			if (mypart.PdgCode() == 2112 && mypart.StatusCode() == 1 )
				_true_neutrons_E += mypart.Trajectory().front().E() - 0.938;

		_true_nu_E = ev_mctruth->at(0).GetNeutrino().Nu().Trajectory().front().E();
		auto const& nu_vtx = ev_mctruth->at(0).GetNeutrino().Nu().Trajectory().front().Position().Vect();
		_true_nu_x = nu_vtx.X();
		_true_nu_y = nu_vtx.Y();
		_true_nu_z = nu_vtx.Z();

		// Make a small 3cm geosphere around neutrino vertex
		_nu_sphere.Center(nu_vtx);
		_nu_sphere.Radius(3.);

		// Check if there is a reco vertex near true vertex
		for (auto const& myvtx : *ev_vtx) {
			auto geovtx = ::geoalgo::Vector(myvtx.X(), myvtx.Y(), myvtx.Z());
			if ( _nu_sphere.Contain(geovtx) )
				_is_reco_vtx_near_true_vtx = true;
		}


		// Store the muon and the proton mctracks for later investigation
		bool _found_muon = false;
		mctrack mu_mctrack;
		bool _found_proton = false;
		mctrack p_mctrack;

		// std::cout<<"Loop over mctracks"<<std::endl;
		for (auto const& mct : *ev_mctrack) {

			if (!mct.size()) continue;

			// Only count mctracks that stat within 1cm of the true neutrino vertex
			if ( !_nu_sphere.Contain(mct.front().Position().Vect()) )
				continue;

			_n_mc_tracks_fromvtx++;

			if (mct.PdgCode() == 13 ) {
				//if you've already found a muon mctrack, return false
				if (_found_muon) return false;
				_found_muon = true;
				mu_mctrack = mct;
				// std::cout<<"FRACITON: found muon. start total energy is "<<mu_mctrack.front().E()<<std::endl;
			}
			else if (mct.PdgCode() == 2212) {
				// if you've already found a proton mctrack, return false;
				if (_found_proton) return false;
				_found_proton = true;
				p_mctrack = mct;
				// std::cout<<"FRACTION: found proton. start total energy is "<<p_mctrack.front().E()<<std::endl;
			}
			else {

				// If other events exist in this event ... if you are using
				// MC_1mu_1pNn0else_Filter_MCTracks, this might happen for example when you have a high
				// energy proton from the neutrino that interacts to create a charged pion.
				// For now just throw this event out. I'm looking for CLEAN 1mu+1p events
				// with no crazy shit happening downstream.
				// Maybe 10% of the 1mu+1p events have this extra crazy shit

				// print(larlite::msg::kWARNING, "EventsContainedStudy",
				//       Form("Found a mctrack that wasn't p or mu! Its PDG is %d and its energy deposited is %0.1f MEV.",
				//            mct.PdgCode(), mct.front().E() - mct.back().E()));
				// if (std::fabs(mct.PdgCode()) == 211) {
				// std::cout << std::endl; std::cout << std::endl; std::cout << std::endl;
				// std::cout << "FRACTION: Here is the pion start point: ";
				// mct.front().Position().Vect().Print();
				// std::cout<<"FRACTION: Here is its deposited energy: "<<mct.front().E()-mct.back().E()<<std::endl;
				// throw std::exception();
				// std::cout<<"FRACTION: Here is its process string: "<<mct.Process()<<std::endl;
				// }


				return false;
			}
		}


		// If you haven't found a proton and a muon, something has gone wrong.
		if (!_found_proton || !_found_muon) {
			print(larlite::msg::kWARNING, "EventsContainedStudy",
			      Form("Either didn't find a proton, or didn't find a muon!"));
			return false;
		}

		_mu_contained = _fidvolBox.Contain(mu_mctrack.front().Position().Vect()) &&
		                _fidvolBox.Contain(mu_mctrack.back().Position().Vect());

		_p_contained = _fidvolBox.Contain(p_mctrack.front().Position().Vect()) &&
		               _fidvolBox.Contain(p_mctrack.back().Position().Vect());

		_true_mu_E = mu_mctrack.front().E() / 1000.;
		_true_p_E  = p_mctrack.front().E() / 1000.;
		_true_mu_len = (mu_mctrack.front().Position().Vect() - mu_mctrack.back().Position().Vect()).Mag();
		_true_p_len = (p_mctrack.front().Position().Vect() - p_mctrack.back().Position().Vect()).Mag();
		_h_mu_E_all->Fill(_true_mu_E);
		_h_nu_E_all->Fill(_true_nu_E);

		if (_mu_contained && _p_contained) {
			_h_mu_E_contained->Fill(_true_mu_E);
			_h_nu_E_contained->Fill(_true_nu_E);
			_evts_fully_contained++;
		}

		// Vector of all tracks starting or (ending & not starting) near MC vertex
		// pair is the (track, whether the track is flipped [ends at vtx instead of starts])
		std::vector<std::pair<larlite::track*, bool> > vtx_tracks;
		// std::cout<<"nu sphere center is "<<_nu_sphere.Center()<<std::endl;
		/// Loop over reco tracks and count the ones coming out of the vertex
		for (auto &trk : *ev_track) {
			if (trk.Length() > 500) {
				// std::cout<<"long-ish Reco track w00t! Length is "<<trk.Length()<<std::endl;
				// std::cout<<"track start is:";
				// trk.Vertex().Print();
				// std::cout<<"track end is: ";
				// trk.End().Print();
			}
			if (_nu_sphere.Contain(trk.Vertex())) {
				// std::cout<<"sphere contains track start! track length is "<<trk.Length()<<std::endl;
				_n_reco_tracks_fromvtx++;
				vtx_tracks.push_back(std::make_pair(&trk, false));
			}
			// If the end of the track is near vtx (and the start isn't, to ignore very short tracks)
			// that track is likely reconstructed backwards.... count those
			if (_nu_sphere.Contain(trk.End()) && !_nu_sphere.Contain(trk.Vertex())) {
				// std::cout<<"sphere contained track end! track length is "<<trk.Length()<<std::endl;
				_n_reco_tracks_endingnearvtx++;
				vtx_tracks.push_back(std::make_pair(&trk, true));
			}

		}

		/// If two reco tracks are coming from the vertex, take the longer one
		/// and assume it is a muon. The shorter is assumed proton.
		/// Use spline to compute the energy of each to store in tree.
		/// Hopefully mu energy + proton energy - proton mass = neutrino energy!
		double muon_mass_GEV = 0.106;
		if (vtx_tracks.size() == 2) {
			double trk0len = std::fabs( (vtx_tracks[0].first->Vertex() - vtx_tracks[0].first->End()).Mag() );
			double trk1len = std::fabs( (vtx_tracks[1].first->Vertex() - vtx_tracks[1].first->End()).Mag() );
			// larlite::track shorter_track = trk0len < trk1len ?
			//                                vtx_tracks[0].first : vtx_tracks[1].first;
			const larlite::track *longer_track = trk0len < trk1len ?
			                                     vtx_tracks[1].first : vtx_tracks[0].first;
			bool longer_track_flipped = trk0len < trk1len ?
			                            vtx_tracks[1].second : vtx_tracks[0].second;

			_reco_mu_len = std::max(trk0len, trk1len);
			_reco_p_len = std::min(trk0len, trk1len);
			_reco_mu_E = myspline->GetMuMomentum(trk0len > trk1len ? trk0len : trk1len) / 1000.;
			_reco_p_E =  myspline->GetPMomentum (trk0len > trk1len ? trk1len : trk0len) / 1000.;
			TVector3 mudir = longer_track_flipped ?
			                 (longer_track->Vertex() - longer_track->End()).Unit()
			                 : (longer_track->End() - longer_track->Vertex()).Unit();

			_reco_CCQE_E = mycalc->ComputeECCQE( (_reco_mu_E + muon_mass_GEV) * 1000., mudir, false);

			/// debug
			if (_true_nu_E > 1.368 && _true_nu_E < 1.37 &&
			        _reco_p_E + _reco_mu_E + 0.106 < 0.528 &&  _reco_p_E + _reco_mu_E + 0.106 > 0.526) {
				std::cout << "\n\n\n\n\nFound weird event!\n\n\n\n\n" << std::endl;
			}
		}
		/// If one reco track from vertex, assume it is muon.
		/// You can still get a CCQE energy (and spline muon energy) from just the muon.
		else if (vtx_tracks.size() == 1) {
			double trklen = std::fabs( (vtx_tracks[0].first->Vertex() - vtx_tracks[0].first->End()).Mag() );
			_reco_mu_len = trklen;
			_reco_mu_E = myspline->GetMuMomentum(trklen) / 1000.;
			TVector3 mudir = vtx_tracks[0].second ?
			                 (vtx_tracks[0].first->Vertex() - vtx_tracks[0].first->End()).Unit()
			                 : (vtx_tracks[0].first->End() - vtx_tracks[0].first->Vertex()).Unit();
			_reco_CCQE_E = mycalc->ComputeECCQE( (_reco_mu_E + muon_mass_GEV) * 1000., mudir, false);
		}
		_tree->Fill();


		return true;
	}

	bool EventsContainedStudy::finalize() {

		std::cout << "Events Contained Study! Input events = "
		          << _evts_pass_filter << ", fully contained events = "
		          << _evts_fully_contained << ", so the fraction of events contained is "
		          << (double)_evts_fully_contained / (double)_evts_pass_filter << std::endl;

		if (_fout) {
			_fout->cd();
			_h_mu_E_contained->Write();
			_h_mu_E_all->Write();
			_h_nu_E_contained->Write();
			_h_nu_E_all->Write();
			_tree->Write();
		}

		return true;
	}

}
#endif
