#ifndef LARLITE_BRYCECHECK_CXX
#define LARLITE_BRYCECHECK_CXX

#include "BryceCheck.h"

namespace larlite {

	bool BryceCheck::initialize() {

		_fidvolBox = FidVolBox();

			//Box here is TPC
		_TPCBox.Min( 0 + 2.,
		                -(::larutil::Geometry::GetME()->DetHalfHeight()) + 2.,
		                0 + 2.);

		_TPCBox.Max( 2 * (::larutil::Geometry::GetME()->DetHalfWidth()) - 2.,
		                ::larutil::Geometry::GetME()->DetHalfHeight() - 2.,
		                ::larutil::Geometry::GetME()->DetLength() - 2.);


		evt_counter = 0;
		CC_counter = 0;
		in_fidvol = 0;
		is_numu = 0;
		one_muon = 0;
		the_muon_exits_TPC = 0;
		the_muon_longenough = 0;

		return true;
	}

	bool BryceCheck::analyze(storage_manager* storage) {
		evt_counter++;
		auto ev_mctruth = storage->get_data<event_mctruth>("generator");
		if (!ev_mctruth) {
			print(larlite::msg::kERROR, __FUNCTION__, Form("Did not find specified data product, mctruth!"));
			return false;
		}
		if (ev_mctruth->size() != 1) {
			// print(larlite::msg::kERROR, __FUNCTION__, Form("MCTruth size is not equal to 1... it equals %lu!", ev_mctruth->size()));
			return false;
		}


		auto const &mctruth = ev_mctruth->at(0);

		//Enforce CC interaction channel
		if ( mctruth.GetNeutrino().CCNC() != 0 ) return false;
		CC_counter++;

		// If neutrino interacts outside of fiducial volume, skip event
		auto const nu_vtx = mctruth.GetNeutrino().Nu().Trajectory().back().Position().Vect();
		if (!_fidvolBox.Contain(nu_vtx)) return false;
		in_fidvol++;

		// If neutrino is not a numu, skip event
		if (mctruth.GetNeutrino().Nu().PdgCode() != 14) return false;
		is_numu++;


		// Now we make sure the muon mctrack from the interaction is EXITING the fiducial volume
		// Grab the MCTracks
		auto ev_mctrack = storage->get_data<event_mctrack>("mcreco");
		if (!ev_mctrack) {
			print(larlite::msg::kERROR, __FUNCTION__, Form("Did not find specified data product, mctrack!"));
			return false;
		}
		if (!ev_mctrack->size()) {
			print(larlite::msg::kERROR, __FUNCTION__, Form("Zero mctracks in this event?! Skipping!"));
			return false;
		}


		// Find the mctrack that is a muon starting from the neutrino interaction
		// We also require it is 1m long at least... might want to remove this if you are using this
		// filter for some reason other than a multiple coloumb scattering analysis which requires 1m
		size_t n_found_muons = 0;
		mctrack the_muon;
		for (auto const& mct : *ev_mctrack) {

			// Sometimes mctracks have zero size. No idea why. Skip them.
			if ( mct.size() < 3 ) continue;

			// origin == 1 means comes from neutrino interaction (IE not cosmic)
			if (mct.Origin() != 1 ) continue;

			// MCTrack has to be truly a muon
			if (mct.PdgCode() != 13 ) continue;

			//MCTrack should start VERY CLOSE to nu vtx
			// (we're talking like 10^-28 or smaller ... below floating point precision)
			if ( (mct.front().Position().Vect() - nu_vtx).Mag2() > 0.0001) continue;


			// Require the muon exits the actual TPC
			// this is done by requiring the energy of the last step in the MCTrack std::vector
			// (which corresponds to the last energy deposition point inside of the TPC)
			// has more energy than the rest mass energy (105.658 MeV)
			// if ((mct.back().E() < 106.)) continue;

			// At this point you have found a viable muon!
			n_found_muons++;
			the_muon = mct;

		}

		// If you didn't find any viable muons, skip the event
		if (!n_found_muons) return false;


		// If somehow you found more than one muon, something has gone wrong. Skip the event.
		if (n_found_muons > 1) {
			print(larlite::msg::kWARNING, __FUNCTION__, Form("More than one viable muon in the event?!"));
			return false;
		}
		one_muon++;


		// Enforce the muon is exiting the TPC.
		if (_TPCBox.Contain(the_muon.back().Position().Vect())) return false;
		the_muon_exits_TPC++;

		// Require the muon is at least one meter in length (in the TPC).
		if ((the_muon.back().Position().Vect() - the_muon.front().Position().Vect()).Mag() < 100.) return false;
		the_muon_longenough++;

		return true;
	}

	bool BryceCheck::finalize() {

		std::cout << "FINALIZING:" << std::endl;
		std::cout << "evt_counter is " << evt_counter << std::endl;
		std::cout << "CC_counter is " << CC_counter << std::endl;
		std::cout << "in_fidvol is " << in_fidvol << std::endl;
		std::cout << "is_numu is " << is_numu << std::endl;
		std::cout << "one_muon is " << one_muon << std::endl;
		std::cout << "the_muon_exits_TPC is " << the_muon_exits_TPC << std::endl;
		std::cout << "the_muon_longenough is " << the_muon_longenough << std::endl;
		std::cout << "DONE FINALIZING" << std::endl;

		return true;
	}

}
#endif
