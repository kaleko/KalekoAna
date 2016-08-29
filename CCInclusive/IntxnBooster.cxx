#ifndef LARLITE_INTXNBOOSTER_CXX
#define LARLITE_INTXNBOOSTER_CXX

#include "IntxnBooster.h"
#include <set>

namespace larlite {

	void IntxnBooster::BoostIntxn( larlite::KalekoNuItxn & itxn, const larlite::event_track * ev_track) {

		// Loop through all tracks in event_track that are already associated with the vertex.
		// For each, loop over all OTHER tracks in the event and decide if any should
		// be added as an "Additional track" to the interaction

		// First get a list of all vertex-track indices
		std::set<size_t> vtx_trk_ids;
		for (auto const vtxtrk : itxn.Tracks())
			vtx_trk_ids.insert(vtxtrk.ID());

		auto const &geovtx = ::geoalgo::Vector(itxn.Vertex().X(),
		                                       itxn.Vertex().Y(),
		                                       itxn.Vertex().Z());

		// Loop over vertex tracks and compare each to all other tracks
		// Keep track of track indices already added to other tracks
		std::set<size_t> added_trk_ids;
		for (auto const vtxtrk : itxn.Tracks()) {

			auto vtxtrk_end =  ::geoalgo::Vector(vtxtrk.Vertex()).SqDist(geovtx) <
			                   ::geoalgo::Vector(vtxtrk.End()).SqDist(geovtx) ?
			                   ::geoalgo::Vector(vtxtrk.End()) :
			                   ::geoalgo::Vector(vtxtrk.Vertex());

			for (auto const &trk : *ev_track) {

				// If this track is already in the list of vertex tracks, skip it
				if (vtx_trk_ids.find(trk.ID()) != vtx_trk_ids.end()) continue;
				// If this track has already been added to the itxn as an "Extra track", skip it
				if (added_trk_ids.find(trk.ID()) != added_trk_ids.end()) continue;

				auto trk_start = ::geoalgo::Vector(trk.Vertex());
				auto trk_end = ::geoalgo::Vector(trk.End());

				if (trk_start.Dist(vtxtrk_end) < 2. || trk_end.Dist(vtxtrk_end) < 2.) {
					added_trk_ids.insert(trk.ID());
					itxn.AddExtraTrack(trk);
				}

			}

		}
		return;
	}


	void IntxnBooster::AddKalmanTracks( larlite::KalekoNuItxn & itxn, const larlite::event_track * ev_track) {

		// Loop through all kalman tracks in event_track.
		// Pick ones that pass close to the 3D vertex.
		// Make sure they aren't in the same 3D direction as any existing tracks in the interaction
		// Add them to the interaction
		
		auto const &tvecvtx = TVector3(itxn.Vertex().X(),
		                               itxn.Vertex().Y(),
		                               itxn.Vertex().Z());

		// Loop over the kalman tracks
		for (auto const &kaltrk : *ev_track) {

			bool add_this_kaltrk = true;

			// Compute how close this track comes to the vertex
			double mindist2 = 99999.;
			for (size_t i = 0; i < kaltrk.NumberTrajectoryPoints(); ++i) {
				double dist2 = (kaltrk.LocationAtPoint(i) - tvecvtx).Mag2();
				if ( dist2 < mindist2 ) mindist2 = dist2;
			}

			if (mindist2 > 5.) continue;

			// If kalman track comes within 5cm of vertex
			// Loop over tracks associated with the vertex to make sure
			// this track isn't pointed in the same direction
			TVector3 kaltrk_dir = (kaltrk.End() - kaltrk.Vertex()).Unit();
			for (auto const vtxtrk : itxn.Tracks()) {
				TVector3 vtxtrk_dir = (vtxtrk.End() - vtxtrk.Vertex()).Unit();
				if (fabs(kaltrk_dir.Dot(vtxtrk_dir)) > 0.9) add_this_kaltrk = false;
			}

			if(add_this_kaltrk){
				itxn.AddTrack(kaltrk);
				itxn.AddCalo(larlite::calorimetry());
				itxn.AddPID(kKALEKO_PID_MAX);
			}
		}// End loop over kalman tracks
		
		return;

	}
}//end namespace larlite

#endif
