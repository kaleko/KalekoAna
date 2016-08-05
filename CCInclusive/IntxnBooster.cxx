#ifndef LARLITE_INTXNBOOSTER_CXX
#define LARLITE_INTXNBOOSTER_CXX

#include "IntxnBooster.h"
#include <set>

namespace larlite {

	void IntxnBooster::BoostIntxn( larlite::KalekoNuItxn & itxn, const larlite::event_track * ev_track){

		// Loop through all tracks in event_track that are already associated with the vertex.
		// For each, loop over all OTHER tracks in the event and decide if any should
		// be added as an "Additional track" to the interaction
		
		// First get a list of all vertex-track indices 
		std::set<size_t> vtx_trk_ids;
		for(auto const vtxtrk : itxn.Tracks())
			vtx_trk_ids.insert(vtxtrk.ID());

		// Loop over vertex tracks and compare each to all other tracks
		for(auto const vtxtrk : itxn.Tracks()){
			for (auto const &trk : *ev_track){
				// If this track is already in the list of vertex tracks, skip it
				if(vtx_trk_ids.find(trk.ID()) != vtx_trk_ids.end()) continue;
			}

		}
		return;
	}

}//end namespace larlite

#endif
