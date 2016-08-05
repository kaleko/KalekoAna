#ifndef LARLITE_KALEKONUITXN_CXX
#define LARLITE_KALEKONUITXN_CXX

#include "KalekoNuItxn.h"

namespace larlite {


	larlite::vertex     KalekoNuItxn::Vertex() const { return fVertex;    }
	std::vector<larlite::track>       KalekoNuItxn::Tracks() const { return fTracks;}
	std::vector<larlite::KalekoPID_t> KalekoNuItxn::PIDs()   const { return fPIDs;}
	std::vector<larlite::calorimetry> KalekoNuItxn::Calos()  const { return fCalos;}
	std::vector<larlite::track>       KalekoNuItxn::ExtraTracks() const { return fExtraTracks;}


	void KalekoNuItxn::printInfo() {
		std::cout << " --- KalekoNuItxn printing info --- " << std::endl;
		std::cout << "\t Vertex is at ("<<fVertex.X()<<", "<<fVertex.Y()<<", "<<fVertex.Z()<<")."<<std::endl;
		std::cout << "\t # tracks associated with vertex: " << fTracks.size() << std::endl;
		std::cout << "\t IDs of those tracks: ";
		for (auto const& trk : fTracks)
			std::cout << trk.ID() << ", ";
		std::cout << std::endl;

		std::cout << "\t # extra tracks: " << fExtraTracks.size() << std::endl;
		std::cout << "\t IDs of those tracks: ";
		for (auto const& trk : fExtraTracks)
			std::cout << trk.ID() << ", ";
		std::cout << std::endl;
	}
}

#endif
