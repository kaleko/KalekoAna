#ifndef LARLITE_KALEKONUITXN_CXX
#define LARLITE_KALEKONUITXN_CXX

#include "KalekoNuItxn.h"

namespace larlite {


  larlite::vertex     KalekoNuItxn::Vertex() const { return fVertex;    }
  std::vector<larlite::track>       KalekoNuItxn::Tracks() const { return fTracks;}
  std::vector<larlite::KalekoPID_t> KalekoNuItxn::PIDs()   const { return fPIDs;}
  std::vector<larlite::calorimetry> KalekoNuItxn::Calos()  const { return fCalos;}

}

#endif
