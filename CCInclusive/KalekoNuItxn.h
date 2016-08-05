/**
 * \file KalekoNuItxn.h
 *
 * \ingroup
 *
 * \brief Class def header for KalekoNuItxn data container
 *
 * @author Kaleko
 */

/** \addtogroup awfeawf

@{*/

#ifndef LARLITE_KALEKONUITXN_H
#define LARLITE_KALEKONUITXN_H

#include "CCInclusiveConstants.h"
#include "DataFormat/vertex.h"
#include "DataFormat/track.h"
#include "DataFormat/calorimetry.h"

namespace larlite {

  class KalekoNuItxn {

  public:

    /// Default constructor
    KalekoNuItxn() {
      fTracks.clear();
      fPIDs.clear();
      fCalos.clear();
    }

    /// Default destructor
    virtual ~KalekoNuItxn() {}

    //--- Getters ---//
    larlite::vertex                   Vertex() const;
    std::vector<larlite::track>       Tracks() const;
    std::vector<larlite::KalekoPID_t> PIDs() const;
    std::vector<larlite::calorimetry> Calos() const;
    std::vector<larlite::track>       ExtraTracks() const;

    //--- Setters ---//
    void Vertex   ( larlite::vertex v )       { fVertex = v; }
    void AddTrack ( larlite::track t  )       { fTracks.push_back(t); }
    void AddPID   ( larlite::KalekoPID_t p  ) { fPIDs.push_back(p); }
    void AddCalo  ( larlite::calorimetry c  ) { fCalos.push_back(c); }
    void ChangePID( size_t PID_idx, larlite::KalekoPID_t p ) { fPIDs.at(PID_idx) = p; }
    void AddExtraTrack( larlite::track t )    { fExtraTracks.push_back(t); }


    void printInfo();
    
  protected:
//      typedef std::pair<larlite::vertex, std::vector< std::pair<larlite::KalekoPID_t, larlite::track> > > CCNuItxn_t;

    larlite::vertex fVertex;  /// The reconstructed vertex

    std::vector<larlite::track>       fTracks; /// A vector (length n) of n tracks associated with vertex
    std::vector<larlite::KalekoPID_t> fPIDs;   /// A vector (length n) of PID decisions associated per vertex track
    std::vector<larlite::calorimetry> fCalos;  /// A vector (length n) of calo objects associated per vertex track

    std::vector<larlite::track>       fExtraTracks; /// Possible additional tracks NOT associated with vertex
                                                    /// E.G. a pion from vtx stops, decays into a muon, which is a separate track

  };

}
#endif

/** @} */ // end of doxygen group
