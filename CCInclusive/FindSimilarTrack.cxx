#ifndef LARLITE_FINDSIMILARTRACK_CXX
#define LARLITE_FINDSIMILARTRACK_CXX

#include "FindSimilarTrack.h"

namespace larlite {

    larlite::track FindSimilarTrack::findSimilarTrack(const larlite::track &trk, const larlite::event_track &ev_track) {

        TVector3 trk_dir = trk.End() - trk.Vertex();

        size_t best_trk_idx = 9999;
        double max_dot_prod = -9999.;
        for(size_t i = 0; i < ev_track.size(); ++i){
            auto const & itrk = ev_track.at(i);
            TVector3 itrk_dir = itrk.End() - itrk.Vertex();
            bool flip_itrk = (itrk.Vertex() - trk.Vertex()).Mag2() > (itrk.End() - trk.Vertex()).Mag2();
            double idotprod = flip_itrk ? (itrk_dir * -1.).Unit().Dot(trk_dir.Unit()) : itrk_dir.Unit().Dot(trk_dir.Unit());
            double idist = std::min( (itrk.Vertex() - trk.Vertex()).Mag(), (itrk.End() - trk.Vertex()).Mag() );
            bool close_enough = idist < 5.;
            if( idotprod > max_dot_prod && close_enough){
                max_dot_prod = idotprod;
                best_trk_idx = i;
            }
        }

        /// If no matching track was found, throw exception
        if(best_trk_idx == 9999)
            throw std::exception();
        
        return ev_track.at(best_trk_idx);

    }

}
#endif
