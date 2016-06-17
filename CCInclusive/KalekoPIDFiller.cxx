#ifndef LARLITE_KALEKOPIDFILLER_CXX
#define LARLITE_KALEKOPIDFILLER_CXX

#include "KalekoPIDFiller.h"

namespace larlite {

    bool KalekoPIDFiller::fillKalekoPIDs(KalekoNuItxn_t &kaleko_itxn) {

        // Loop over associated track lengths to find the longest one.
        // Set the longest one's PID value to muon.
        double longest_trklen = -999.;
        size_t longest_itrk   = 0;
        for (size_t itrk = 0; itrk < kaleko_itxn.second.size(); ++itrk){
            double itrklen = kaleko_itxn.second.at(itrk).second.Length();
            if (itrklen > longest_trklen){
                longest_itrk = itrk;
                longest_trklen = itrklen;
            }
        }
        kaleko_itxn.second.at(longest_itrk).first = kKalekoMuon;

        // Loop over the other associated tracks and if length is > 20cm call it pion,
        // else call it proton (this is SUPER hacky and a very quick effort, should
        // be improved with calorimetry etc. stuff soon)
        for (size_t itrk = 0; itrk < kaleko_itxn.second.size(); ++itrk){
            if(itrk == longest_itrk) continue;
            double itrklen = kaleko_itxn.second.at(itrk).second.Length();
            if(itrklen > 20.) kaleko_itxn.second.at(itrk).first = kKalekoChargedPion;
            else kaleko_itxn.second.at(itrk).first = kKalekoProton;
        }
        
        return true;
    }

}
#endif
