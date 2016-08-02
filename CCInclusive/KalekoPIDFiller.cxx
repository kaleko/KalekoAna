#ifndef LARLITE_KALEKOPIDFILLER_CXX
#define LARLITE_KALEKOPIDFILLER_CXX

#include "KalekoPIDFiller.h"

namespace larlite {

    bool KalekoPIDFiller::fillKalekoPIDs(KalekoNuItxn &kaleko_itxn) {

        // Loop over associated track lengths to find the longest one.
        // Set the longest one's PID value to muon.
        double longest_trklen = -999.;
        size_t longest_itrk   = 0;
        for (size_t itrk = 0; itrk < kaleko_itxn.Tracks().size(); ++itrk){
            double itrklen = kaleko_itxn.Tracks().at(itrk).Length();
            if (itrklen > longest_trklen){
                longest_itrk = itrk;
                longest_trklen = itrklen;
            }
        }
        kaleko_itxn.ChangePID(longest_itrk, kKalekoMuon);

        // Loop over the other associated tracks and if length is > 20cm call it pion,
        // else call it proton (this is SUPER hacky and a very quick effort, should
        // be improved with calorimetry etc. stuff soon)
        for (size_t itrk = 0; itrk < kaleko_itxn.Tracks().size(); ++itrk){
            if(itrk == longest_itrk) continue;
            double itrklen = kaleko_itxn.Tracks().at(itrk).Length();
            if(itrklen > 20.) kaleko_itxn.ChangePID(itrk,kKalekoChargedPion);
            else kaleko_itxn.ChangePID(itrk,kKalekoProton);
        }
        
        return true;
    }

}
#endif
