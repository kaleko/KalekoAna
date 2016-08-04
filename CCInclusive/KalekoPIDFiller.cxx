#ifndef LARLITE_KALEKOPIDFILLER_CXX
#define LARLITE_KALEKOPIDFILLER_CXX

#include "KalekoPIDFiller.h"

namespace larlite {

    bool KalekoPIDFiller::fillKalekoPIDs_simple(KalekoNuItxn &kaleko_itxn) {

        // Loop over associated track lengths to find the longest one.
        // Set the longest one's PID value to muon.
        double longest_trklen = -999.;
        size_t longest_itrk   = 0;
        for (size_t itrk = 0; itrk < kaleko_itxn.Tracks().size(); ++itrk) {
            double itrklen = kaleko_itxn.Tracks().at(itrk).Length();
            if (itrklen > longest_trklen) {
                longest_itrk = itrk;
                longest_trklen = itrklen;
            }
        }
        kaleko_itxn.ChangePID(longest_itrk, kKalekoMuon);

        // Loop over the other associated tracks and if length is > 20cm call it pion,
        // else call it proton (this is SUPER hacky and a very quick effort, should
        // be improved with calorimetry etc. stuff soon)
        for (size_t itrk = 0; itrk < kaleko_itxn.Tracks().size(); ++itrk) {
            if (itrk == longest_itrk) continue;
            double itrklen = kaleko_itxn.Tracks().at(itrk).Length();
            if (itrklen > 20.) kaleko_itxn.ChangePID(itrk, kKalekoChargedPion);
            else kaleko_itxn.ChangePID(itrk, kKalekoProton);
        }

        return true;
    }


    bool KalekoPIDFiller::fillKalekoPIDs(KalekoNuItxn &kaleko_itxn) {

        // Loop over associated track lengths to find the longest one.
        // Set the longest one's PID value to muon.
        double longest_trklen = -999.;
        size_t longest_itrk   = 0;
        for (size_t itrk = 0; itrk < kaleko_itxn.Tracks().size(); ++itrk) {
            double itrklen = kaleko_itxn.Tracks().at(itrk).Length();
            if (itrklen > longest_trklen) {
                longest_itrk = itrk;
                longest_trklen = itrklen;
            }
        }
        kaleko_itxn.ChangePID(longest_itrk, kKalekoMuon);

        // Loop over the other associated tracks use calorimetry cut on 4
        // (average dEdx over the track)
        for (size_t itrk = 0; itrk < kaleko_itxn.Tracks().size(); ++itrk) {
            if (itrk == longest_itrk) continue;
            auto thecalo = kaleko_itxn.Calos().at(itrk);

            // Compute average calo
            double avg_dedx = 0;
            for (size_t j = 0; j < thecalo.dEdx().size(); ++j)
                avg_dedx += thecalo.dEdx().at(j);
            avg_dedx /= thecalo.dEdx().size();

            // if avg dedx is > 4, it's proton
            if (avg_dedx >= 4.) kaleko_itxn.ChangePID(itrk, kKalekoProton);
            // if avg dedx is between 1.5 and 4, it's charged pion
            else if (avg_dedx > 1.5 && avg_dedx < 4.) kaleko_itxn.ChangePID(itrk, kKalekoChargedPion);
            // if avg dedx is less than 1.5, dedx is not well reconstructed and make decision based on length
            else {
                auto thetrack = kaleko_itxn.Tracks().at(itrk);
                double thetracklen = (thetrack.End() - thetrack.Vertex()).Mag();
                if (thetracklen > 20)
                    kaleko_itxn.ChangePID(itrk, kKalekoChargedPion);
                else
                    kaleko_itxn.ChangePID(itrk, kKalekoProton);
            }
        }

        return true;
    }

}
#endif
