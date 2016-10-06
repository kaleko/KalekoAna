#ifndef LARLITE_TESTTRKHITASSN_CXX
#define LARLITE_TESTTRKHITASSN_CXX

#include "TestTrkHitAssn.h"
#include "DataFormat/track.h"
#include "DataFormat/hit.h"
#include "DataFormat/mctrack.h"

namespace larlite {

    bool TestTrkHitAssn::initialize() {

        _tree = new TTree("tree","tree");
        _tree->Branch("sum_hit_ADC",&_sum_hit_ADC,"sum_hit_ADC/D");
        _tree->Branch("mct_depE",&_mct_depE,"mct_depE/D");
        _tree->Branch("n_asstd_hits",&_n_asstd_hits,"n_asstd_hits/I");
        _tree->Branch("n_tothits",&_n_tothits,"n_tothits/I");

        return true;
    }

    bool TestTrkHitAssn::analyze(storage_manager* storage) {


        // auto ev_mctrack = storage->get_data<event_mctrack>("mcreco");
        // if (!ev_mctrack) {
        //     print(larlite::msg::kERROR, __FUNCTION__, Form("Did not find specified data product, mctrack!"));
        //     return false;
        // }
        // if (ev_mctrack->size() != 1) {
        //     print(larlite::msg::kERROR, __FUNCTION__, Form("# of MCTracks in this event is not 1. Run this on single muons!"));
        //     return false;
        // }

        auto ev_track = storage->get_data<event_track>("trackkalmanhit");
        if (!ev_track) {
            print(larlite::msg::kERROR, __FUNCTION__, Form("Did not find specified data product, track!"));
            return false;
        }
        if (!ev_track->size()) {
            print(larlite::msg::kERROR, __FUNCTION__, Form("Zero reconstructed tracks in this event!"));
            return false;
        }

        event_hit* ev_hit = nullptr;
        auto const& ass_hit_v = storage->find_one_ass(ev_track->id(),
                                ev_hit);

        if (!ev_hit) {
            std::cerr << "no event_hit!" << std::endl;
            return false;
        }

        if (ev_hit->empty()) {
            std::cout << "ev_hit empty" << std::endl;
            return false;
        }

        _n_tothits = ev_hit->size();
        // std::cout << "ev_hit size is " << _n_tothits << std::endl;
        // std::cout << "ev_hit id is std::pair<"
        //           << ev_hit->id().first << ","
        //           << ev_hit->id().second.c_str()
        //           << ">" << std::endl;


        size_t reco_track_idx = 0; //since there is only one track, otherwise choose the index in ev_track you care about

        // get the hits associated with this track
        _n_asstd_hits = ass_hit_v.at(reco_track_idx).size();
        _sum_hit_ADC = 0;
        // std::cout << "ass_hit_v has size " << ass_hit_v.size() << std::endl;
        // std::cout << "ass_hit_v.at(0) has size " << _n_asstd_hits << std::endl;

        for (auto const& asstd_hit_idx : ass_hit_v.at(reco_track_idx))
        {
            const larlite::hit ihit = ev_hit->at(asstd_hit_idx);
            _sum_hit_ADC += ihit.SummedADC();
        }

        // _mct_depE = ev_mctrack->at(0).front().E() - ev_mctrack->at(0).back().E();

        // std::cout << "MCTrack truth start kinetic energy is roughly : " << ev_mctrack->at(0).front().E() - 106. << std::endl;
        // std::cout << "MCTrack deposited energy is : " << ev_mctrack->at(0).front().E() - ev_mctrack->at(0).back().E() << " MEV." << std::endl;
        // std::cout << "Summed ADC of associated hits is : "<<_sum_hit_ADC<<std::endl;

        _tree->Fill();
        
        return true;
    }

    bool TestTrkHitAssn::finalize() {

        if(_fout){
            _fout->cd();
            _tree->Write();
        }

        return true;
    }

}
#endif
