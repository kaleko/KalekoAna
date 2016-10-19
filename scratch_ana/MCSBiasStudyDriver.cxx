#ifndef LARLITE_MCSBIASSTUDYDRIVER_CXX
#define LARLITE_MCSBIASSTUDYDRIVER_CXX

#include "MCSBiasStudyDriver.h"

namespace larlite {

    bool MCSBiasStudyDriver::initialize() {

        if (!_module)
            _module = new MCSBiasStudy();

        return true;
    }

    bool MCSBiasStudyDriver::analyze(storage_manager* storage) {

        auto ev_track = storage->get_data<event_track>("pandoraNuPMA");
        if (!ev_track) {
            print(larlite::msg::kERROR, __FUNCTION__, Form("Did not find specified data product, track!"));
            return false;
        }
        if (!ev_track->size())
            return false;
        // if (ev_track->size() != 1)
        //     return false;

        // auto ev_mctrack = storage->get_data<event_mctrack>("mcreco");
        // if (!ev_mctrack) {
        //     print(larlite::msg::kERROR, __FUNCTION__, Form("Did not find specified data product, mctrack!"));
        //     return false;
        // }
        // if (!ev_mctrack->size())
        //     return false;
        // if (ev_mctrack->size() != 1)
        //     return false;


        //  if ( !_fidvol.Contain(ev_mctrack->at(0).front().Position().Vect()) || !_fidvol.Contain(ev_mctrack->at(0).back().Position().Vect()) ) {
        //     // print(larlite::msg::kERROR, __FUNCTION__, Form("MC track not fully contained in fiducial volume."));
        //     return false;
        // }
        // if ( !_fidvol.Contain(ev_track->at(0).Vertex()) || !_fidvol.Contain(ev_track->at(0).End()) ) {
        //     // print(larlite::msg::kERROR, __FUNCTION__, Form("Reco track not fully contained in fiducial volume."));
        //     return false;
        // }

        // if ((ev_track->at(0).End() - ev_track->at(0).Vertex()).Mag() < 100.) return false;

        for (auto const &track : *ev_track)
            // if (_fidvol.Contain(track.Vertex()) && _fidvol.Contain(track.End()))
            // looking at stopping cosmic muons right now
            if(track.Vertex().Y() > 100 && _fidvol.Contain(track.End())){
                // std::cout<<"analyzing track!!!!"<<std::endl;
                _module->AnalyzeTrack(track);
            }

        // for (auto const &mct : *ev_mctrack){
        //     if(mct.size() < 3) continue;
        //     // if (_fidvol.Contain(track.Vertex()) && _fidvol.Contain(track.End()))
        //     // looking at stopping cosmic muons right now
        //     if(mct.front().Y() > 100 && _fidvol.Contain(mct.back().Position().Vect()))
        //         _module->AnalyzeTrack(mct);
        // }
         // _module->AnalyzeTrack(ev_track->at(0));


        // for (auto const &mct : *ev_mctrack){
        //     if(mct.size() < 3) continue;
        //     // if (_fidvol.Contain(track.Vertex()) && _fidvol.Contain(track.End()))
        //     // looking at stopping cosmic muons right now
        //     if(mct.front().Position().Vect().Y() > 100 && _fidvol.Contain(mct.back().Position().Vect()))
        //         _module->AnalyzeTrack(mct);
        // }




        return true;
    }

    bool MCSBiasStudyDriver::finalize() {

        if (_fout) {
            _fout->cd();
            _module->GetTree()->Write();
            _module->GetSegTree()->Write();
            _module->GetTMCTree()->Write();
        }

        return true;
    }

}
#endif
