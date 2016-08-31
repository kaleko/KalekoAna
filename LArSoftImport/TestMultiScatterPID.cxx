#ifndef TESTMULTISCATTERPID_CXX
#define TESTMULTISCATTERPID_CXX

#include "TestMultiScatterPID.h"
#include "DataFormat/mctrack.h"

namespace larlite {

  bool TestMultiScatterPID::initialize() {

    _tmc = kaleko::TrackMomentumCalculator();

    if(!_ana_tree){
      _ana_tree = new TTree("ana_tree","ana_tree");
      _ana_tree->Branch("llhd",&_llhd,"llhd/D");
      _ana_tree->Branch("len",&_length,"len/D");

    }

    return true;
  }
  
  bool TestMultiScatterPID::analyze(storage_manager* storage) {
  
    //Get the MCTracks
    auto ev_mctrack = storage->get_data<event_mctrack>("mcreco");    
    if(!ev_mctrack) {
      print(larlite::msg::kERROR,__FUNCTION__,Form("Did not find specified data product, mctrack!"));
      return false;
    }  

    //Loop over the MCTracks
    for(auto const& mct : *ev_mctrack){
      _llhd = _tmc.GetMuMultiScatterLLHD2(mct);
      _length = (mct.End().Position().Vect() - mct.Start().Position().Vect()).Mag();
      _ana_tree->Fill();
    }  


    return true;
  }

  bool TestMultiScatterPID::finalize() {

    if(_fout) { _fout->cd(); _ana_tree->Write(); }
    
    else 
      print(larlite::msg::kERROR,__FUNCTION__,"Did not find an output file pointer!!! File not opened?");
    

    if(_ana_tree)
      delete _ana_tree;

    return true;
  }

}
#endif
