#ifndef MCTRACKSCONTAINEDFILTER_CXX
#define MCTRACKSCONTAINEDFILTER_CXX

#include "MCTracksContainedFilter.h"
#include "FidVolBox.h"

namespace larlite {

  bool MCTracksContainedFilter::initialize() {

    _myGeoAABox = FidVolBox();

    // //Set DistToBoxWall's "box" to be TPC 
    // _myGeoAABox.Min( 0,
		  //    -(::larutil::Geometry::GetME()->DetHalfHeight()),
		  //    0);
    
    // _myGeoAABox.Max( 2*(::larutil::Geometry::GetME()->DetHalfWidth()),
		  //    ::larutil::Geometry::GetME()->DetHalfHeight(),
		  //    ::larutil::Geometry::GetME()->DetLength());
    
    kept_events = 0;
    total_events = 0;

    return true;
  }
  
  bool MCTracksContainedFilter::analyze(storage_manager* storage) {

    //Grab the MCTracks
    auto ev_mctrack = storage->get_data<event_mctrack>("mcreco");    
    if(!ev_mctrack) {
      print(larlite::msg::kERROR,__FUNCTION__,Form("Did not find specified data product, mctrack!"));
      return false;
    }

    total_events++;

    //Make sure all MC tracks are fully contained
    for (auto &track : *ev_mctrack){
      if(track.Origin() != 1) continue;
      if ( !isFullyContained(track) ) return false;
    }
    
    //Looks like you want to keep this event.
    kept_events++;
    return true;

  }

  bool MCTracksContainedFilter::finalize() {

    std::cout<<"MCTracksContainedFilter: "<<total_events<<" total events, "<<kept_events<<" kept events."<<std::endl;
    
    return true;
  }

  bool MCTracksContainedFilter::isFullyContained(mctrack const & mytrack){

    //Make sure track MC start point and MC end point are in active volume

    //   Point_t mypoint(mytrack.Start().Position());

    if(_myGeoAABox.Contain(mytrack.Start().Position()) > 0 &&
       _myGeoAABox.Contain(mytrack.End().Position()) > 0)

      return true;

    return false;
  }
}
#endif
