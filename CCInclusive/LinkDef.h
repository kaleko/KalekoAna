 //
// cint script to generate libraries
// Declaire namespace & classes you defined
// #pragma statement: order matters! Google it ;)
//

#ifdef __CINT__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ enum larlite::InputFileType_t+;

#pragma link C++ class larlite::EfficiencyStudy+;
#pragma link C++ class larlite::TrackMomentumSplines+;
#pragma link C++ class larlite::NuEnergyCalc+;
#pragma link C++ class larlite::TestCCQECalc+;
#pragma link C++ class larlite::TestSpline+;
#pragma link C++ class larlite::EventSelector+;
#pragma link C++ class larlite::TestKinematicNuReco+;
#pragma link C++ class larlite::XiaoEventAna+;
#pragma link C++ class larlite::MC_1mu1pNn0else_Filter+;
#pragma link C++ class larlite::MC_1mu1pNn0else_Filter_MCTracks+;
#pragma link C++ class larlite::TestMCTruth+;
#pragma link C++ class larlite::EventsContainedStudy+;
#pragma link C++ class larlite::NuMuCCFilter+;
#pragma link C++ class larlite::XiaoNuFinder+;

#pragma link C++ class larlite::TestPPiCalo+;
#pragma link C++ class larlite::QuickMikePlot+;
#pragma link C++ class larlite::TestTrkHitAssn+;
#pragma link C++ class larlite::KalekoPIDFiller+;
#pragma link C++ class larlite::FidVolBox+;
#pragma link C++ class larlite::KalekoTrackStitcher+;

#pragma link C++ class larlite::KalekoNuItxn+;
#pragma link C++ class larlite::FluxAna+;
//ADD_NEW_CLASS ... do not change this line
#endif




















