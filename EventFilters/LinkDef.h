//
// cint script to generate libraries
// Declaire namespace & classes you defined
// #pragma statement: order matters! Google it ;)
//

#ifdef __CINT__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class larlite::kaleko::MCTracksContainedFilter+;
#pragma link C++ class larlite::OnlyOneHighEMCTrackFilter+;

#pragma link C++ class larlite::CCQEMCFilter+;
#pragma link C++ class larlite::MBLEESelectionFilter+;
#pragma link C++ class larlite::MC_1eNpNn0else_Filter+;
#pragma link C++ class larlite::PerfectReco_1eNp0else_Filter+;

//#pragma link C++ class ertool::ERAnaLowEExcess+;
//#pragma link C++ class ertool::ERAnaSimpleAna+;

#pragma link C++ class larlite::MBLEESelectionFilter+;


#pragma link C++ class larlite::CleanPi0Filter+;
#pragma link C++ class larlite::MCShowersContainedFilter+;
#pragma link C++ class larlite::MC_CCnumu_MuonMCTrackContainedFilter+;
//ADD_NEW_CLASS ... do not change this line

#endif







