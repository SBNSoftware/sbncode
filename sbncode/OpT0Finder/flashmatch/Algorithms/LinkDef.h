//
// cint script to generate libraries
// Declaire namespace & classes you defined
// #pragma statement: order matters! Google it ;)
//

#ifdef __CINT__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

//#pragma link C++ class flashmatch::QWeightPoint+;
//#pragma link C++ class flashmatch::CommonAmps+;
#pragma link C++ class flashmatch::TimeCompatMatch+;
#pragma link C++ class flashmatch::MaxNPEWindow+;
#pragma link C++ class flashmatch::TimeRange+;
#pragma link C++ class flashmatch::TimeRangeSet+;
//#pragma link C++ class flashmatch::FilterArray+;
//#pragma link C++ class flashmatch::NPtFilter+;
#pragma link C++ class flashmatch::QLLMatch+;

#pragma link C++ class flashmatch::PhotonLibHypothesis+;
//#pragma link C++ class flashmatch::ChargeAnalytical+;
#pragma link C++ class flashmatch::LightPath+;
//ADD_NEW_CLASS ... do not change this line
#endif
