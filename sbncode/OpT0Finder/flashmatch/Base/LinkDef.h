//
// cint script to generate libraries
// Declaire namespace & classes you defined
// #pragma statement: order matters! Google it ;)
//

#ifdef __CLING__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ namespace flashmatch;
#pragma link C++ namespace flashmatch::msg;
#pragma link C++ class flashmatch::DetectorSpecs+;
#pragma link C++ class flashmatch::QPoint_t+;
#pragma link C++ class flashmatch::QCluster_t+;
#pragma link C++ class flashmatch::Flash_t+;
#pragma link C++ class flashmatch::QPoint_t+;
#pragma link C++ class flashmatch::QCluster_t+;
#pragma link C++ class flashmatch::FlashMatch_t+;
#pragma link C++ class std::vector<flashmatch::Flash_t>+;
#pragma link C++ class std::vector<flashmatch::QPoint_t>+;
#pragma link C++ class std::vector<flashmatch::QCluster_t>+;
#pragma link C++ class std::vector<flashmatch::FlashMatch_t>+;
#pragma link C++ class flashmatch::QClusterArray_t+;
#pragma link C++ class flashmatch::FlashArray_t+;
#pragma link C++ class flashmatch::FlashMatchManager+;
#pragma link C++ class flashmatch::BaseAlgorithm+;
#pragma link C++ class flashmatch::BaseProhibitAlgo+;
#pragma link C++ class flashmatch::BaseTPCFilter+;
#pragma link C++ class flashmatch::BaseFlashFilter+;
#pragma link C++ class flashmatch::BaseFlashMatch+;
#pragma link C++ class flashmatch::BaseFlashHypothesis+;
//ADD_NEW_CLASS ... do not change this line
#endif










