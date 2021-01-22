//
// cint script to generate libraries
// Declaire namespace & classes you defined
// #pragma statement: order matters! Google it ;)
//

#ifdef __CINT__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ namespace flashmatch+;
#pragma link C++ class flashmatch::PSet+;
/*
#pragma link C++ function flashmatch::PSet::get< string > (const string&)+;
#pragma link C++ function flashmatch::PSet::get< double > (const string&)+;
#pragma link C++ function flashmatch::PSet::get< float > (const string&)+;
#pragma link C++ function flashmatch::PSet::get< int > (const string&)+;
#pragma link C++ function flashmatch::PSet::get< short > (const string&)+;
#pragma link C++ function flashmatch::PSet::get< unsigned int > (const string&)+;
#pragma link C++ function flashmatch::PSet::get< unsigned short > (const string&)+;
#pragma link C++ function flashmatch::PSet::get< size_t > (const string&)+;
*/
#pragma link C++ function flashmatch::ConfigFile2String(const string)+;
#pragma link C++ function flashmatch::CreatePSetFromFile(const string)+;
#pragma link C++ class flashmatch::ConfigManager+;
#pragma link C++ namespace phot+;
#pragma link C++ namespace sim+;
#pragma link C++ class sim::PhotonVoxel+;
#pragma link C++ class sim::PhotonVoxelDef+;
#pragma link C++ class phot::PhotonVisibilityService+;
#pragma link C++ class phot::PhotonLibrary+;

//ADD_NEW_CLASS ... do not change this line
#endif
