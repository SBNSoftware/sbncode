//
// cint script to generate libraries
// Declaire namespace & classes you defined
// #pragma statement: order matters! Google it ;)
//

#ifdef __CINT__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ namespace geoalgo+;
#pragma link C++ class geoalgo::Vector+;
#pragma link C++ class std::vector<geoalgo::Vector>+;
#pragma link C++ class std::vector<std::vector<geoalgo::Vector> >+;
#pragma link C++ class std::map<geoalgo::Vector,string>+;
#pragma link C++ class geoalgo::Trajectory+;
#pragma link C++ class std::vector<geoalgo::Trajectory>+;
#pragma link C++ class geoalgo::HalfLine+;
#pragma link C++ class std::vector<geoalgo::HalfLine>+;
#pragma link C++ class geoalgo::Line+;
#pragma link C++ class std::vector<geoalgo::Line>+;
#pragma link C++ class geoalgo::DirectedLine+;
#pragma link C++ class std::vector<geoalgo::DirectedLine>+;
#pragma link C++ class geoalgo::LineSegment+;
#pragma link C++ class std::vector<geoalgo::LineSegment>+;
#pragma link C++ class geoalgo::AABox+;
#pragma link C++ class std::vector<geoalgo::AABox>+;
#pragma link C++ class geoalgo::Cylinder+;
#pragma link C++ class std::vector<geoalgo::Cylinder>+;
#pragma link C++ class geoalgo::Cone+;
#pragma link C++ class std::vector<geoalgo::Cone>+;
#pragma link C++ class geoalgo::Sphere+;
#pragma link C++ class std::vector<geoalgo::Sphere>+;
#pragma link C++ class std::pair<geoalgo::Vector,string>+;
#pragma link C++ class std::map<geoalgo::Vector,string>+;

#pragma link C++ class geoalgo::GeoAlgo+;
#pragma link C++ class geoalgo::GeoObjCollection+;
//ADD_NEW_CLASS ... do not change this line

#endif
















