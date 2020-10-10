#ifndef BASICTOOL_GEOOBJCOLLECTION_CXX
#define BASICTOOL_GEOOBJCOLLECTION_CXX

#include "GeoObjCollection.h"

namespace geoalgo {

  GeoObjCollection::GeoObjCollection() 
  {}

  void GeoObjCollection::Clear()
  {
    _pt_v.clear();
    _box_v.clear();
    _seg_v.clear();
    _trj_v.clear();
    _lin_v.clear();
    _cone_v.clear();
    //_sphere_v.clear();
    _pt_col.clear();
    _box_col.clear();
    _seg_col.clear();
    _trj_col.clear();
    _lin_col.clear();
    _cone_col.clear();
    _sphere_col.clear();
    _labels.clear();
  }

  void GeoObjCollection::_AddLabel_(const Point_t& pt, std::string label)
  {
    auto const iter = _labels.find(pt);
    if(iter==_labels.end()) 

      _labels.insert(std::make_pair(pt,Form("%s\n",label.c_str())));

    else {

      label = Form("%s%s\n",(*iter).second.c_str(),label.c_str());

      _labels[pt] = label;
    }
  }

  void GeoObjCollection::Add(const Point_t& pt, std::string name, std::string c)
  {
    if(name.empty()) name = Form("Pt (%zu)",_pt_v.size());
    _AddLabel_(pt,name);
    _pt_v.push_back(pt);
    _pt_col.push_back(c);
  }

  void GeoObjCollection::Add(const AABox_t& box, std::string name, std::string c)
  {
    if(name.empty()) name = Form("AABox (%zu)",_box_v.size());
    _AddLabel_(box.Min() + (box.Max() - box.Min())/2.,name);
    _box_v.push_back(box);
    _box_col.push_back(c);
  }

  void GeoObjCollection::Add(const LineSegment_t& seg, std::string name, std::string c)
  {
    if(name.empty()) name = Form("LSeg (%zu)",_seg_v.size());
    _AddLabel_(seg.End(),name);
    _seg_v.push_back(seg);
    _seg_col.push_back(c);
  }


  void GeoObjCollection::Add(const HalfLine_t& lin, std::string name, std::string c)
  {
    if(name.empty()) name = Form("Line (%zu)",_lin_v.size());
    _AddLabel_(lin.Start()+lin.Start()*10,name);
    _lin_v.push_back(lin);
    _lin_col.push_back(c);
  }

  void GeoObjCollection::Add(const Trajectory_t& trj, std::string name, std::string c)
  {
    if(trj.size()<2)
      throw GeoAlgoException("Trajectory size cannot be smaller than 2!");
    if(name.empty()) name = Form("Trj (%zu)",_trj_v.size());
    _AddLabel_(trj.back(),name);
    _trj_v.push_back(trj);
    _trj_col.push_back(c);

    return;
  }

  void GeoObjCollection::Add(const Cone_t& cone, std::string name, std::string c)
  {
    if(name.empty()) name = Form("Cone (%zu)",_cone_v.size());
    _AddLabel_(cone.Start()+cone.Dir()*cone.Length(),name);
    _cone_v.push_back(cone);
    _cone_col.push_back(c);
  }


  void GeoObjCollection::Add(const Sphere_t& sphere, std::string name, std::string c)
  {
    if(name.empty()) name = Form("Sphere (%zu)",_sphere_v.size());
    _AddLabel_(sphere.Center(),name);
    _sphere_v.push_back(sphere);
    _sphere_col.push_back(c);
  }

}

#endif
