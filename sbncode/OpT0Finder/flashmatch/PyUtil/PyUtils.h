#ifndef __OPT0FINDER_PYUTILS_H__
#define __OPT0FINDER_PYUTILS_H__

struct _object;
typedef _object PyObject;

#include <Python.h>
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
//#include <numpy/ndarrayobject.h>
#include "numpy/arrayobject.h"
#include <vector>
#include "flashmatch/Base/OpT0FinderTypes.h"
#include "flashmatch/GeoAlgo/GeoTrajectory.h"
namespace flashmatch {
  /// Utility function: call one-time-only numpy module initialization (you don't have to call)
  void SetPyUtil();

  ::geoalgo::Trajectory as_geoalgo_trajectory(PyObject* pyarray);
  PyObject* as_ndarray(const QCluster_t& traj);
  PyObject* as_ndarray(const Flash_t& traj);
  PyObject* as_ndarray(const ::geoalgo::Trajectory& traj);
  
  /// convert vectors into np array
  template <class T>
  PyObject* _as_ndarray(const std::vector<T>& data);
  PyObject* as_ndarray(const std::vector< short              > &data);
  PyObject* as_ndarray(const std::vector< unsigned short     > &data);
  PyObject* as_ndarray(const std::vector< int                > &data);
  PyObject* as_ndarray(const std::vector< unsigned int       > &data);
  PyObject* as_ndarray(const std::vector< long long          > &data);
  PyObject* as_ndarray(const std::vector< unsigned long long > &data);
  PyObject* as_ndarray(const std::vector< float              > &data);
  PyObject* as_ndarray(const std::vector< double             > &data);
  
  /// copy array
  template <class T>
  void _copy_array(PyObject *arrayin, const std::vector<T> &cvec);
  void copy_array(PyObject *arrayin, const std::vector< unsigned short > &cvec);
  void copy_array(PyObject *arrayin, const std::vector< unsigned int   > &cvec);
  void copy_array(PyObject *arrayin, const std::vector< short          > &cvec);
  void copy_array(PyObject *arrayin, const std::vector< int            > &cvec);
  void copy_array(PyObject *arrayin, const std::vector< long long      > &cvec);
  void copy_array(PyObject *arrayin, const std::vector< float          > &cvec);
  void copy_array(PyObject *arrayin, const std::vector< double         > &cvec);
  
  template <class T> int ctype_to_numpy();
  template<> int ctype_to_numpy<short>();
  template<> int ctype_to_numpy<unsigned short>();
  template<> int ctype_to_numpy<int>();
  template<> int ctype_to_numpy<unsigned int>();
  template<> int ctype_to_numpy<long long>();
  template<> int ctype_to_numpy<unsigned long long>();
  template<> int ctype_to_numpy<float>();
  template<> int ctype_to_numpy<double>();
  template <class T> PyObject* numpy_array(std::vector<size_t> dims);
}

#endif
