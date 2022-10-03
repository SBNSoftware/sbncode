#ifndef __OPT0FINDER_PYUTILS_CXX__
#define __OPT0FINDER_PYUTILS_CXX__

#include "PyUtils.h"
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
//#include <numpy/ndarrayobject.h>
#include "numpy/arrayobject.h"
#include "flashmatch/Base/FMWKInterface.h"
#include <cassert>
#include <iostream>

namespace flashmatch {

  void SetPyUtil() {
    static bool once = false;
    if (!once) {
      _import_array();
      once = true;
    }
  }

  PyObject* as_ndarray(const QCluster_t& traj) {
    SetPyUtil();
    std::vector<size_t> dims(2,0);
    dims[0] = traj.size();
    dims[1] = 5;
    auto pyarray = numpy_array<double>(dims);

    double **carray;
    const int dtype = NPY_DOUBLE;
    PyArray_Descr *descr = PyArray_DescrFromType(dtype);
    npy_intp pydims[2];
    if (PyArray_AsCArray(&pyarray, (void **)&carray, pydims, 2, descr) < 0) {
      std::cerr<<"Failed to create 2D numpy array"<<std::endl;
      throw std::exception();
    }
    assert(pydims[0] == ((int)(dims[0])) && pydims[1] == ((int)(dims[1])));

    for(size_t idx=0; idx<dims[0]; ++idx) {
      carray[idx][0] = traj[idx].x;
      carray[idx][1] = traj[idx].y;
      carray[idx][2] = traj[idx].z; 
      carray[idx][3] = traj[idx].q;
      carray[idx][4] = traj[idx].trk;
    }
    PyArray_Free(pyarray,  (void *)carray);
    return pyarray;
  }

  PyObject* as_ndarray(const ::geoalgo::Trajectory& traj) {
    SetPyUtil();
    std::vector<size_t> dims(2,0);
    dims[0] = traj.size();
    dims[1] = 3;
    auto pyarray = numpy_array<double>(dims);

    double **carray;
    const int dtype = NPY_DOUBLE;
    PyArray_Descr *descr = PyArray_DescrFromType(dtype);
    npy_intp pydims[2];
    if (PyArray_AsCArray(&pyarray, (void **)&carray, pydims, 2, descr) < 0) {
      std::cerr<<"Failed to create 2D numpy array"<<std::endl;
      throw std::exception();
    }
    assert(pydims[0] == ((int)(dims[0])) && pydims[1] == ((int)(dims[1])));

    for(size_t idx=0; idx<dims[0]; ++idx) {
      carray[idx][0] = traj[idx][0];
      carray[idx][1] = traj[idx][1];
      carray[idx][2] = traj[idx][2];
    }
    PyArray_Free(pyarray,  (void *)carray);
    return pyarray;
  }

  PyObject* as_ndarray(const Flash_t& flash) {
    SetPyUtil();
    std::vector<size_t> dims(1);
    dims[0] = flash.pe_v.size();
    auto pyarray = numpy_array<double>(dims);

    double *carray;
    const int dtype = NPY_DOUBLE;
    PyArray_Descr *descr = PyArray_DescrFromType(dtype);
    npy_intp pydims[1];
    if (PyArray_AsCArray(&pyarray, (void **)&carray, pydims, 1, descr) < 0) {
      std::cerr<<"Failed to create 2D numpy array"<<std::endl;
      throw std::exception();
    }
    assert(pydims[0] == ((int)(dims[0])));

    for(size_t idx=0; idx<dims[0]; ++idx)
      carray[idx] = flash.pe_v[idx];
    PyArray_Free(pyarray,  (void *)carray);
    return pyarray;
  }

  /*
    void copy_array(PyObject *arrayin, const std::vector<float> &cvec) {
    SetPyUtil();
    PyArrayObject *ptr = (PyArrayObject *)(arrayin);
    
    //std::cout<< PyArray_NDIM(ptr) << std::endl
    //         << PyArray_DIM(ptr,0)<<std::endl
    //         << PyArray_SIZE(ptr) << std::endl;
    
    // Check dimension size is 1:
    if (PyArray_NDIM(ptr) != 1){
    throw std::exception();
    }
    
    if ((long)(cvec.size()) != PyArray_SIZE(ptr))
    throw std::exception();
    npy_intp loc[1];
    loc[0] = 0;
    auto fptr = (float *)(PyArray_GetPtr(ptr, loc));
    for (size_t i = 0; i < size_t(PyArray_SIZE(ptr)); ++i) {
    // std::cout << fptr[i] << std::endl;
    fptr[i] = cvec[i];
    };
    }
  */
  
  template<class T>
  void _copy_array(PyObject *arrayin, const std::vector<T> &cvec) {
    SetPyUtil();
    PyArrayObject *ptr = (PyArrayObject *)(arrayin);
    
    //std::cout<< PyArray_NDIM(ptr) << std::endl
    //         << PyArray_DIM(ptr,0)<<std::endl
    //         << PyArray_SIZE(ptr) << std::endl;
    
    // Check dimension size is 1:
    if (PyArray_NDIM(ptr) != 1){
      throw std::exception();
    }
    
    if ((long)(cvec.size()) != PyArray_SIZE(ptr))
      throw std::exception();
    npy_intp loc[1];
    loc[0] = 0;
    auto fptr = (T *)(PyArray_GetPtr(ptr, loc));
    for (size_t i = 0; i < size_t(PyArray_SIZE(ptr)); ++i) {
      // std::cout << fptr[i] << std::endl;
      fptr[i] = cvec[i];
    };
  }
  
  template void _copy_array< unsigned short >(PyObject *arrayin, const std::vector< unsigned short > &cvec);
  template void _copy_array< unsigned int   >(PyObject *arrayin, const std::vector< unsigned int   > &cvec);
  template void _copy_array< short          >(PyObject *arrayin, const std::vector< short          > &cvec);
  template void _copy_array< int            >(PyObject *arrayin, const std::vector< int            > &cvec);
  template void _copy_array< long long      >(PyObject *arrayin, const std::vector< long long      > &cvec);
  template void _copy_array< float          >(PyObject *arrayin, const std::vector< float          > &cvec);
  template void _copy_array< double         >(PyObject *arrayin, const std::vector< double         > &cvec);
  
  void copy_array(PyObject *arrayin, const std::vector< unsigned short > &cvec) { _copy_array(arrayin, cvec); }
  void copy_array(PyObject *arrayin, const std::vector< unsigned int   > &cvec) { _copy_array(arrayin, cvec); }
  void copy_array(PyObject *arrayin, const std::vector< short          > &cvec) { _copy_array(arrayin, cvec); }
  void copy_array(PyObject *arrayin, const std::vector< int            > &cvec) { _copy_array(arrayin, cvec); }
  void copy_array(PyObject *arrayin, const std::vector< long long      > &cvec) { _copy_array(arrayin, cvec); }
  void copy_array(PyObject *arrayin, const std::vector< float          > &cvec) { _copy_array(arrayin, cvec); }
  void copy_array(PyObject *arrayin, const std::vector< double         > &cvec) { _copy_array(arrayin, cvec); }
  
  template<> int ctype_to_numpy<short>() { SetPyUtil(); SetPyUtil(); return NPY_INT16; }
  template<> int ctype_to_numpy<unsigned short>() { SetPyUtil(); return NPY_UINT16; }
  template<> int ctype_to_numpy<int>() { SetPyUtil(); return NPY_INT32; }
  template<> int ctype_to_numpy<unsigned int>() { SetPyUtil(); return NPY_UINT32; }
  template<> int ctype_to_numpy<long long>() { SetPyUtil(); return NPY_INT64; }
  template<> int ctype_to_numpy<unsigned long long>() { SetPyUtil(); return NPY_UINT64; }
  template<> int ctype_to_numpy<float>() { SetPyUtil(); return NPY_FLOAT32; }
  template<> int ctype_to_numpy<double>() { SetPyUtil(); return NPY_FLOAT64; }
  
  /*
    PyObject *as_ndarray(const std::vector<float> &vec) {
    SetPyUtil();
    
    if (vec.size() >= INT_MAX) {
    LARCV_CRITICAL() << "Length of data vector too long to specify ndarray. "
    "Use by batch call."
    << std::endl;
    throw larbys();
    }
    int nd = 1;
    npy_intp dims[1];
    dims[0] = (int)vec.size();
    PyArrayObject *array = (PyArrayObject *)PyArray_SimpleNewFromData(
    nd, dims, NPY_FLOAT, (char *)&(vec[0]));
    return PyArray_Return(array);
    }
  */
  
  template <class T>
  PyObject *_as_ndarray(const std::vector<T> &vec) {
    SetPyUtil();
    
    if (vec.size() >= INT_MAX) {
      std::cerr << "Length of data vector too long to specify ndarray. " << std::endl;
      throw std::exception();
    }
    int nd = 1;
    npy_intp dims[1];
    dims[0] = (int)vec.size();
    PyArrayObject *array = (PyArrayObject *)PyArray_SimpleNewFromData(
								      nd, dims, ctype_to_numpy<T>(), (char *)&(vec[0]));
    return PyArray_Return(array);
  }
  
  template PyObject* _as_ndarray< short              > (const std::vector< short              >& vec);
  template PyObject* _as_ndarray< unsigned short     > (const std::vector< unsigned short     >& vec);
  template PyObject* _as_ndarray< int                > (const std::vector< int                >& vec);
  template PyObject* _as_ndarray< unsigned int       > (const std::vector< unsigned int       >& vec);
  template PyObject* _as_ndarray< long long          > (const std::vector< long long          >& vec);
  template PyObject* _as_ndarray< unsigned long long > (const std::vector< unsigned long long >& vec);
  template PyObject* _as_ndarray< float              > (const std::vector< float              >& vec);
  template PyObject* _as_ndarray< double             > (const std::vector< double             >& vec);
  
  PyObject* as_ndarray(const std::vector< short              >& vec) { return _as_ndarray< short              >(vec); }
  PyObject* as_ndarray(const std::vector< unsigned short     >& vec) { return _as_ndarray< unsigned short     >(vec); }
  PyObject* as_ndarray(const std::vector< int                >& vec) { return _as_ndarray< int                >(vec); }
  PyObject* as_ndarray(const std::vector< unsigned int       >& vec) { return _as_ndarray< unsigned int       >(vec); }
  PyObject* as_ndarray(const std::vector< long long          >& vec) { return _as_ndarray< long long          >(vec); }
  PyObject* as_ndarray(const std::vector< unsigned long long >& vec) { return _as_ndarray< unsigned long long >(vec); }
  PyObject* as_ndarray(const std::vector< float              >& vec) { return _as_ndarray< float              >(vec); }
  PyObject* as_ndarray(const std::vector< double             >& vec) { return _as_ndarray< double             >(vec); }
  
  template<class T>
  PyObject* numpy_array(std::vector<size_t> dims)
  {
    SetPyUtil();
    int nd_ = dims.size();
    npy_intp dims_[nd_];
    for(size_t i=0; i<dims.size(); ++i) dims_[i] = dims[i];
    /*
      std::cout<<"NUMPY TYPE " << ctype_to_numpy<T>() << std::endl;
      std::cout<<"DIMS " << dims.size() << std::endl;
      std::cout<<"ND " << nd_ << std::endl;
      std::cout<<"Shape ";
      for(size_t i=0; i<dims.size(); ++i) std::cout<< " " << dims_[i];
      std::cout<<std::endl;
    */
    //PyObject* res = PyArray_SimpleNew(nd_,dims_,ctype_to_numpy<T>());
    PyObject* res = PyArray_ZEROS(nd_,dims_,ctype_to_numpy<T>(),0);
    //Py_INCREF(res);
    /*
      std::cout<<PyArray_NDIM((PyArrayObject*)res) << std::endl;
      std::cout<<PyArray_SIZE((PyArrayObject*)res) << std::endl;
    */
    PyArrayObject *ptr = (PyArrayObject*)(res);
    // Check dimension size is 1:
    //std::cout<<"ndim " << PyArray_NDIM(ptr) << std::endl;
    //size_t len = PyArray_SIZE(ptr);
    //std::cout<<"len " << len <<std::endl;
    //npy_intp loc[1];
    //loc[0] = 0;
    //auto fptr = (T *)(PyArray_GetPtr(ptr, loc));
    /*
      std::cout<<"fptr " << fptr << std::endl;
      for (size_t i = 0; i < len; ++i) {
      std::cout << fptr[i] << std::endl;
      fptr[i] = T(1);
      }
    */
    PyArray_INCREF(ptr);
    
    return res;
  }
  
  template PyObject* numpy_array<float>(std::vector<size_t>dims);
  
}
#endif
  
