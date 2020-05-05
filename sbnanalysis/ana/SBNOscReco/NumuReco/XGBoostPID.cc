#include "XGBoostPID.h"
#include <iostream>
#include <cmath>
#include <cassert>
#include <vector>


ana::SBNOsc::XGBoostPID::XGBoostPID():
  fReady(false)
{
  XGBoosterCreate(NULL, 0, &fHandle);
  XGBoosterSetParam(fHandle, "seed", "0");
  // XGBoosterSetParam(fHandle, "seed", "");
}

void ana::SBNOsc::XGBoostPID::SetModelFile(const char *modelFile) {
  std::cout << "Loading model file: " << modelFile << std::endl;
  XGBoosterLoadModel(fHandle, modelFile);
  fReady = true;
}

float ana::SBNOsc::XGBoostPID::PredictOne(const std::map<std::string, float> &data) const {
  DMatrixHandle Matrix;

  std::vector<float> mat {
    data.at("phi"),
    data.at("theta"),
    data.at("length"),
    data.at("crt_hit_distance"),
    data.at("chi2_muon"),
    data.at("chi2_diff"),
    data.at("chi2_proton"),
    data.at("contained")
  };

  std::cout << "Inp: ";
  for (float f: mat)  std::cout << f << " ";
  std::cout << std::endl;

  int fail;
  fail = XGDMatrixCreateFromMat(&mat[0], 1, mat.size(), NAN, &Matrix);
  std::cout << "Matrix failed: " << fail << std::endl;

  const float *f = NULL;
  bst_ulong out_len;
  std::cout << "Score ptr: " << f << std::endl;
  XGBoosterPredict(fHandle, Matrix, 0, 0, 0, &out_len, &f);
  std::cout << "BDT failed: " << fail << std::endl;
  std::cout << "Score ptr: " << f << std::endl;
  std::cout << "Score: " << *f << std::endl;
  assert(out_len == 1);
  XGDMatrixFree(Matrix);

  return *f;
}

ana::SBNOsc::XGBoostPID::~XGBoostPID() {
  XGBoosterFree(fHandle);
}

