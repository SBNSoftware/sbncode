#include "CAFAna/Prediction/IncrementalCholeskyDecomp.h"

#include "CAFAna/Core/MathUtil.h"

#include <cmath>
#include <iostream>

namespace ana
{
  // --------------------------------------------------------------------------
  void IncrementalCholeskyDecomp::Extend()
  {
    // First up, expand the decomp matrix
    for(std::vector<double>& row: fUpper) row.push_back(0); // extend each row
    fUpper.emplace_back(fUpper.size()+1, 0); // and add one more
  }

  // --------------------------------------------------------------------------
  void IncrementalCholeskyDecomp::SetLastRow(const std::vector<double>& row)
  {
    // Implementation based on https://en.wikipedia.org/wiki/Cholesky_decomposition#Adding_and_removing_rows_and_columns

    const unsigned int N = fUpper.size();

    // S11 = L11  -- true automatically

    // S12 = L11^T \ A12
    for(unsigned int i = 0; i+1 < N; ++i){ // which row we're solving
      double res = row[i];
      for(unsigned int j = 0; j < i; ++j) res -= fUpper[j][i] * fUpper[j][N-1];
      if(res != 0) fUpper[i][N-1] = res/fUpper[i][i];
    }

    // S22 = chol(A22 - S12T S12)
    double dot = 0;
    for(unsigned int i = 0; i+1 < N; ++i) dot += util::sqr(fUpper[i][N-1]);
    fUpper[N-1][N-1] = sqrt(row[N-1] - dot);
  }

  // --------------------------------------------------------------------------
  void IncrementalCholeskyDecomp::SetLastRow(const std::vector<double>& bigrow,
                                             const std::vector<int>& idxs)
  {
    std::vector<double> row;
    row.reserve(idxs.size());
    for(int i: idxs) row.push_back(bigrow[i]);
    SetLastRow(row);
  }

  // --------------------------------------------------------------------------
  std::vector<double> IncrementalCholeskyDecomp::Solve(const std::vector<double>& b) const
  {
    const unsigned int N = fUpper.size();

    // Solve Ly = b
    std::vector<double> y(N);

    for(unsigned int i = 0; i < N; ++i){ // which row we're solving
      double res = b[i];
      for(unsigned int j = 0; j < i; ++j) res -= fUpper[j][i] * y[j];
      if(res != 0) y[i] = res/fUpper[i][i];
    }

    // And now LTx = y
    std::vector<double> x(N);

    for(int i = N-1; i >= 0; --i){
      double res = y[i];
      for(unsigned int j = i+1; j < N; ++j) res -= fUpper[i][j] * x[j];
      if(res != 0) x[i] = res/fUpper[i][i];
    }

    return x;
  }

  // --------------------------------------------------------------------------
  void IncrementalCholeskyDecomp::Print() const
  {
    for(const std::vector<double>& row: fUpper){
      for(double v: row) std::cout << v << "\t";
      std::cout << std::endl;
    }
  }
}
