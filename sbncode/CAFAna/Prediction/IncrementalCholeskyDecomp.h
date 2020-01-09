#pragma once

#include <vector>

namespace ana
{
  // Based on https://en.wikipedia.org/wiki/Cholesky_decomposition#Adding_and_removing_rows_and_columns
  class IncrementalCholeskyDecomp
  {
  public:
    /// Initialize for a 0x0 matrix
    IncrementalCholeskyDecomp() {}

    /// Add one row and column of zeros
    void Extend();

    /// This is of course also the right-most column...
    void SetLastRow(const std::vector<double>& row);

    void AddRow(const std::vector<double>& row)
    {
      Extend();
      SetLastRow(row);
    }

    std::vector<double> Solve(const std::vector<double>& b) const;

    void Print() const;

  protected:
    std::vector<std::vector<double>> fUpper;
  };
}
