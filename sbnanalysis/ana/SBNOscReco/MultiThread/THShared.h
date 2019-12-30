#ifndef _THShared_hh_
#define _THShared_hh_

#include <memory>
#include "HistoManager.h"

class TH1;
class TH2;
class TH3;

class TH1Shared {
public:
  TH1Shared(): hist(NULL) {}
  TH1Shared(TH1 *_hist);
  void Fill(double x, double scale=1.);
  TH1 *Get();
  const TH1 *Get() const;
private:
  std::shared_ptr<TH1> hist;
};

class TH2Shared {
public:
  TH2Shared(): hist(NULL) {}
  TH2Shared(TH2 *_hist);
  void Fill(double x, double y, double scale=1.);
  TH2 *Get();
  const TH2 *Get() const;
private:
  std::shared_ptr<TH2> hist;
};

class TH3Shared {
public:
  TH3Shared(): hist(NULL) {}
  TH3Shared(TH3 *_hist);
  void Fill(double x, double y, double z,double scale=1.);
  TH3 *Get();
  const TH3 *Get() const;
private:
  std::shared_ptr<TH3> hist;
};
#endif
