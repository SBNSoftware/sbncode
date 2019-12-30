#ifndef _HISTOManager_hh_
#define _HISTOManager_hh_

#include <vector>
#include <thread>

#include "concurrentqueue/concurrentqueue.h"

#include "TH1.h"
#include "TH2.h"
#include "TH3.h"

class HistoManager {
public:
  ~HistoManager();

  static void Fill(TH1 *hist, double x, double scale=1.);
  static void Fill(TH2 *hist, double x, double y, double scale=1.);
  static void Fill(TH3 *hist, double x, double y, double z, double scale=1.);

  static unsigned SizeHint();
  static void Start(unsigned _n_instances=0);
  static void Finish();

private:
  explicit HistoManager(unsigned n_instance=1, bool start=false);

  static HistoManager &Instance();

  void DoStart(unsigned _n_instances=0);
  void DoFinish();

  enum FillType {
   k1D, k2D, k3D
  };

  struct FillData {
    FillType type; 
    TH1 *ptr;
    double x;
    double y;
    double z;
    double scale;
  };

  template <typename T>
  unsigned Hash(T *ptr) {
    return ((size_t)ptr / std::alignment_of_v<T>) % n_instances;
  }

  unsigned n_instances;
  std::vector<moodycamel::ConcurrentQueue<FillData>> queues;
  bool finished;
  bool cleaned;
  std::vector<std::thread> workers;
  friend void WorkDeque(moodycamel::ConcurrentQueue<FillData> &queue, bool &finished);
};

void WorkDeque(moodycamel::ConcurrentQueue<HistoManager::FillData> &queue, bool &finished);


#endif
