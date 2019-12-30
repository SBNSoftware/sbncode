#include "HistoManager.h"

#include "TH1.h"
#include "TH2.h"
#include "TH3.h"

#include <unistd.h>

HistoManager::HistoManager(unsigned _n_instances, bool start):
  n_instances(_n_instances),
  queues(_n_instances)
{
  finished = false;
  cleaned = false;
  if (start) DoStart();
}

HistoManager& HistoManager::Instance() {
  static HistoManager manager(1, true);
  return manager;
}

void HistoManager::Start(unsigned _n_instances) {
  HistoManager &manager = HistoManager::Instance();
  manager.DoStart(_n_instances);
  return;

}
void HistoManager::DoStart(unsigned _n_instances) {
  if (_n_instances != 0 && _n_instances != n_instances) {
    n_instances = _n_instances;
    queues = std::vector<moodycamel::ConcurrentQueue<FillData>>(n_instances);
  }
  for (unsigned i = 0; i < n_instances; i++) {
    workers.emplace_back(WorkDeque, std::ref(queues[i]), std::ref(finished));
  }
  return;
}

void HistoManager::Finish() {
  HistoManager &manager = HistoManager::Instance();
  manager.DoFinish();
  return;
}

void HistoManager::DoFinish() {
  if (cleaned) return;
  finished = true; 
  for (unsigned i = 0; i < n_instances; i++) {
    workers[i].join();
  }
  cleaned = true;
  return;
}

unsigned HistoManager::SizeHint() {
  HistoManager &manager = HistoManager::Instance();
  unsigned ret = 0;
  for (const moodycamel::ConcurrentQueue<FillData> &q: manager.queues) {
    ret += q.size_approx();
  }
  return ret;
}


HistoManager::~HistoManager() {
  DoFinish();
}

void HistoManager::Fill(TH1 *hist, double x, double scale) {
  HistoManager &manager = HistoManager::Instance();

  HistoManager::FillData f;
  f.type = HistoManager::k1D;
  f.ptr = hist;
  f.x = x;
  f.scale = scale;
  unsigned index = manager.Hash(hist);
  manager.queues[index].enqueue(f);
  return;
}

void HistoManager::Fill(TH2 *hist, double x, double y, double scale) {
  HistoManager &manager = HistoManager::Instance();

  HistoManager::FillData f;
  f.type = HistoManager::k2D;
  f.ptr = hist;
  f.x = x;
  f.y = y;
  f.scale = scale;
  unsigned index = manager.Hash(hist);
  manager.queues[index].enqueue(f);
  return;
}

void HistoManager::Fill(TH3 *hist, double x, double y, double z, double scale) {
  HistoManager &manager = HistoManager::Instance();

  HistoManager::FillData f;
  f.type = HistoManager::k3D;
  f.ptr = hist;
  f.x = x;
  f.y = y;
  f.z = z;
  f.scale = scale;
  unsigned index = manager.Hash(hist);
  manager.queues[index].enqueue(f);
  return;
}

void WorkDeque(moodycamel::ConcurrentQueue<HistoManager::FillData> &queue, bool &finished) {
  HistoManager::FillData f;

  while (true) {
    bool success = queue.try_dequeue(f);
    if (success) {
      switch (f.type) {
        case HistoManager::k1D:
          f.ptr->Fill(f.x, f.scale);
          break;
        case HistoManager::k2D:
          ((TH2*)f.ptr)->Fill(f.x, f.y, f.scale);
          break;
        case HistoManager::k3D:
          ((TH3*)f.ptr)->Fill(f.x, f.y, f.z, f.scale);
          break;
      }
    }
    else {
      if (finished) break;
      usleep(1000); // if there's nothing there, wait a millisecond
    }
  }

  return;
}


