#include "CAFAna/Core/ThreadPool.h"

#include "CAFAna/Core/Progress.h"

#include <cassert>
#include <unistd.h>

namespace ana
{
  //----------------------------------------------------------------------
  ThreadPool::ThreadPool(unsigned int maxThreads)
    : fMaxThreads(maxThreads), fNumLiveThreads(0),
      fTasksCompleted(0), fTotalTasks(0), fProgress(0)
  {
    if(maxThreads == 0){
      fMaxThreads = sysconf(_SC_NPROCESSORS_ONLN);
    }

    pthread_mutex_init(&fTasksLock, 0);
    pthread_mutex_init(&fThreadsLock, 0);
    pthread_mutex_init(&fProgressLock, 0);
  }

  //----------------------------------------------------------------------
  ThreadPool::~ThreadPool()
  {
    // This is safe, even if it's already been explicitly called
    Finish();

    pthread_mutex_destroy(&fTasksLock);
    pthread_mutex_destroy(&fThreadsLock);
    pthread_mutex_destroy(&fProgressLock);

    delete fProgress;
  }

  //----------------------------------------------------------------------
  void ThreadPool::ShowProgress(const std::string& title)
  {
    if(!fProgress) fProgress = new Progress(title);
  }

  //----------------------------------------------------------------------
  void ThreadPool::Finish()
  {
    void* junk;
    for(pthread_t& th: fThreads) pthread_join(th, &junk);

    assert(fNumLiveThreads == 0);

    fThreads.clear(); // Make it safe to call a second time

    if(fProgress) fProgress->Done();
  }

  //----------------------------------------------------------------------
  void* ThreadPool::WorkerFunc(void* arg)
  {
    // We smuggle essentially the this pointer in through the argument
    ThreadPool* pool = (ThreadPool*)arg;

    while(true){
      pthread_mutex_lock(&pool->fTasksLock);
      if(pool->fTasks.empty()){
        // Nothing to do, commit suicide
        pthread_mutex_unlock(&pool->fTasksLock);

        pthread_mutex_lock(&pool->fThreadsLock);
        --pool->fNumLiveThreads;
        pthread_mutex_unlock(&pool->fThreadsLock);

        return 0;
      }
      // Get a task to do
      func_t task = pool->fTasks.front();
      pool->fTasks.pop_front();
      pthread_mutex_unlock(&pool->fTasksLock);

      // Actually do the user's work
      task();

      pthread_mutex_lock(&pool->fProgressLock);
      ++pool->fTasksCompleted;

      if(pool->fProgress)
        pool->fProgress->SetProgress(pool->fTasksCompleted/double(pool->fTotalTasks));
      pthread_mutex_unlock(&pool->fProgressLock);
    }
  }

  //----------------------------------------------------------------------
  void ThreadPool::AddTaskHelper(func_t func)
  {
    pthread_mutex_lock(&fTasksLock);
    fTasks.push_back(func);
    pthread_mutex_unlock(&fTasksLock);

    pthread_mutex_lock(&fProgressLock);
    ++fTotalTasks;
    pthread_mutex_unlock(&fProgressLock);

    pthread_mutex_lock(&fThreadsLock);
    if(fNumLiveThreads < fMaxThreads){
      fThreads.push_back(pthread_t());
      pthread_create(&fThreads.back(), 0, WorkerFunc, this);
      ++fNumLiveThreads;
    }
    pthread_mutex_unlock(&fThreadsLock);
  }
}
