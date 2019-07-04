#pragma once

#include <deque>
#include <functional>
#include <vector>

#include "pthread.h"

namespace ana
{
  class Progress;

  /// \brief A very simple thread pool for use by \ref Surface
  ///
  /// With great power comes great responsibility. Use caution on the grid or
  /// on shared interactive machines before you get yelled at or banned.
  class ThreadPool
  {
  public:
    /// \param maxThreads Maximum number of threads to use at one time.
    ///                   If unspecified, uses number of cores in machine.
    ///                   Use great caution on the grid or on shared
    ///                   interactive machines.
    explicit ThreadPool(unsigned int maxThreads = 0);
    virtual ~ThreadPool();

    void ShowProgress(const std::string& title);

    /// Wait for all threads to complete before returning
    void Finish();

    /// Add task with arguments
    template<class F, class... A> void AddTask(F func, A... args);

    /// Add member function task, with arguments
    template<class T, class M, class... A> void AddMemberTask(T* obj, M meth, A... args);


  protected:
    /// \brief The type of the user's worker functions
    ///
    /// Use std::bind etc to pass arguments
    typedef std::function<void(void)> func_t;

    void AddTaskHelper(func_t func);


    static void* WorkerFunc(void* arg);

    unsigned int fMaxThreads;

    struct Task{func_t func; void* ctx;};

    pthread_mutex_t fTasksLock;
    std::deque<func_t> fTasks;

    ///< Actually, this is protecting \ref fNumLiveThreads
    pthread_mutex_t fThreadsLock;
    std::vector<pthread_t> fThreads; ///< All threads we ever created
    unsigned int fNumLiveThreads; ///< Number of threads that are running

    /// Protects \ref fTasksCompleted and \ref fTotalTasks too
    pthread_mutex_t fProgressLock;
    int fTasksCompleted;
    int fTotalTasks; ///< How many tasks have we ever seen?
    Progress* fProgress;
  };
}


/// Add task with arguments
template<class F, class... A> inline void ana::ThreadPool::AddTask(F func, A... args)
{
  AddTaskHelper(std::bind(func, args...));
}

/// Add member function task, with arguments
template<class T, class M, class... A> inline void ana::ThreadPool::AddMemberTask(T* obj, M meth, A... args)
{
  AddTaskHelper(std::bind(meth, obj, args...));
}
