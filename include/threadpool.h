#ifndef PIXEL_RAY_THREAD_POOL_H
#define PIXEL_RAY_THREAD_POOL_H

#include <pthread.h>
#include <atomic>

pthread_mutex_t _thread_pool_mutex = PTHREAD_MUTEX_INITIALIZER;

template<int THREADS>
class ThreadPool {
private:
    struct ThreadArg {
        bool flag;                      // Whether there is a task
        void (__cdecl* func)(void*);    // Function that the thread needs to execute
        void* data;                     // Parameter provided to the function
        pthread_cond_t* timeToWork,* finished;
        std::atomic<int>* cntTask;
    } arg[THREADS];
    pthread_t thread[THREADS];
    pthread_cond_t timeToWork, finished;
    std::atomic<int> cntTask; // Count the number of tasks still in progress
    static void* worker(void* _arg);
public:
    ThreadPool();
    void init();
    void setTask(int id, void (__cdecl* func)(void*), void* data);
    void work();
};

#include "threadpool.hpp"

#endif //PIXEL_RAY_THREAD_POOL_H
