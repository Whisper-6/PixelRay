#ifndef THREADPOOL_H
#define THREADPOOL_H

#include <pthread.h>
#include <atomic>

// 公共锁
pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;

template<int THREADS>
class ThreadPool {
private:
    // 每次任务提供给子线程的参数
    struct ThreadArg {
        bool flag;                      // 是否有任务
        void (__cdecl* func)(void*);   // 子线程需要执行的函数
        void* data;                     // 提供给该函数的参数
        pthread_cond_t* timeToWork,* finished;
        std::atomic<int>* cntTask;
    } arg[THREADS];
    pthread_t thread[THREADS];
    pthread_cond_t timeToWork, finished;
    std::atomic<int> cntTask; // 计数还在进行的任务数目
    static void* worker(void* _arg);
public:
    ThreadPool();
    void init();
    void setTask(int id, void (__cdecl* func)(void*), void* data);
    void work();
};

template<int THREADS>
ThreadPool<THREADS>::ThreadPool() {
    timeToWork = finished = PTHREAD_COND_INITIALIZER;
    for (int i = 0; i < THREADS; i++) {
        arg[i].flag = 0;
        arg[i].timeToWork = &timeToWork;
        arg[i].finished = &finished;
        arg[i].cntTask = &cntTask;
    }
}

// 子线程函数
template<int THREADS>
void* ThreadPool<THREADS>::worker(void* _arg) {
    auto arg = (ThreadArg*)_arg;
    while (true) {
        pthread_mutex_lock(&mutex);
        while (!arg->flag || arg->cntTask->load() == 0)
            pthread_cond_wait(arg->timeToWork, &mutex);
        pthread_mutex_unlock(&mutex);
        arg->flag = 0;
        arg->func(arg->data);
        arg->cntTask->fetch_add(-1);
        if (arg->cntTask->load() == 0)
            pthread_cond_signal(arg->finished);
    }
    return nullptr;
}

// 创建所有子线程
template<int THREADS>
void ThreadPool<THREADS>::init() {
    for (int i = 0; i < THREADS; i++) {
        int status = pthread_create(&thread[i], nullptr,
            ThreadPool<THREADS>::worker, (void*)&arg[i]);
        if (status == -1) {
            puts("Failed to create thread.");
            exit(-1);
        }
    }
}

// 给第 id 个子线程布置任务
template<int THREADS>
void ThreadPool<THREADS>::setTask(int id, void (__cdecl* func)(void*), void* data) {
    arg[id].flag = 1;
    arg[id].func = func;
    arg[id].data = data;
}

// 任务布置完毕，让所有子线程开始工作，并等待它们全部完成
template<int THREADS>
void ThreadPool<THREADS>::work() {
    int _cntTask = 0;
    for (int i = 0; i < THREADS; i++)
        _cntTask += arg[i].flag;
    cntTask.store(_cntTask);
    pthread_cond_broadcast(&timeToWork);
    pthread_mutex_lock(&mutex);
    while (cntTask.load() > 0)
        pthread_cond_wait(&finished, &mutex);
    pthread_mutex_unlock(&mutex);
}

#endif