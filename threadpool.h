#ifndef THREADPOOL_H
#define THREADPOOL_H

#include <pthread.h>
#include <atomic>

// ������
pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;

template<int THREADS>
class ThreadPool {
private:
    // ÿ�������ṩ�����̵߳Ĳ���
    struct ThreadArg {
        bool flag;                      // �Ƿ�������
        void (__cdecl* func)(void*);   // ���߳���Ҫִ�еĺ���
        void* data;                     // �ṩ���ú����Ĳ���
        pthread_cond_t* timeToWork,* finished;
        std::atomic<int>* cntTask;
    } arg[THREADS];
    pthread_t thread[THREADS];
    pthread_cond_t timeToWork, finished;
    std::atomic<int> cntTask; // �������ڽ��е�������Ŀ
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

// ���̺߳���
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

// �����������߳�
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

// ���� id �����̲߳�������
template<int THREADS>
void ThreadPool<THREADS>::setTask(int id, void (__cdecl* func)(void*), void* data) {
    arg[id].flag = 1;
    arg[id].func = func;
    arg[id].data = data;
}

// ��������ϣ����������߳̿�ʼ���������ȴ�����ȫ�����
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