#ifndef THREADPOOL
#define THREADPOOL
#include <thread>
#include <iostream>
#include <vector>
#include <queue>
#include <mutex>
#include <condition_variable>
#include <atomic>
#include <functional>
#include <future>
#include "spdlog/spdlog.h"
class ThreadPool{
private:
    std::vector<std::thread> workers;
    std::queue<std::function<void()>> Tasks;
    std::mutex queMutex;
    std::condition_variable condition;
    std::atomic<bool> stop;
public:
    ThreadPool(size_t n);
    

    ~ThreadPool(){ // 获得锁，设置stop变量，通知所有的thread执行自己已经分配的任务
        {
            std::unique_lock<std::mutex> lock(this->queMutex);
            stop = true;
        }
        condition.notify_all();
        for(auto& thread : workers){
            thread.join();
        }
    }

    template<class F , class ... Args>
    auto enqueue(F&& f, Args ...args) -> std::future<typename std::invoke_result<F,Args...>::type>{
       using result_type =typename std::invoke_result<F,Args...>::type;

        auto task = std::make_shared<std::packaged_task<result_type()>>(
            std::bind(std::forward<F>(f), std::forward<Args>(args) ...)
        );
        std::future<result_type> res = task->get_future();
        {
            std::unique_lock<std::mutex> lock(this->queMutex);
            if(stop) throw std::runtime_error("ThreadPool Stop!");
            this->Tasks.emplace([task]{
                (*task)();
            });
            // spdlog::info("this Tasks size() : {}" , this->Tasks.size());
        }
        
        condition.notify_one();
        return res;
    }
};

#endif