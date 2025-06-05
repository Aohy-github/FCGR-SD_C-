#include "../include/ThreadPool.h"

ThreadPool::ThreadPool(size_t n): stop(false){
    for(size_t i = 0 ; i < n ; i++){
        workers.emplace_back([this] (){
            while(true){
                std::function<void()> taks;
                {
                    std::unique_lock<std::mutex> lock(this->queMutex);
                    condition.wait(lock , [this](){
                        return this->stop || !this->Tasks.empty();
                    });
                    if(this->stop && this->Tasks.empty()) return;
                    taks = std::move(this->Tasks.front());
                    this->Tasks.pop();
                }
                taks();
            }
        });
    }
}