#ifndef HTS_THREADUTILS_H
#define HTS_THREADUTILS_H

#include <atomic>
#include <condition_variable>
#include <functional>
#include <future>
#include <mutex>
#include <queue>
#include <thread>
#include <vector>

// inspired from C++ Concurency in Action by Anythony Williams

template <typename T>
class threadsafe_queue {
private:
    std::atomic_bool done;
    mutable std::mutex mut;
    std::queue<std::shared_ptr<T> > data_queue;
    std::condition_variable data_cond;
public:
    threadsafe_queue(): done(false)
    {}

    bool is_done() {
        return done;
    }

    void set_done() {
        if (!done) {
            done = true;
            data_cond.notify_all();
        }
    }
    
    void wait_and_pop(T& value) {
        std::unique_lock<std::mutex> lk(mut);
        data_cond.wait(lk, [this] { return !data_queue.empty() || done; });
        if (done) {
            return;
        }
        value = std::move(*data_queue.front());
        data_queue.pop();
    }

    bool try_pop(T& value) {
        std::lock_guard<std::mutex> lk(mut);
        if(data_queue.empty() || done) {
            return false;
        }
        value = std::move(*data_queue.front());
        data_queue.pop();
        return true;
    }
    
    std::shared_ptr<T> wait_and_pop() {
        std::unique_lock<std::mutex> lk(mut);
        data_cond.wait(lk, [this] { return !data_queue.empty() || done; });
        if (done) {
            return std::shared_ptr<T>();
        }
        std::shared_ptr<T> res=data_queue.front();
        data_queue.pop();
        return res;
    }

    std::shared_ptr<T> try_pop() {
        std::lock_guard<std::mutex> lk(mut);
        if(data_queue.empty() || done) {
            return std::shared_ptr<T>();
        }

        std::shared_ptr<T> res=data_queue.front();
        data_queue.pop();
        return res;
    }

    void push(T&& new_value) {
        std::shared_ptr<T> data( std::make_shared<T>(std::forward<T>(new_value)) );
        std::lock_guard<std::mutex> lk(mut);
        data_queue.push(data);
        data_cond.notify_one();
    }

    bool empty() const {
        std::lock_guard<std::mutex> lk(mut);
        return data_queue.empty();
    }
};

class thread_guard
{
    std::thread& t;
    const std::function<void()> *func;
public:
    explicit thread_guard(std::thread& t_, const std::function<void()> &func_):
        t(t_), func(&func_)
    {}
    explicit thread_guard(std::thread& t_):
        t(t_), func(NULL)
    {}


    ~thread_guard() {
        if(func) {
            (*func)();
        }
        if(t.joinable()) {
            t.join();
        }
    }
    thread_guard(thread_guard const&) = delete;
    thread_guard& operator=(thread_guard const&) = delete;
};

class function_wrapper
{
    struct impl_base {
        virtual void call() = 0;
        virtual ~impl_base() {}
    };
    std::unique_ptr<impl_base> impl;

    template<typename F>
    struct impl_type: impl_base {
        F f;
        impl_type(F&& f_): f(std::move(f_)) {}
        void call() { f(); }
    };

public:
    template<typename F>
    function_wrapper(F&& f):
        impl(new impl_type<F>(std::move(f)))
    {}

    void operator()() { impl->call(); }

    function_wrapper() = default;

    function_wrapper(function_wrapper&& other):
        impl(std::move(other.impl))
    {}

    function_wrapper& operator=(function_wrapper&& other) {
        impl = std::move(other.impl);
        return *this;
    }

    function_wrapper(const function_wrapper&) = delete;
    function_wrapper(function_wrapper&) = delete;
    function_wrapper& operator=(const function_wrapper&) = delete;
};

class join_threads {
    std::vector<std::thread>& threads;
    threadsafe_queue<function_wrapper> & queue;
    
public:
    explicit join_threads(std::vector<std::thread>& threads_, threadsafe_queue<function_wrapper>& queue_):
        threads(threads_), queue(queue_)
    {}
    
    ~join_threads() {
        queue.set_done();
        for(unsigned long i=0;i<threads.size();++i)
        {
            if(threads[i].joinable())
                threads[i].join();
        }
    }
};

class thread_pool
{
    std::atomic_bool done;
    threadsafe_queue<function_wrapper> work_queue;
    std::vector<std::thread> threads;

    // note joiner must be the last member variable so it is destroyed first
    join_threads joiner;

    void worker_thread() {
        while(!done) {
            function_wrapper task;
            work_queue.wait_and_pop(task);
            if (!done) {
                task();
            }
        }
    }

public:
    template<typename FunctionType> std::future<typename std::result_of<FunctionType()>::type> submit(FunctionType f) {
        typedef typename std::result_of<FunctionType()>::type result_type;
        std::packaged_task<result_type()> task(std::move(f));
        std::future<result_type> res(task.get_future());
        work_queue.push(std::move(task));
        return res;
    }

    thread_pool(): done(false), joiner(threads, work_queue)
    {
        // todo add this as prameter
        unsigned const thread_count = 4;//std::thread::hardware_concurrency();

        try {
            for (unsigned i = 0; i < thread_count; ++i) {
                threads.push_back(std::thread(&thread_pool::worker_thread, this));
            }
        } catch (...) {
            done = true;
            throw;
        }
    }

    ~thread_pool() {
        done = true;
    }
};
        
#endif
