#ifndef UTILS_H
#define UTILS_H

#include "ioHandler.h"
#include <unordered_map>
#include <string>
#include <iostream>
#include <fstream>
#include <memory>
#include <vector>
#include <pthread.h>
#include <boost/program_options.hpp>
#include <boost/bind.hpp>
#include <boost/thread/thread.hpp>
#include "version.h"
#include "typedefs.h"
#include <queue>
#include <future>
namespace po = boost::program_options;


void setupCounter(Counter &c);
void write_stats(const std::string &statsFile, const bool &appendStats, const Counter &c, const std::string &program_name);

void outputWriters(std::shared_ptr<OutputWriter> &pe, std::shared_ptr<OutputWriter> &se, bool fastq_out, bool tab_out, bool interleaved_out, bool unmapped_out,  bool force, bool gzip_out, bool std_out, std::string &prefix);

void version_or_help(std::string program_name, po::options_description &desc, po::variables_map vm);

void setDefaultParams(po::options_description &desc, std::string program_name);
void setDefaultParamsCutting(po::options_description &desc);
void setDefaultParamsTrim(po::options_description &desc);

template <typename READ, typename FUNC, typename WRITER>
class Thread_Pool {
public:

    Thread_Pool(int _thread_count, FUNC _worker, WRITER _write) : thread_count(_thread_count), worker_function(_worker), writer_function(_write) { 
        finished = false;
        for (int i = 0; i < _thread_count; ++i) {
            boost::thread *t = new boost::thread(&Thread_Pool::worker_helper_thread, this); 
            g.add_thread(t);
        }
        boost::thread *t = new boost::thread(&Thread_Pool::writer_helper_thread, this); 
        g.add_thread(t);
    }

    void worker_helper_thread() {
        while (!finished || !protected_empty()) {
             
            boost::mutex::scoped_lock lock(queue_protect);
            if ( worker_queue.empty() ) {
                avail_data.wait(lock);
            } else {
                READ r = worker_queue.front();
                worker_queue.pop();
                lock.unlock();
                avail_data.notify_one();
                auto i = worker_function(r);
                
                {
                    boost::mutex::scoped_lock lock_io(io_mutex);
                    writer_queue.push(boost::bind(writer_function, i, r));
                    lock_io.unlock();
                    avail_io.notify_one();
                }
            }
        }   
    }

    void writer_helper_thread() {
        while (!finished || !protected_empty_io() || !protected_empty()  ) {
            boost::mutex::scoped_lock lock_io( io_mutex);
            if (writer_queue.empty()) {
                avail_io.wait(lock_io);
            } else {
                (writer_queue.front())();
                writer_queue.pop();
                lock_io.unlock();
            }
        }
    } 

    bool protected_empty() {
        boost::mutex::scoped_lock lock(queue_protect);
        return worker_queue.empty();
    }

    bool protected_empty_io() {
        boost::mutex::scoped_lock lock_io(io_mutex);
        return writer_queue.empty();
    }

    void push(READ r) {
        boost::mutex::scoped_lock lock(io_mutex);
        worker_queue.push( r );
        lock.unlock();
        avail_data.notify_one(); 
    }

    void wait() {
        finished = true;
        while ( !protected_empty() ) {
            avail_data.notify_all();
        }
        while ( !protected_empty_io() ) {
            avail_io.notify_all();
        }
        g.join_all();
    }

    ~Thread_Pool() {
        finished = true;
        wait();
    }

private:
    std::atomic<bool> finished;
    int thread_count;
    boost::mutex io_mutex;
    boost::mutex queue_protect;
     
    boost::condition_variable avail_data;
    boost::condition_variable avail_io;
    boost::thread_group g;
    std::queue<READ> worker_queue;
    std::queue<boost::function<void()> > writer_queue;
    FUNC worker_function;
    WRITER writer_function;
};






#endif
