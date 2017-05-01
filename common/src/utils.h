#ifndef UTILS_H
#define UTILS_H

#include "ioHandler.h"
#include <unordered_map>
#include <string>
#include <iostream>
#include <fstream>
#include <memory>
#include <vector>
#include <boost/program_options.hpp>
#include <boost/bind.hpp>
#include <boost/thread/thread.hpp>
#include "version.h"
#include "typedefs.h"

namespace po = boost::program_options;


void setupCounter(Counter &c);
void write_stats(const std::string &statsFile, const bool &appendStats, const Counter &c, const std::string &program_name);

void outputWriters(std::shared_ptr<OutputWriter> &pe, std::shared_ptr<OutputWriter> &se, bool fastq_out, bool tab_out, bool interleaved_out, bool unmapped_out,  bool force, bool gzip_out, bool std_out, std::string &prefix);

void version_or_help(std::string program_name, po::options_description &desc, po::variables_map vm);

void setDefaultParams(po::options_description &desc, std::string program_name);
void setDefaultParamsCutting(po::options_description &desc);
void setDefaultParamsTrim(po::options_description &desc);

class Thread_Pool {
public:
    Thread_Pool(int _thread_count) : thread_count(_thread_count) { }

    void done() {
        ++thread_count;
        thread_limit_cv.notify_one();
    }

    template <typename FUNC, typename WRITER>
    void worker_function(FUNC f, WRITER w) {

        while(thread_count <= 0) {
            boost::mutex::scoped_lock scoped_lock(thread_limit_mutex);
            if (scoped_lock) {
                thread_limit_cv.wait(scoped_lock);
            }
        }
        --thread_count;
        boost::thread *t = new boost::thread(  [this, f, w] {
            auto i = f();
            {
                boost::unique_lock<boost::shared_mutex> lock(io_mutex);
                if (lock) {
                    w(i);
                    done();
                }
            }
        } ) ;
        g.add_thread(t);
    }

    void wait() {
        g.join_all();
    }

    ~Thread_Pool() {
        wait();
    }

private:
    int thread_count = 0;
    boost::mutex thread_limit_mutex;
    boost::shared_mutex io_mutex;
    boost::condition_variable thread_limit_cv;
    boost::thread_group g;

};





#endif
