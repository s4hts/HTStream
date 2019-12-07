#include <gtest/gtest.h>
#include "threadutils.h"

#include <functional>
#include <iostream>
#include <numeric>
#include <sstream>


class ThreadTest : public ::testing::Test {
public:
    
};

TEST_F(ThreadTest, accumulate) {
    long max = 10000000;
    long block_size = 10000;
    thread_pool threads;

    
    std::vector<long> data_vec;
    for (long i = 0; i < max; ++i) {
        data_vec.push_back(i);
    }

    std::vector<std::future<long> > futures;
    for (auto i = data_vec.begin(); i < data_vec.end(); i+=block_size) {
        futures.push_back(threads.submit([i,block_size]() {
                    long b = 0;
                    // this loop is just to work the cpus
                    for (auto x = 0; x <= block_size; ++x) {
                        b+=x/std::accumulate(i, i+block_size, 0l);
                    }
                    return std::accumulate(i, i+block_size, b); }));
    }

    long sum = 0;
    for (auto i = futures.begin(); i < futures.end(); ++i) {
        sum += i->get();
    }

    ASSERT_EQ(sum, std::accumulate(data_vec.begin(), data_vec.end(), 0l));
}
              
TEST_F(ThreadTest, accumulate_queue_max) {
    long max = 10000;
    long block_size = 100;
    thread_pool threads(1);

    
    std::vector<long> data_vec;
    for (long i = 0; i < max; ++i) {
        data_vec.push_back(i);
    }

    std::vector<std::future<long> > futures;
    for (auto i = data_vec.begin(); i < data_vec.end(); i+=block_size) {
        futures.push_back(threads.submit([i,block_size]() { return static_cast<long>(std::accumulate(i, i+block_size, 0)); }));
    }

    long sum = 0;
    for (auto i = futures.begin(); i < futures.end(); ++i) {
        sum += i->get();
    }

    ASSERT_EQ(sum, std::accumulate(data_vec.begin(), data_vec.end(), 0));
}

              
