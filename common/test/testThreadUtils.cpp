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
    long max = 10000;
    long block_size = 100;
    thread_pool threads;

    
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

              
