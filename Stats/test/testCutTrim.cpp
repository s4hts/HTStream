#include "gtest/gtest.h"
#include <sstream>
#include <iostream>
#include "stats.h"

class Stats : public ::testing::Test {
    public:
        const std::string readData_1 = "@Read1\nACTG\n+\n####\n";
        const std::string readData_2 = "@Read2\nACTG\n+\n####\n";
};

TEST_F(Stats, BasicTrim) {
    std::istringstream in1(readData_1);
    std::istringstream in2(readData_2);

    InputReader<PairedEndRead, PairedEndReadFastqImpl> ifp(in1, in2);
    StatsCounters counters;
    while(ifp.has_next()) {
        auto i = ifp.next();
        PairedEndRead *per = dynamic_cast<PairedEndRead*>(i.get());
        counters.input(*per);
        counters.output(*per);
    }
    ASSERT_EQ(2, counters.c["A"]);
    ASSERT_EQ(2, counters.c["T"]);
    ASSERT_EQ(2, counters.c["C"]);
    ASSERT_EQ(2, counters.c["G"]);
    ASSERT_EQ(8, counters.c["R1_Q_Sum"]);
    ASSERT_EQ(8, counters.c["R2_Q_Sum"]);
};

