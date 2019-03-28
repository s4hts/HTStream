#include "gtest/gtest.h"
#include <sstream>
#include <iostream>
#include "hts_Stats.h"

class Stats : public ::testing::Test {
    public:
        const std::string readData_1 = "@Read1\nACTG\n+\nI#I#\n";
        const std::string readData_2 = "@Read2\nACTG\n+\nI#II\n";
};

TEST_F(Stats, BasicTrim) {
    std::istringstream in1(readData_1);
    std::istringstream in2(readData_2);

    InputReader<PairedEndRead, PairedEndReadFastqImpl> ifp(in1, in2);
    StatsCounters counters("/dev/null", true, false, "stats", "");
    while(ifp.has_next()) {
        auto i = ifp.next();
        PairedEndRead *per = dynamic_cast<PairedEndRead*>(i.get());
        counters.input(*per);
        counters.output(*per);
    }
    ASSERT_EQ(2, counters.A);
    ASSERT_EQ(2, counters.T);
    ASSERT_EQ(2, counters.C);
    ASSERT_EQ(2, counters.G);
    ASSERT_EQ(2, counters.R1_bQ30);
    ASSERT_EQ(3, counters.R2_bQ30);
};
