#include "gtest/gtest.h"
#include <sstream>
#include <iostream>
#include <boost/program_options.hpp>
#include "hts_Stats.h"

namespace po = boost::program_options;

class StatsTest : public ::testing::Test {
    public:
        po::variables_map vm;
        const std::string readData_1 = "@Read1\nACTGAC\n+\nI#I#AH\n";
        const std::string readData_2 = "@Read2\nACTGAC\n+\nI#IIDH\n";
        const std::string readData_empty = "@Read0\n\n+\n\n";
};

TEST_F(StatsTest, BasicTrim) {
    std::istringstream in1(readData_1);
    std::istringstream in2(readData_2);

    InputReader<PairedEndRead, PairedEndReadFastqImpl> ifp(in1, in2);
    StatsCounters counters("hts_Stats", vm);
    while(ifp.has_next()) {
        auto i = ifp.next();
        PairedEndRead *per = dynamic_cast<PairedEndRead*>(i.get());
        counters.input(*per);
        counters.output(*per);
    }
    ASSERT_EQ(4u, counters.A);
    ASSERT_EQ(2u, counters.T);
    ASSERT_EQ(4u, counters.C);
    ASSERT_EQ(2u, counters.G);
    ASSERT_EQ(4u, counters.R1_bQ30);
    ASSERT_EQ(5u, counters.R2_bQ30);
    ASSERT_EQ(1u, counters.R1_bases[1][1]);
    ASSERT_EQ(0u, counters.R1_bases[1][0]);
    ASSERT_EQ(1u, counters.R1_bases[4][0]);
    ASSERT_EQ(0u, counters.R1_bases[4][4]);
    ASSERT_EQ(1u, counters.R1_qualities[1][2]);
    ASSERT_EQ(0u, counters.R1_qualities[1][42]);
    ASSERT_EQ(1u, counters.R1_qualities[4][32]);
    ASSERT_EQ(0u, counters.R1_qualities[4][42]);
};

TEST_F(StatsTest, ZeroLengthRead) {
    std::istringstream in1(readData_empty);

    InputReader<SingleEndRead, SingleEndReadFastqImpl> ifs(in1);
    StatsCounters counters("hts_Stats", vm);

    ASSERT_TRUE(ifs.has_next());
    auto i = ifs.next();
    SingleEndRead *ser = dynamic_cast<SingleEndRead*>(i.get());
    ASSERT_NE(nullptr, ser);

    counters.input(*ser);
    counters.output(*ser);

    ASSERT_EQ(1u, counters.SE_In);
    ASSERT_EQ(1u, counters.SE_Out);
    ASSERT_EQ(1u, counters.TotalFragmentsInput);
    ASSERT_EQ(1u, counters.TotalFragmentsOutput);
    ASSERT_EQ(0u, counters.SE_BpLen_In);
    ASSERT_EQ(0u, counters.TotalBasepairsInput);
    ASSERT_EQ(1u, counters.SE_Length.size());
    ASSERT_EQ(1u, counters.SE_Length[0]);
    ASSERT_TRUE(counters.SE_bases.empty());
    ASSERT_TRUE(counters.SE_qualities.empty());
}
