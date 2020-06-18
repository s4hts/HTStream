#include "gtest/gtest.h"
#include <sstream>
#include <iostream>
#include "hts_MinScreener.h"

class MinScreenerTest : public ::testing::Test {
public:
    const std::string phixTest = "ACTGACTGACTGACTGACTGACTGACTG";
    const std::string readData_1 = "@R1\nAAAAACTGACTGACTGTTTT\n+\nAAAAACTGACTGACTGTTTT\n";
    const size_t lookup_kmer_test = 2;
};

TEST_F(MinScreenerTest, check_fasta) {
}
