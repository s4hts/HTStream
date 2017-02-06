#include "gtest/gtest.h"
#include <sstream>
#include <iostream>
#include "nremover.h"

class TrimN : public ::testing::Test {
    public:
        const std::string readData_1 = "@Read1\nTTTTTNGAAAAAAAAAGNTTTTT\n+\n#######################\n";
        const std::string readData_2 = "@Read1\nAAAAAAAAAAAAAAAAAAAAAAAA\n+\n########################\n";
        const std::string readData_3 = "@Read1\nTTTTTNGGNTTTTTTTTTTTTTN\n+\n#######################\n";
        const std::string readData_4 = "@Read1\nNTTTTAGGATTTTTTTTTTTTTN\n+\n#######################\n";
        const std::string readData_5 = "@Read1\nATTTTAGGATTTTTTTTTTTTTN\n+\n#######################\n";
        const std::string readData_6 = "@Read1\nNTTTTAGGATTTTTTTTTTTTTA\n+\n#######################\n";
        size_t min_trim = 5;
        size_t min_length = 5;
        size_t max_mismatch = 3;
};

TEST_F(TrimN, EdgeRightN) {
    std::istringstream in1(readData_5);
    std::istringstream in2(readData_5);

    InputReader<PairedEndRead, PairedEndReadFastqImpl> ifp(in1, in2);

    while(ifp.has_next()) {
        auto i = ifp.next();
        PairedEndRead *per = dynamic_cast<PairedEndRead*>(i.get());
        trim_n(per->non_const_read_one(), min_trim);
        ASSERT_EQ("ATTTTAGGATTTTTTTTTTTTT", (per->non_const_read_one()).get_sub_seq());
    }
};


TEST_F(TrimN, EdgeLeftN) {
    std::istringstream in1(readData_6);
    std::istringstream in2(readData_6);

    InputReader<PairedEndRead, PairedEndReadFastqImpl> ifp(in1, in2);

    while(ifp.has_next()) {
        auto i = ifp.next();
        PairedEndRead *per = dynamic_cast<PairedEndRead*>(i.get());
        trim_n(per->non_const_read_one(), min_trim);
        ASSERT_EQ("TTTTAGGATTTTTTTTTTTTTA", (per->non_const_read_one()).get_sub_seq());
    }
};


TEST_F(TrimN, EdgeCaseNonEnds) {
    std::istringstream in1(readData_4);
    std::istringstream in2(readData_4);

    InputReader<PairedEndRead, PairedEndReadFastqImpl> ifp(in1, in2);

    while(ifp.has_next()) {
        auto i = ifp.next();
        PairedEndRead *per = dynamic_cast<PairedEndRead*>(i.get());
        trim_n(per->non_const_read_one(), min_trim);
        ASSERT_EQ("TTTTAGGATTTTTTTTTTTTT", (per->non_const_read_one()).get_sub_seq());
    }
};

TEST_F(TrimN, BasicTrim) {
    std::istringstream in1(readData_1);
    std::istringstream in2(readData_2);

    InputReader<PairedEndRead, PairedEndReadFastqImpl> ifp(in1, in2);

    while(ifp.has_next()) {
        auto i = ifp.next();
        PairedEndRead *per = dynamic_cast<PairedEndRead*>(i.get());
        trim_n(per->non_const_read_one(), min_trim);
        ASSERT_EQ("GAAAAAAAAAG", (per->non_const_read_one()).get_sub_seq());
    }
};

TEST_F(TrimN, NoTrim) {
    std::istringstream in1(readData_1);
    std::istringstream in2(readData_2);

    InputReader<PairedEndRead, PairedEndReadFastqImpl> ifp(in1, in2);

    while(ifp.has_next()) {
        auto i = ifp.next();
        PairedEndRead *per = dynamic_cast<PairedEndRead*>(i.get());
        trim_n(per->non_const_read_two(), min_trim);
        ASSERT_EQ("AAAAAAAAAAAAAAAAAAAAAAAA", (per->non_const_read_two()).get_sub_seq());
    }
};

TEST_F(TrimN, TwoBasicTrim) {
    std::istringstream in1(readData_3);
    std::istringstream in2(readData_2);

    InputReader<PairedEndRead, PairedEndReadFastqImpl> ifp(in1, in2);

    while(ifp.has_next()) {
        auto i = ifp.next();
        PairedEndRead *per = dynamic_cast<PairedEndRead*>(i.get());
        trim_n(per->non_const_read_one(), min_trim);
        ASSERT_EQ("TTTTTTTTTTTTT", (per->non_const_read_one()).get_sub_seq());
    }
};

