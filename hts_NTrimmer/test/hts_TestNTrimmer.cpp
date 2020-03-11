#include "gtest/gtest.h"
#include <sstream>
#include <iostream>
#include "hts_NTrimmer.h"

class TrimN : public ::testing::Test {
public:
    const std::string readData_1 = "@Read1\nTTTTTNGAAAAAAAAAGNTTTTT\n+\n#######################\n";
    const std::string readData_2 = "@Read1\nAAAAAAAAAAAAAAAAAAAAAAAA\n+\n########################\n";
    const std::string readData_3 = "@Read1\nTTTTTNGGNTTTTTTTTTTTTTN\n+\n#######################\n";
    const std::string readData_4 = "@Read1\nNTTTTAGGATTTTTTTTTTTTTN\n+\n#######################\n";
    const std::string readData_5 = "@Read1\nATTTTAGGATTTTTTTTTTTTTN\n+\n#######################\n";
    const std::string readData_6 = "@Read1\nNTTTTAGGATTTTTTTTTTTTTA\n+\n#######################\n";
    const std::string readData_7 = "@Read1\nGTTTTAGGATTNTTTTTTTTTTA\n+\n#######################\n";
    const std::string readData_8 = "@Read1\nNNNNNNNNNNNNNNNNNNNNNNN\n+\n#######################\n";
    const std::string readData_9a = "@Read1\nCTGACTGACTGANNNACTGACTGACTGNCTGACTG\n+\n###################################\n";
    const std::string readData_9b = "@Read1\nCTGACTGACTGANNNACTGACTGACTGANTGACTG\n+\n###################################\n";
    NTrimmer nt;
};

TEST_F(TrimN, Exclude) {
    std::istringstream in1(readData_1);
    std::istringstream in2(readData_5);

    InputReader<PairedEndRead, PairedEndReadFastqImpl> ifp(in1, in2);

    while(ifp.has_next()) {
        auto i = ifp.next();
        PairedEndRead *per = dynamic_cast<PairedEndRead*>(i.get());
        nt.trim_n(per->non_const_read_one(), true);
        ASSERT_EQ("", (per->non_const_read_one()).get_sub_seq());
    }
};

TEST_F(TrimN, EdgeRightN) {
    std::istringstream in1(readData_5);
    std::istringstream in2(readData_5);

    InputReader<PairedEndRead, PairedEndReadFastqImpl> ifp(in1, in2);

    while(ifp.has_next()) {
        auto i = ifp.next();
        PairedEndRead *per = dynamic_cast<PairedEndRead*>(i.get());
        nt.trim_n(per->non_const_read_one(), false);
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
        nt.trim_n(per->non_const_read_one(), false);
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
        nt.trim_n(per->non_const_read_one(), false);
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
        nt.trim_n(per->non_const_read_one(), false);
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
        nt.trim_n(per->non_const_read_two(), false);
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
        nt.trim_n(per->non_const_read_one(), false);
        ASSERT_EQ("TTTTTTTTTTTTT", (per->non_const_read_one()).get_sub_seq());
    }
};

TEST_F(TrimN, equalTrim) {
    std::istringstream in1(readData_7);

    InputReader<SingleEndRead, SingleEndReadFastqImpl> ifp(in1);

    while(ifp.has_next()) {
        auto i = ifp.next();
        SingleEndRead *per = dynamic_cast<SingleEndRead*>(i.get());
        nt.trim_n(per->non_const_read_one(), false);
        ASSERT_EQ("GTTTTAGGATT", (per->non_const_read_one()).get_sub_seq());
    }
};

TEST_F(TrimN, allN) {
    std::istringstream in1(readData_8);

    InputReader<SingleEndRead, SingleEndReadFastqImpl> ifp(in1);

    while(ifp.has_next()) {
        auto i = ifp.next();
        SingleEndRead *per = dynamic_cast<SingleEndRead*>(i.get());
        nt.trim_n(per->non_const_read_one(), false);
        ASSERT_EQ("N", (per->non_const_read_one()).get_sub_seq());
    }
};

TEST_F(TrimN, longN) {
    std::istringstream in1(readData_9a);
    std::istringstream in2(readData_9b);

    InputReader<PairedEndRead, PairedEndReadFastqImpl> ifp(in1, in2);

    while(ifp.has_next()) {
        auto i = ifp.next();
        PairedEndRead *per = dynamic_cast<PairedEndRead*>(i.get());
        nt.trim_n(per->non_const_read_one(), false);
        nt.trim_n(per->non_const_read_two(), false);
        ASSERT_EQ("CTGACTGACTGA", (per->non_const_read_one()).get_sub_seq());
        ASSERT_EQ("ACTGACTGACTGA", (per->non_const_read_two()).get_sub_seq());
    }
};
