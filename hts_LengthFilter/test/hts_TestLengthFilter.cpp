#include "gtest/gtest.h"
#include <sstream>
#include <iostream>
#include "hts_LengthFilter.h"

class LengthFilterTest : public ::testing::Test {
public:
    const std::string readData_1 = "@Read1\nTGACTTGACATTAAGCAAGTACCAGTACCGATACCATAGGACCCAAGGTA\n+\n##################################################\n";
    const std::string readData_2 = "@Read2\nGTCCTATGGT\n+\n##########\n";

    size_t min_length = 5;
    LengthFilter lf;
};


TEST_F(LengthFilterTest, MaxLength) {
    std::istringstream in1(readData_1);
    std::istringstream in2(readData_2);

    InputReader<PairedEndRead, PairedEndReadFastqImpl> ifp(in1, in2);
    std::shared_ptr<std::ostringstream> out1(new std::ostringstream);
    std::shared_ptr<HtsOfstream> hts_of(new HtsOfstream(out1));
    std::shared_ptr<OutputWriter> tab(new ReadBaseOutTab(hts_of));

    WriterHelper writer(tab, tab, false, false);
    while(ifp.has_next()) {
        auto i = ifp.next();
        PairedEndRead *per = dynamic_cast<PairedEndRead*>(i.get());
        lf.length_filter(per->non_const_read_one(), 0, 8);
        lf.length_filter(per->non_const_read_two(), 0, 11);
        writer(*per);
    }
    ASSERT_EQ("Read2\tGTCCTATGGT\t##########\n", out1->str());
};


TEST_F(LengthFilterTest, MinLength) {
    std::istringstream in1(readData_1);
    std::istringstream in2(readData_2);

    InputReader<PairedEndRead, PairedEndReadFastqImpl> ifp(in1, in2);
    std::shared_ptr<std::ostringstream> out1(new std::ostringstream);
    std::shared_ptr<HtsOfstream> hts_of(new HtsOfstream(out1));
    std::shared_ptr<OutputWriter> tab(new ReadBaseOutTab(hts_of));

    WriterHelper writer(tab, tab, false, false);
    while(ifp.has_next()) {
        auto i = ifp.next();
        PairedEndRead *per = dynamic_cast<PairedEndRead*>(i.get());
        lf.length_filter(per->non_const_read_one(), 11, 100);
        lf.length_filter(per->non_const_read_two(), 11, 100);
        writer(*per);
    }
    ASSERT_EQ("Read1\tTGACTTGACATTAAGCAAGTACCAGTACCGATACCATAGGACCCAAGGTA\t##################################################\n", out1->str());
};



TEST_F(LengthFilterTest, MaxLengthStranded) {
    std::istringstream in1(readData_1);
    std::istringstream in2(readData_2);

    InputReader<PairedEndRead, PairedEndReadFastqImpl> ifp(in1, in2);
    std::shared_ptr<std::ostringstream> out1(new std::ostringstream);
    std::shared_ptr<HtsOfstream> hts_of(new HtsOfstream(out1));
    std::shared_ptr<OutputWriter> tab(new ReadBaseOutTab(hts_of));
    WriterHelper writer(tab, tab, true, false);
    while(ifp.has_next()) {
        auto i = ifp.next();
        PairedEndRead *per = dynamic_cast<PairedEndRead*>(i.get());
        lf.length_filter(per->non_const_read_one(), 0, 11);
        lf.length_filter(per->non_const_read_two(), 0, 11);
        writer(*per);
    }
    ASSERT_EQ("Read2\tACCATAGGAC\t##########\n", out1->str());
};


TEST_F(LengthFilterTest, MaxLengthOrphan) {
    std::istringstream in1(readData_1);
    std::istringstream in2(readData_2);

    InputReader<PairedEndRead, PairedEndReadFastqImpl> ifp(in1, in2);
    std::shared_ptr<std::ostringstream> out1(new std::ostringstream);
    std::shared_ptr<HtsOfstream> hts_of(new HtsOfstream(out1));
    std::shared_ptr<OutputWriter> tab(new ReadBaseOutTab(hts_of));
    WriterHelper writer(tab, tab, false, true);
    while(ifp.has_next()) {
        auto i = ifp.next();
        PairedEndRead *per = dynamic_cast<PairedEndRead*>(i.get());
        lf.length_filter(per->non_const_read_one(), 0, 11);
        lf.length_filter(per->non_const_read_two(), 0, 11);
        writer(*per);
    }
    ASSERT_EQ("", out1->str());
};
