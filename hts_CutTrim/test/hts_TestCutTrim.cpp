#include "gtest/gtest.h"
#include <sstream>
#include <iostream>
#include "hts_CutTrim.h"

class CutTrim : public ::testing::Test {
    public:
        const std::string readData_1 = "@Read1\nTGACTTGACATTAAGCAAGTACCAGTACCGATACCATAGGACCCAAGGTA\n+\n##################################################\n";
        const std::string readData_2 = "@Read2\nGCTACCTTGGGTCCTATGGTATCGGTACTGGTACTTGCTTAATGTCAAGT\n+\n##################################################\n";

        size_t min_length = 5;
};

TEST_F(CutTrim, BasicTrim) {
    std::istringstream in1(readData_1);
    std::istringstream in2(readData_2);

    InputReader<PairedEndRead, PairedEndReadFastqImpl> ifp(in1, in2);

    while(ifp.has_next()) {
        auto i = ifp.next();
        PairedEndRead *per = dynamic_cast<PairedEndRead*>(i.get());
        cut_trim(per->non_const_read_one(), 5, 15, 0);
        ASSERT_EQ("TGACATTAAGCAAGTACCAGTACCGATACC", (per->non_const_read_one()).get_sub_seq());
    }
};

TEST_F(CutTrim, Stranded) {
    std::istringstream in1(readData_2); //Reverse these two so R1 is discarded and R2 is RCed
    std::istringstream in2(readData_1);

    InputReader<PairedEndRead, PairedEndReadFastqImpl> ifp(in1, in2);
    std::shared_ptr<std::ostringstream> out1(new std::ostringstream);
    {
        std::shared_ptr<HtsOfstream> hts_of(new HtsOfstream(out1));
        std::shared_ptr<OutputWriter> tab(new ReadBaseOutTab(hts_of));
        while(ifp.has_next()) {
            auto i = ifp.next();
            PairedEndRead *per = dynamic_cast<PairedEndRead*>(i.get());
            per->checkDiscarded(min_length);
            cut_trim(per->non_const_read_one(), 25, 26, 0);
            cut_trim(per->non_const_read_two(), 5, 15, 0);
            writer_helper(per, tab, tab, true);
        }
    }
    ASSERT_EQ("Read1\tGGTATCGGTACTGGTACTTGCTTAATGTCA\t##############################\n", out1->str());
};

TEST_F(CutTrim, MaxLength) {
    std::istringstream in1(readData_1);
    std::istringstream in2(readData_2);

    InputReader<PairedEndRead, PairedEndReadFastqImpl> ifp(in1, in2);
    std::shared_ptr<std::ostringstream> out1(new std::ostringstream);
    {
        std::shared_ptr<HtsOfstream> hts_of(new HtsOfstream(out1));
        std::shared_ptr<OutputWriter> tab(new ReadBaseOutTab(hts_of));
        while(ifp.has_next()) {
            auto i = ifp.next();
            PairedEndRead *per = dynamic_cast<PairedEndRead*>(i.get());
            per->checkDiscarded(min_length);
            cut_trim(per->non_const_read_one(), 0, 0, 10);
            cut_trim(per->non_const_read_two(), 0, 0, 10);
            writer_helper(per, tab, tab);
        }
    }
    ASSERT_EQ("Read1\tTGACTTGACA\t##########\tRead2\tGCTACCTTGG\t##########\n", out1->str());
};

TEST_F(CutTrim, Both) {
    std::istringstream in1(readData_1);
    std::istringstream in2(readData_2);

    InputReader<PairedEndRead, PairedEndReadFastqImpl> ifp(in1, in2);
    std::shared_ptr<std::ostringstream> out1(new std::ostringstream);
    {
        std::shared_ptr<HtsOfstream> hts_of(new HtsOfstream(out1));
        std::shared_ptr<OutputWriter> tab(new ReadBaseOutTab(hts_of));
        while(ifp.has_next()) {
            auto i = ifp.next();
            PairedEndRead *per = dynamic_cast<PairedEndRead*>(i.get());
            per->checkDiscarded(min_length);
            cut_trim(per->non_const_read_one(), 5, 10, 10);
            cut_trim(per->non_const_read_two(), 5, 10, 10);
            writer_helper(per, tab, tab);
        }
    }
    ASSERT_EQ("Read1\tTGACATTAAG\t##########\tRead2\tCTTGGGTCCT\t##########\n", out1->str());
};
