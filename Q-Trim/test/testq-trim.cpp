#include "gtest/gtest.h"
#include <sstream>
#include <iostream>
#include "q_trim.h"

class SDTest : public ::testing::Test {
    public:
        const std::string readData_1 = "@Read1\nTTTTTGGAAAAAAAAAGTCTTTGTTG\n+\n#####ZDDDDDDDDDDDDDZ######\n";
        const std::string readData_2 = "@Read1\nAAAAAAAAAAAAAAAAAAAAAAAA\n+\n########################\n";
        size_t min_length = 5;
        //+ 33 handled in main
        size_t min_qual = 20 + 33;
};

TEST_F(SDTest, BasicTrim) {
    std::istringstream in1(readData_1);
    std::istringstream in2(readData_2);

    InputReader<PairedEndRead, PairedEndReadFastqImpl> ifp(in1, in2);

    while(ifp.has_next()) {
        auto i = ifp.next();
        PairedEndRead *per = dynamic_cast<PairedEndRead*>(i.get());
        trim_left(per->non_const_read_one(), min_qual);
        trim_right(per->non_const_read_one(), min_qual);
        ASSERT_EQ("ZDDDDDDDDDDDDDZ", (per->non_const_read_one()).get_sub_qual());
    }
};

TEST_F(SDTest, AllTrim) {
    std::istringstream in1(readData_1);
    std::istringstream in2(readData_2);

    InputReader<PairedEndRead, PairedEndReadFastqImpl> ifp(in1, in2);

    while(ifp.has_next()) {
        auto i = ifp.next();
        PairedEndRead *per = dynamic_cast<PairedEndRead*>(i.get());
        trim_left(per->non_const_read_two(), min_qual);
        trim_right(per->non_const_read_two(), min_qual);
        ASSERT_EQ("", (per->non_const_read_two()).get_sub_seq());
    }
};

TEST_F(SDTest, Stranded) {
    std::istringstream in1(readData_2); //REverse these two so R1 is discarded and R2 is RCed
    std::istringstream in2(readData_1);

    InputReader<PairedEndRead, PairedEndReadFastqImpl> ifp(in1, in2);
    std::shared_ptr<std::ostringstream> out1(new std::ostringstream);

    {
        std::shared_ptr<HtsOfstream> hts_of(new HtsOfstream(out1));
        std::shared_ptr<OutputWriter> tab(new ReadBaseOutTab(hts_of));
        while(ifp.has_next()) {
            auto i = ifp.next();
            PairedEndRead *per = dynamic_cast<PairedEndRead*>(i.get());
            trim_left(per->non_const_read_one(), min_qual);
            trim_right(per->non_const_read_one(), min_qual);
            trim_left(per->non_const_read_two(), min_qual);
            trim_right(per->non_const_read_two(), min_qual);
            per->checkDiscarded(min_length);
            //per->checkDiscarded(min_length);
            //stranded (R2 will be RCed)
            writer_helper(per, tab, tab, true);
        }
    }
    ASSERT_EQ("Read1\tAGACTTTTTTTTTCC\tZDDDDDDDDDDDDDZ\n", out1->str());
};

