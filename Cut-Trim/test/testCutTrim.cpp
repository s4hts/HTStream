#include "gtest/gtest.h"
#include <sstream>
#include <iostream>
#include "cut_trim.h"

class CutTrim : public ::testing::Test {
    public:
        const std::string readData_1 = "@Read1\nGGGGGGGGGGCAAAAAAAACGGGGGGGGGG\n+\n##############################\n";
        const std::string readData_2 = "@Read1\nAAAAAAAAAAAAAAAAAAAAAAAA\n+\n########################\n";
        size_t min_length = 5;
        size_t mega_min_length = 10000;
        size_t cut_size = 10;
};

TEST_F(CutTrim, BasicTrim) {
    std::istringstream in1(readData_1);
    std::istringstream in2(readData_2);

    InputReader<PairedEndRead, PairedEndReadFastqImpl> ifp(in1, in2);

    while(ifp.has_next()) {
        auto i = ifp.next();
        PairedEndRead *per = dynamic_cast<PairedEndRead*>(i.get());
        cut_trim(per->non_const_read_one(), false, false, cut_size);
        ASSERT_EQ("CAAAAAAAAC", (per->non_const_read_one()).get_sub_seq());
    }
};

TEST_F(CutTrim, Stranded) {
    std::istringstream in1(readData_2); //REverse these two so R1 is discarded and R2 is RCed
    std::istringstream in2(readData_1);

    InputReader<PairedEndRead, PairedEndReadFastqImpl> ifp(in1, in2);
    std::shared_ptr<std::ostringstream> out1(new std::ostringstream);
    Counter c;
    {
        std::shared_ptr<HtsOfstream> hts_of(new HtsOfstream(out1));
        std::shared_ptr<OutputWriter> tab(new ReadBaseOutTab(hts_of));
        while(ifp.has_next()) {
            auto i = ifp.next();
            PairedEndRead *per = dynamic_cast<PairedEndRead*>(i.get());
            per->checkDiscarded(min_length);
            cut_trim(per->non_const_read_one(), false, false, cut_size);
            cut_trim(per->non_const_read_two(), false, false, cut_size);
            writer_helper(per, tab, tab, true, c);
        }
    }
    ASSERT_EQ("Read1\tGTTTTTTTTG\t##########\n", out1->str());
};

