/*#include "gtest/gtest.h"
#include <sstream>
#include <iostream>
#include "phix_remover.h"

class SDTest : public ::testing::Test {
    public:
        const std::string readData_1 = "@Read1\nGGGGGGGGGGCAAAAAAAACGGGGGGGGGG\n+\n##############################\n";
        const std::string readData_2 = "@Read1\nAAAAAAAAAAAAAAAAAAAAAAAA\n+\n########################\n";
        size_t min_length = 5;
        size_t mega_min_length = 10000;
        size_t cut_size = 10;
};

TEST_F(SDTest, BasicTrim) {
    std::istringstream in1(readData_1);
    std::istringstream in2(readData_2);

    InputReader<PairedEndRead, PairedEndReadFastqImpl> ifp(in1, in2);

    while(ifp.has_next()) {
        auto i = ifp.next();
        PairedEndRead *per = dynamic_cast<PairedEndRead*>(i.get());
        Read &rb1 = per->non_const_read_one();
        Read &rb2 = per->non_const_read_two();
        rb1.setLCut(cut_size);
        rb1.setRCut(rb1.getLength() - cut_size);
        ASSERT_EQ("CAAAAAAAAC", (per->non_const_read_one()).get_sub_seq());
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
            Read &rb1 = per->non_const_read_one();
            Read &rb2 = per->non_const_read_two();

            rb1.setLCut(cut_size);
            rb1.setRCut(rb1.getLength() - cut_size);
            rb2.setLCut(cut_size);
            rb2.setRCut(rb2.getLength() - cut_size);
            per->checkDiscarded(min_length);
            writer_helper(per, tab, tab, true);
        }
    }
    ASSERT_EQ("Read1\tGTTTTTTTTG\t##########\n", out1->str());
};
*/
