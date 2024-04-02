#include "gtest/gtest.h"
#include <sstream>
#include <iostream>
#include "hts_ExtractUMI.h"

class ExtractUMITest : public ::testing::Test {
public:

    const std::string readData_1 = "@Read1\nNAAAAAGACATTAAGCAA\n+\n!!!!!!############\n";
    const std::string readData_2 = "@Read2\nTTTTTTGACATTAAGCAA\n+\n!!!!!!############\n";

    ExtractUMI eu;
};

TEST_F(ExtractUMITest, BasicExtract) { // SE extract test
    std::istringstream in1(readData_1);

    InputReader<SingleEndRead, SingleEndReadFastqImpl> ifp(in1);

    auto i = ifp.next();
    SingleEndRead *ser = dynamic_cast<SingleEndRead*>(i.get());
    eu.extract_umi(ser->non_const_read_one(), 6, 0, 0, 33, false, false);
    ASSERT_EQ("Read1_NAAAAA", (ser->non_const_read_one()).get_id_first());
};


TEST_F(ExtractUMITest, ReadOneExtract) { // R1 extract test
    std::istringstream in1(readData_1); 
    std::istringstream in2(readData_2);

    InputReader<PairedEndRead, PairedEndReadFastqImpl> ifp(in1, in2);
    std::shared_ptr<std::ostringstream> out1(new std::ostringstream);
    std::shared_ptr<HtsOfstream> hts_of(new HtsOfstream(out1));
    std::shared_ptr<OutputWriter> tab(new ReadBaseOutTab(hts_of));

    WriterHelper writer(tab, tab);
    while(ifp.has_next()) {
        auto i = ifp.next();
        PairedEndRead *per = dynamic_cast<PairedEndRead*>(i.get());
        eu.extract_umi(per->non_const_read_one(), 6, 0, 0, 33, false, false);
        eu.extract_umi(per->non_const_read_two(), 6, 0, 0, 33, false, false);
        writer(*per);
    }
    ASSERT_EQ("Read1_NAAAAA\tGACATTAAGCAA\t############\tRead2_NAAAAA\tTTTTTTGACATTAAGCAA\t!!!!!!############\n", out1->str());
};


TEST_F(ExtractUMITest, ReadTwoExtract) { // R2 extract test
    std::istringstream in1(readData_1); 
    std::istringstream in2(readData_2);

    InputReader<PairedEndRead, PairedEndReadFastqImpl> ifp(in1, in2);
    std::shared_ptr<std::ostringstream> out1(new std::ostringstream);
    std::shared_ptr<HtsOfstream> hts_of(new HtsOfstream(out1));
    std::shared_ptr<OutputWriter> tab(new ReadBaseOutTab(hts_of));

    WriterHelper writer(tab, tab);
    while(ifp.has_next()) {
        auto i = ifp.next();
        PairedEndRead *per = dynamic_cast<PairedEndRead*>(i.get());
        eu.extract_umi(per->non_const_read_two(), 6, 0, 0, 33, false, false);
        eu.extract_umi(per->non_const_read_one(), 6, 0, 0, 33, false, false);
        writer(*per);
    }
    ASSERT_EQ("Read1_TTTTTT\tNAAAAAGACATTAAGCAA\t!!!!!!############\tRead2_TTTTTT\tGACATTAAGCAA\t############\n", out1->str());
};


TEST_F(ExtractUMITest, QualFilt) { // Hard Filter extract test
    std::istringstream in1(readData_1);

    InputReader<SingleEndRead, SingleEndReadFastqImpl> ifp(in1);

    auto i = ifp.next();
    SingleEndRead *ser = dynamic_cast<SingleEndRead*>(i.get());
    eu.extract_umi(ser->non_const_read_one(), 6, 20, 0, 33, false, false);
    ASSERT_EQ(true, (ser->non_const_read_one()).getDiscard());
};


TEST_F(ExtractUMITest, AvgQualFilt) { // Average Filter extract test
    std::istringstream in1(readData_1);

    InputReader<SingleEndRead, SingleEndReadFastqImpl> ifp(in1);

    auto i = ifp.next();
    SingleEndRead *ser = dynamic_cast<SingleEndRead*>(i.get());
    eu.extract_umi(ser->non_const_read_one(), 6, 0, 20, 33, false, false);
    ASSERT_EQ(true, (ser->non_const_read_one()).getDiscard());
};


TEST_F(ExtractUMITest, HomopolymerFilt) { // Homopolymer filter test
    std::istringstream in2(readData_2);

    InputReader<SingleEndRead, SingleEndReadFastqImpl> ifp(in2);

    auto i = ifp.next();
    SingleEndRead *ser = dynamic_cast<SingleEndRead*>(i.get());
    eu.extract_umi(ser->non_const_read_one(), 6, 0, 0, 33, true, false);
    ASSERT_EQ(true, (ser->non_const_read_one()).getDiscard());
};


TEST_F(ExtractUMITest, NFilt) { // N filter test
    std::istringstream in1(readData_1);

    InputReader<SingleEndRead, SingleEndReadFastqImpl> ifp(in1);

    auto i = ifp.next();
    SingleEndRead *ser = dynamic_cast<SingleEndRead*>(i.get());
    eu.extract_umi(ser->non_const_read_one(), 6, 0, 0, 33, false, true);
    ASSERT_EQ(true, (ser->non_const_read_one()).getDiscard());
};
