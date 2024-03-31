#include "gtest/gtest.h"
#include <sstream>
#include <iostream>
#include "hts_ExtractUMI.h"

class ExtractUMITest : public ::testing::Test {
public:

    const std::string readData_1 = "@Read1\nAAAAAAGACATTAAGCAA\n+\n!!!!!!############\n";
    const std::string readData_2 = "@Read2\nTTTTTTGACATTAAGCAA\n+\n!!!!!!############\n";

    ExtractUMI eu;
};

TEST_F(ExtractUMITest, BasicExtract) { // SE extract test
    std::istringstream in1(readData_1);

    InputReader<SingleEndRead, SingleEndReadFastqImpl> ifp(in1);

    auto i = ifp.next();
    SingleEndRead *ser = dynamic_cast<SingleEndRead*>(i.get());
    eu.extract_umi(ser->non_const_read_one(), 6, 0, 33);
    ASSERT_EQ("Read1_AAAAAA", (ser->non_const_read_one()).get_id_first());
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
        eu.extract_umi(per->non_const_read_one(), 6, 0, 33);
        eu.extract_umi(per->non_const_read_two(), 6, 0, 33);
        writer(*per);
    }
    ASSERT_EQ("Read1_AAAAAA\tGACATTAAGCAA\t############\tRead2_AAAAAA\tTTTTTTGACATTAAGCAA\t!!!!!!############\n", out1->str());
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
        eu.extract_umi(per->non_const_read_two(), 6, 0, 33);
        eu.extract_umi(per->non_const_read_one(), 6, 0, 33);
        writer(*per);
    }
    ASSERT_EQ("Read1_TTTTTT\tAAAAAAGACATTAAGCAA\t!!!!!!############\tRead2_TTTTTT\tGACATTAAGCAA\t############\n", out1->str());
};


TEST_F(ExtractUMITest, FilteredExtract) { // Filtered extract teset
    std::istringstream in1(readData_1);

    InputReader<SingleEndRead, SingleEndReadFastqImpl> ifp(in1);

    auto i = ifp.next();
    SingleEndRead *ser = dynamic_cast<SingleEndRead*>(i.get());
    eu.extract_umi(ser->non_const_read_one(), 6, 20, 33);
    ASSERT_EQ(true, (ser->non_const_read_one()).getDiscard());
};
