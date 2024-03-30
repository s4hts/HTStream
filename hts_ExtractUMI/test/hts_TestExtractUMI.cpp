#include "gtest/gtest.h"
#include <sstream>
#include <iostream>
#include "hts_ExtractUMI.h"

class ExtractUMITest : public ::testing::Test {
public:
    const std::string readData_1 = "@Read1\nTGACTTGACATTAAGCAAGTACCAGTACCGATACCATAGGACCCAAGGTA\n+\n##################################################\n";

    ExtractUMI eu;
};

TEST_F(ExtractUMITest, BasicExtract) {
    std::istringstream in1(readData_1);

    InputReader<SingleEndRead, SingleEndReadFastqImpl> ifp(in1);

    auto i = ifp.next();
    PairedEndRead *per = dynamic_cast<PairedEndRead*>(i.get());
    eu.extract_umi(per->non_const_read_one(), 6);
    ASSERT_EQ("@Read1_TGACTT", (per->non_const_read_one()).get_id_first());
};
