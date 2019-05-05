#include "gtest/gtest.h"
#include <sstream>
#include <iostream>
#include "hts_Primers.h"

class Primer : public ::testing::Test {
    public:
        const std::string readData_1 = "@Read1\nTTTTTNGAAAAAAAAAGNTTTTT\n+\n#######################\n";
};

TEST_F(Primer, stub) {
    std::istringstream in1(readData_1);
    std::istringstream in2(readData_1);

    InputReader<PairedEndRead, PairedEndReadFastqImpl> ifp(in1, in2);

    while(ifp.has_next()) {
        auto i = ifp.next();
        PairedEndRead *per = dynamic_cast<PairedEndRead*>(i.get());
    }
};
