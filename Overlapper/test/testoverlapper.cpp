#include "gtest/gtest.h"
#include <sstream>
#include <iostream>
#include "overlapper.h"

class Overlapper : public ::testing::Test {
    public:
        
        const std::string readData_perfect_overlap_R1 = "@R1_perfect\nACTTGACATTAAGCAAGTACCAGTACCGATACCATAGGACCCAAGGTA\n+\n111111111111111111111111111111111111111111111111\n";
        const std::string readData_perfect_overlap_R2 = "@R2_perfect\nTACCTTGGGTCCTATGGTATCGGTACTGGTACTTGCTTAATGTCAAGT\n+\n111111111111111111111111111111111111111111111111\n";
        
        const std::string readData_normal_overlap_R1 = "@R1_\nGGTAAACCATTAACCAAATTGGACATATCTCATCTATCCTTATACTTA\n+\n111111111111111111111111111111111111111111111111\n";
        const std::string readData_normal_overlap_R2 = "@R2_\nCTCTCGGTCTCCTCTCGTTTCTCGTTCGCGCTAAGTATAAGGATAGA\n+\n11111111111111111111111111111111111111111111111\n";
        
        const std::string readData_trim_overlap_R1 = "@R1_\nACCATAACATAAACCATTAACCAAATTGGACATATCTCATCTATCCTTATACTTA\n+\n1111111111111111111111111111111111111111111111111111111\n";
        const std::string readData_trim_overlap_R2 = "@R2_\nGGTTTATGTTATGGTAATATAGTATAGAGTATAGTTGCGTC\n+\n11111111111111111111111111111111111111111\n";

        const std::string readData_engulf_r2_R1 = "@R1_\nGGTAAACCATTAACCAAATTGGACATATCTCATCTATCCTTATACTTA\n+\n111111111111111111111111111111111111111111111111\n";
        const std::string readData_engulf_r2_R2 = "@R2_\nAGATGAGATATGTCCAATTTGGTTAATGGTT\n+\n1111111111111111111111111111111";
        
        const std::string readData_engulf_r1_R1 = "@R1_\nAACCATTAACCAAATTGGACATATCTCATCT\n+\n1111111111111111111111111111111\n";
        const std::string readData_engulf_r1_R2 = "@R2_\nTAAGTATAAGGATAGATGAGATATGTCCAATTTGGTTAATGGTTTACC\n+\n111111111111111111111111111111111111111111111111\n";
};

TEST_F(Overlapper, perfectOverlap) {
    std::istringstream in1(readData_perfect_overlap_R1);
    std::istringstream in2(readData_perfect_overlap_R2);

    InputReader<PairedEndRead, PairedEndReadFastqImpl> ifp(in1, in2);
    spReadBase rb;

    while(ifp.has_next()) {
        auto i = ifp.next();
        PairedEndRead *per = dynamic_cast<PairedEndRead*>(i.get());
        rb = check_read(*per, 100, 100);
    }
    if (rb != nullptr) {
        std::cout << "HERE\n";
        const Read &r = rb->get_read();
        std::cout << "HERE\n";
        std::cout << r.get_id() << '\n';
        std::cout << "HERE\n";
    }
};

