#include "gtest/gtest.h"
#include <sstream>
#include <iostream>
#include "overlapper.h"

class Overlapper : public ::testing::Test {
    public:
        const double misDensity = 0.25;        
        const size_t kmer = 8;
        const size_t kmerOffset = 2;
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

TEST_F(Overlapper, engulfR1) {
    std::istringstream in1(readData_engulf_r1_R1);
    std::istringstream in2(readData_engulf_r1_R2);

    InputReader<PairedEndRead, PairedEndReadFastqImpl> ifp(in1, in2);
    spReadBase rb;
    histVec tmp = nullptr;
	Counter counters;	
    setupCounter(counters);

    while(ifp.has_next()) {
        auto i = ifp.next();
        PairedEndRead *per = dynamic_cast<PairedEndRead*>(i.get());
        rb = check_read(*per, counters, misDensity, 10, tmp, false, 10, false, kmer, kmerOffset);
    }
    const Read &r = rb->get_read();
    ASSERT_EQ(r.get_seq(), "TAAGTATAAGGATAGATGAGATATGTCCAATTTGGTTAATGGTT");
};

TEST_F(Overlapper, engulfR2) {
    std::istringstream in1(readData_engulf_r2_R1);
    std::istringstream in2(readData_engulf_r2_R2);

    InputReader<PairedEndRead, PairedEndReadFastqImpl> ifp(in1, in2);
    spReadBase rb;
    histVec tmp = nullptr;
	Counter counters;	
    setupCounter(counters);


    while(ifp.has_next()) {
        auto i = ifp.next();
        PairedEndRead *per = dynamic_cast<PairedEndRead*>(i.get());
        rb = check_read(*per, counters, misDensity, 10, tmp, false, 10, false, kmer, kmerOffset);
    }
    const Read &r = rb->get_read();
    ASSERT_EQ(r.get_seq(), "GGTAAACCATTAACCAAATTGGACATATCTCATCT");
};


TEST_F(Overlapper, trim) {
    std::istringstream in1(readData_trim_overlap_R1);
    std::istringstream in2(readData_trim_overlap_R2);

    InputReader<PairedEndRead, PairedEndReadFastqImpl> ifp(in1, in2);
    spReadBase rb;
    histVec tmp = nullptr;
	Counter counters;	
    setupCounter(counters);


    while(ifp.has_next()) {
        auto i = ifp.next();
        PairedEndRead *per = dynamic_cast<PairedEndRead*>(i.get());
        rb = check_read(*per, counters, misDensity, 10, tmp, false, 10, false, kmer, kmerOffset);
    }
    const Read &r = rb->get_read();
    ASSERT_EQ(r.get_seq(), "ACCATAACATAAACC");
};


TEST_F(Overlapper, normal) {
    std::istringstream in1(readData_normal_overlap_R1);
    std::istringstream in2(readData_normal_overlap_R2);

    InputReader<PairedEndRead, PairedEndReadFastqImpl> ifp(in1, in2);
    spReadBase rb;
    histVec tmp = nullptr;
	Counter counters;	
    setupCounter(counters);


    while(ifp.has_next()) {
        auto i = ifp.next();
        PairedEndRead *per = dynamic_cast<PairedEndRead*>(i.get());
        rb = check_read(*per, counters, misDensity, 10, tmp, false, 10, false, kmer, kmerOffset);
    }
    const Read &r = rb->get_read();
    ASSERT_EQ(r.get_seq(), "GGTAAACCATTAACCAAATTGGACATATCTCATCTATCCTTATACTTAGCGCGAACGAGAAACGAGAGGAGACCGAGAG");
};


TEST_F(Overlapper, perfectOverlap) {
    std::istringstream in1(readData_perfect_overlap_R1);
    std::istringstream in2(readData_perfect_overlap_R2);

    InputReader<PairedEndRead, PairedEndReadFastqImpl> ifp(in1, in2);
    spReadBase rb;
    histVec tmp = nullptr;
	Counter counters;	
    setupCounter(counters);


    while(ifp.has_next()) {
        auto i = ifp.next();
        PairedEndRead *per = dynamic_cast<PairedEndRead*>(i.get());
        rb = check_read(*per, counters, misDensity, 10, tmp, false, 10, false, kmer, kmerOffset);
    }
    const Read &r = rb->get_read();
    ASSERT_EQ(r.get_seq(), "ACTTGACATTAAGCAAGTACCAGTACCGATACCATAGGACCCAAGGTA");
};

TEST_F(Overlapper, testHist) {
    std::istringstream in1(readData_perfect_overlap_R1);
    std::istringstream in2(readData_perfect_overlap_R2);

    InputReader<PairedEndRead, PairedEndReadFastqImpl> ifp(in1, in2);
    spReadBase rb;
    histVec tmp = histVec(new std::vector<unsigned long long int>);
    Counter counters;	
    setupCounter(counters);

    while(ifp.has_next()) {
        auto i = ifp.next();

     PairedEndRead *per = dynamic_cast<PairedEndRead*>(i.get());
        rb = check_read(*per, counters, misDensity, 10, tmp, false, 10, false, kmer, kmerOffset);
    }
    const Read &r = rb->get_read();
    std::string overlapOutput =  "ACTTGACATTAAGCAAGTACCAGTACCGATACCATAGGACCCAAGGTA";

    ASSERT_EQ((*tmp)[overlapOutput.length()], 1);
    ASSERT_EQ(r.get_seq(), overlapOutput);
};

