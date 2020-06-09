#include "gtest/gtest.h"
#include <sstream>
#include <iostream>
#include "hts_Overlapper.h"

class OverlapperTest : public ::testing::Test {
    public:
        const double misDensity = 0.25;        
        const size_t kmer = 8;
        const size_t kmerOffset = 1;
        const std::string readData_perfect_overlap_R1 = "@R1_perfect\nACTTGACATTAAGCAAGTACCAGTACCGATACCATAGGACCCAAGGTA\n+\n111111111111111111111111111111111111111111111111\n";
        const std::string readData_perfect_overlap_R1_short_R1 = "@R1_perfect\nACTTGACATTAAGCAAGTACCAGTACCGATACCATAGGACCCAAG\n+\n111111111111111111111111111111111111111111111\n";
        const std::string readData_perfect_overlap_R2 = "@R2_perfect\nTACCTTGGGTCCTATGGTATCGGTACTGGTACTTGCTTAATGTCAAGT\n+\n111111111111111111111111111111111111111111111111\n";
        
        const std::string readData_normal_overlap_R1 = "@R1_\nGGTAAACCATTAACCAAATTGGACATATCTCATCTATCCTTATACTTA\n+\n111111111111111111111111111111111111111111111111\n";
        const std::string readData_normal_overlap_R1_short_R1 = "@R1_\nGGTAAACCATTAACCAAATTGGACATATCTCATCTATCCTTATA\n+\n11111111111111111111111111111111111111111111\n";
        const std::string readData_normal_overlap_R2 = "@R2_\nCTCTCGGTCTCCTCTCGTTTCTCGTTCGCGCTAAGTATAAGGATAGA\n+\n11111111111111111111111111111111111111111111111\n";
        
        const std::string readData_trim_overlap_R1 = "@R1_\nACCATAACATAAACCATTAACCAAATTGGACATATCTCATCTATCCTTATACTTA\n+\n1111111111111111111111111111111111111111111111111111111\n";
        const std::string readData_trim_overlap_R1_short_R1 = "@R1_\nACCATAACATAAACCATTAACCAAATTGGACATATCTCATCTATCCTTATAC\n+\n1111111111111111111111111111111111111111111111111111\n";
        const std::string readData_trim_overlap_R2 = "@R2_\nGGTTTATGTTATGGTAATATAGTATAGAGTATAGTTGCGTC\n+\n11111111111111111111111111111111111111111\n";

        const std::string readData_engulf_r2_R1 = "@R1_\nGGTAAACCATTAACCAAATTGGACATATCTCATCTATCCTTATACTTA\n+\n111111111111111111111111111111111111111111111111\n";
        const std::string readData_engulf_r2_R2 = "@R2_\nAGATGAGATATGTCCAATTTGGTTAATGGTT\n+\n1111111111111111111111111111111";
        
        const std::string readData_engulf_r1_R1 = "@R1_\nAACCATTAACCAAATTGGACATATCTCATCT\n+\n1111111111111111111111111111111\n";
        const std::string readData_engulf_r1_R2 = "@R2_\nTAAGTATAAGGATAGATGAGATATGTCCAATTTGGTTAATGGTTTACC\n+\n111111111111111111111111111111111111111111111111\n";

        //error read
        const std::string readData_test_R1 = "@R1_\nCGCCCAGCTAGTTCTTCAAATCTTTAACAAATACTTCAGCAACTCACCTTCTGGGTCTCTTTCCCCCCTATCTGTTGCAATAAACACAGGGGTGTCCAAG\n+\nAAFFFJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ\n";
        const std::string readData_test_R2 = "@R2_\nAAAAAAAAAAAA\n+\n<A-7-<7FJF<<";
    Overlapper ol;
};

TEST_F(OverlapperTest, test_check_length) {
    std::istringstream in1(readData_test_R1);
    std::istringstream in2(readData_test_R2);

    InputReader<PairedEndRead, PairedEndReadFastqImpl> ifp(in1, in2);
    SingleEndReadPtr rb;
        
    while(ifp.has_next()) {
        auto i = ifp.next();
        PairedEndReadPtr per = std::move(i);
        rb = ol.check_read(per, misDensity, 100, 8, 20, kmer, kmerOffset).second;
        ASSERT_EQ(rb, nullptr);
        const std::string s1 =  "CGCCCAGCTAGTTCTTCAAATCTTTAACAAATACTTCAGCAACTCACCTTCTGGGTCTCTTTCCCCCCTATCTGTTGCAATAAACACAGGGGTGTCCAAG";
        const std::string s2 =  "AAAAAAAAAAAA";
        ASSERT_EQ( (per->non_const_read_one()).get_sub_seq(), s1);
        ASSERT_EQ( (per->non_const_read_two()).get_sub_seq(), s2);
    }
};


TEST_F(OverlapperTest, engulfR1) {
    std::istringstream in1(readData_engulf_r1_R1);
    std::istringstream in2(readData_engulf_r1_R2);

    InputReader<PairedEndRead, PairedEndReadFastqImpl> ifp(in1, in2);
    SingleEndReadPtr rb;
        
    while(ifp.has_next()) {
        auto i = ifp.next();
        PairedEndReadPtr per = std::move(i);
        rb = ol.check_read(per, misDensity, 100, 10, 10, kmer, kmerOffset).second;
    }
    const Read &r = rb->get_read();
    ASSERT_EQ(r.get_seq(), "AACCATTAACCAAATTGGACATATCTCATCTATCCTTATACTTA");
};

TEST_F(OverlapperTest, engulfR2) {
    std::istringstream in1(readData_engulf_r2_R1);
    std::istringstream in2(readData_engulf_r2_R2);

    InputReader<PairedEndRead, PairedEndReadFastqImpl> ifp(in1, in2);
    SingleEndReadPtr rb;
    
    while(ifp.has_next()) {
        auto i = ifp.next();
        PairedEndReadPtr per = std::move(i);
        rb = ol.check_read(per, misDensity, 100, 10, 10, kmer, kmerOffset).second;
    }
    const Read &r = rb->get_read();
    ASSERT_EQ(r.get_seq(), "GGTAAACCATTAACCAAATTGGACATATCTCATCT");
};


TEST_F(OverlapperTest, trim) {
    std::istringstream in1(readData_trim_overlap_R1);
    std::istringstream in2(readData_trim_overlap_R2);

    InputReader<PairedEndRead, PairedEndReadFastqImpl> ifp(in1, in2);
    SingleEndReadPtr rb;
    
    while(ifp.has_next()) {
        auto i = ifp.next();
        PairedEndReadPtr per = std::move(i);
        rb = ol.check_read(per, misDensity, 100, 10, 10, kmer, kmerOffset).second;
    }
    const Read &r = rb->get_read();
    ASSERT_EQ(r.get_seq(), "ACCATAACATAAACC");
};

TEST_F(OverlapperTest, trim_short_R1) {
    std::istringstream in1(readData_trim_overlap_R1_short_R1);
    std::istringstream in2(readData_trim_overlap_R2);

    InputReader<PairedEndRead, PairedEndReadFastqImpl> ifp(in1, in2);
    SingleEndReadPtr rb;
    
    while(ifp.has_next()) {
        auto i = ifp.next();
        PairedEndReadPtr per = std::move(i);
        rb = ol.check_read(per, misDensity, 100, 10, 10, kmer, kmerOffset).second;
    }
    const Read &r = rb->get_read();
    ASSERT_EQ(r.get_seq(), "ACCATAACATAAACC");
};

TEST_F(OverlapperTest, normal) {
    std::istringstream in1(readData_normal_overlap_R1);
    std::istringstream in2(readData_normal_overlap_R2);

    InputReader<PairedEndRead, PairedEndReadFastqImpl> ifp(in1, in2);
    SingleEndReadPtr rb;
    
    while(ifp.has_next()) {
        auto i = ifp.next();
        PairedEndReadPtr per = std::move(i);
        rb = ol.check_read(per, misDensity, 100, 10, 10, kmer, kmerOffset).second;
    }
    const Read &r = rb->get_read();
    ASSERT_EQ(r.get_seq(), "GGTAAACCATTAACCAAATTGGACATATCTCATCTATCCTTATACTTAGCGCGAACGAGAAACGAGAGGAGACCGAGAG");
};

TEST_F(OverlapperTest, normal_short_R1) {
    std::istringstream in1(readData_normal_overlap_R1_short_R1);
    std::istringstream in2(readData_normal_overlap_R2);

    InputReader<PairedEndRead, PairedEndReadFastqImpl> ifp(in1, in2);
    SingleEndReadPtr rb;
    
    while(ifp.has_next()) {
        auto i = ifp.next();
        PairedEndReadPtr per = std::move(i);
        rb = ol.check_read(per, misDensity, 100, 10, 10, kmer, kmerOffset).second;
    }
    const Read &r = rb->get_read();
    ASSERT_EQ(r.get_seq(), "GGTAAACCATTAACCAAATTGGACATATCTCATCTATCCTTATACTTAGCGCGAACGAGAAACGAGAGGAGACCGAGAG");
};

TEST_F(OverlapperTest, perfectOverlap) {
    std::istringstream in1(readData_perfect_overlap_R1);
    std::istringstream in2(readData_perfect_overlap_R2);

    InputReader<PairedEndRead, PairedEndReadFastqImpl> ifp(in1, in2);
    SingleEndReadPtr rb;
    
    while(ifp.has_next()) {
        auto i = ifp.next();
        PairedEndReadPtr per = std::move(i);
        rb = ol.check_read(per, misDensity, 100, 10, 10, kmer, kmerOffset).second;
    }
    const Read &r = rb->get_read();
    ASSERT_EQ(r.get_seq(), "ACTTGACATTAAGCAAGTACCAGTACCGATACCATAGGACCCAAGGTA");
};

TEST_F(OverlapperTest, perfectOverlap_short_R1) {
    std::istringstream in1(readData_perfect_overlap_R1_short_R1);
    std::istringstream in2(readData_perfect_overlap_R2);

    InputReader<PairedEndRead, PairedEndReadFastqImpl> ifp(in1, in2);
    SingleEndReadPtr rb;
    
    while(ifp.has_next()) {
        auto i = ifp.next();
        PairedEndReadPtr per = std::move(i);
        rb = ol.check_read(per, misDensity, 100, 10, 10, kmer, kmerOffset).second;
    }
    const Read &r = rb->get_read();
    ASSERT_EQ(r.get_seq(), "ACTTGACATTAAGCAAGTACCAGTACCGATACCATAGGACCCAAGGTA");
};

TEST_F(OverlapperTest, noOverlap) {
    std::istringstream in1(readData_normal_overlap_R1);
    std::istringstream in2(readData_perfect_overlap_R2);

    InputReader<PairedEndRead, PairedEndReadFastqImpl> ifp(in1, in2);
    SingleEndReadPtr rb;
    
    while(ifp.has_next()) {
        auto i = ifp.next();
        PairedEndReadPtr per = std::move(i);
        rb = ol.check_read(per, misDensity, 100, 10, 10, kmer, kmerOffset).second;
        ASSERT_EQ(rb, nullptr);
        const std::string s1 =  "GGTAAACCATTAACCAAATTGGACATATCTCATCTATCCTTATACTTA";
        const std::string s2 =  "TACCTTGGGTCCTATGGTATCGGTACTGGTACTTGCTTAATGTCAAGT";
        ASSERT_EQ( (per->non_const_read_one()).get_sub_seq(), s1);
        ASSERT_EQ( (per->non_const_read_two()).get_sub_seq(), s2);
    }
};

TEST_F(OverlapperTest, noOverlap_short_R1) {
    std::istringstream in1(readData_normal_overlap_R1_short_R1);
    std::istringstream in2(readData_perfect_overlap_R2);

    InputReader<PairedEndRead, PairedEndReadFastqImpl> ifp(in1, in2);
    SingleEndReadPtr rb;
    
    while(ifp.has_next()) {
        auto i = ifp.next();
        PairedEndReadPtr per = std::move(i);
        rb = ol.check_read(per, misDensity, 100, 10, 10, kmer, kmerOffset).second;
        ASSERT_EQ(rb, nullptr);
        const std::string s1 =  "GGTAAACCATTAACCAAATTGGACATATCTCATCTATCCTTATA";
        const std::string s2 =  "TACCTTGGGTCCTATGGTATCGGTACTGGTACTTGCTTAATGTCAAGT";
        ASSERT_EQ( (per->non_const_read_one()).get_sub_seq(), s1);
        ASSERT_EQ( (per->non_const_read_two()).get_sub_seq(), s2);
    }
};

