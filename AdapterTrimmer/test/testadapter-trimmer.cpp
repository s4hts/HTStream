#include "gtest/gtest.h"
#include <sstream>
#include <iostream>
#include "adapter-trimmer.h"

class Adapter : public ::testing::Test {
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

TEST_F(Adapter, engulfR1) {
    std::istringstream in1(readData_engulf_r1_R1);
    std::istringstream in2(readData_engulf_r1_R2);

    InputReader<PairedEndRead, PairedEndReadFastqImpl> ifp(in1, in2);
    unsigned int rb;
	

    while(ifp.has_next()) {
        auto i = ifp.next();
        PairedEndRead *per = dynamic_cast<PairedEndRead*>(i.get());
        rb = check_read(*per, misDensity, 10, 10, kmer, kmerOffset, 10);
        const std::string s1 =  "TAAGTATAAGGATAGATGAGATATGTCCAATTTGGTTAATGGTTTACC";
        const std::string s2 =  "AACCATTAACCAAATTGGACATATCTCATCT";
        ASSERT_EQ( (per->non_const_read_one()).get_sub_seq(), s1);
        ASSERT_EQ( (per->non_const_read_two()).get_sub_seq(), s2);
    }
};

TEST_F(Adapter, engulfR2) {
    std::istringstream in1(readData_engulf_r2_R1);
    std::istringstream in2(readData_engulf_r2_R2);

    InputReader<PairedEndRead, PairedEndReadFastqImpl> ifp(in1, in2);
    unsigned int rb;

    while(ifp.has_next()) {
        auto i = ifp.next();
        PairedEndRead *per = dynamic_cast<PairedEndRead*>(i.get());
        rb = check_read(*per, misDensity, 10, 10, kmer, kmerOffset, 10);
        const std::string s1 =  "GGTAAACCATTAACCAAATTGGACATATCTCATCTATCCTTATACTTA";
        const std::string s2 =  "AGATGAGATATGTCCAATTTGGTTAATGGTT";
        ASSERT_EQ( (per->non_const_read_one()).get_sub_seq(), s1);
        ASSERT_EQ( (per->non_const_read_two()).get_sub_seq(), s2);

    }
};


TEST_F(Adapter, trim) {
    std::istringstream in1(readData_trim_overlap_R1);
    std::istringstream in2(readData_trim_overlap_R2);

    InputReader<PairedEndRead, PairedEndReadFastqImpl> ifp(in1, in2);
    unsigned int rb;

    while(ifp.has_next()) {
        auto i = ifp.next();
        PairedEndRead *per = dynamic_cast<PairedEndRead*>(i.get());
        rb = check_read(*per, misDensity, 10, 10, kmer, kmerOffset, 10);
        ASSERT_EQ( (per->non_const_read_one()).get_sub_seq(), "ACCATAACATAAACCATTAACCAAATTGGACATATCTCATCTATCCTTATACTTA");
        ASSERT_EQ( (per->non_const_read_two()).get_sub_seq(), "GGTTTATGTTATGGT");

    }
    
};


TEST_F(Adapter, normal) {
    std::istringstream in1(readData_normal_overlap_R1);
    std::istringstream in2(readData_normal_overlap_R2);

    InputReader<PairedEndRead, PairedEndReadFastqImpl> ifp(in1, in2);
    unsigned int rb;

    while(ifp.has_next()) {
        auto i = ifp.next();
        PairedEndRead *per = dynamic_cast<PairedEndRead*>(i.get());
        rb = check_read(*per, misDensity, 10, 10, kmer, kmerOffset, 10);
        ASSERT_EQ( (per->non_const_read_one()).get_sub_seq(), "GGTAAACCATTAACCAAATTGGACATATCTCATCTATCCTTATACTTA");
        ASSERT_EQ( (per->non_const_read_two()).get_sub_seq(), "CTCTCGGTCTCCTCTCGTTTCTCGTTCGCGCTAAGTATAAGGATAGA");
    }
};


TEST_F(Adapter, perfectOverlap) {
    std::istringstream in1(readData_perfect_overlap_R1);
    std::istringstream in2(readData_perfect_overlap_R2);

    InputReader<PairedEndRead, PairedEndReadFastqImpl> ifp(in1, in2);
    unsigned int rb;


    while(ifp.has_next()) {
        auto i = ifp.next();
        PairedEndRead *per = dynamic_cast<PairedEndRead*>(i.get());
        rb = check_read(*per, misDensity, 10, 10, kmer, kmerOffset, 10);
        ASSERT_EQ( (per->non_const_read_one()).get_sub_seq(), "ACTTGACATTAAGCAAGTACCAGTACCGATACCATAGGACCCAAGGTA");
        ASSERT_EQ( (per->non_const_read_two()).get_sub_seq(), "TACCTTGGGTCCTATGGTATCGGTACTGGTACTTGCTTAATGTCAAGT");
    }
};

