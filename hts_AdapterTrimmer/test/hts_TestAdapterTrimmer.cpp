#include "gtest/gtest.h"
#include <sstream>
#include <iostream>
#include "hts_AdapterTrimmer.h"

class Adapter : public ::testing::Test {
    public:
        const double misDensity = 0.25;
        const size_t kmer = 8;
        const size_t kmerOffset = 1;
        const std::string readData_perfect_overlap_R1 = "@R1_perfect\nACTTGACATTAAGCAAGTACCAGTACCGATACCATAGGACCCAAGGTA\n+\n111111111111111111111111111111111111111111111111\n";
        const std::string readData_perfect_overlap_R1_short_R1 = "@R1_perfect\nACTTGACATTAAGCAAGTACCAGTACCGATACCATAGGACCCAA\n+\n11111111111111111111111111111111111111111111\n";
        const std::string readData_perfect_overlap_R2 = "@R2_perfect\nTACCTTGGGTCCTATGGTATCGGTACTGGTACTTGCTTAATGTCAAGT\n+\n111111111111111111111111111111111111111111111111\n";

        const std::string readData_normal_overlap_R1 = "@R1_\nGGTAAACCATTAACCAAATTGGACATATCTCATCTATCCTTATACTTA\n+\n111111111111111111111111111111111111111111111111\n";
        const std::string readData_normal_overlap_R1_short_R1 = "@R1_\nGGTAAACCATTAACCAAATTGGACATATCTCATCTATCCTTATA\n+\n11111111111111111111111111111111111111111111\n";
        const std::string readData_normal_overlap_R2 = "@R2_\nCTCTCGGTCTCCTCTCGTTTCTCGTTCGCGCTAAGTATAAGGATAGA\n+\n11111111111111111111111111111111111111111111111\n";

        const std::string readData_trim_overlap_R1 = "@R1_\nACCATAACATAAACCATTAACCAAATTGGACATATCTCATCTATCCTTATACTTA\n+\n1111111111111111111111111111111111111111111111111111111\n";
        const std::string readData_trim_overlap_R1_short_R1 = "@R1_\nACCATAACATAAACCATTAACCAAATTGGACATATCTCATCTATCCTTATA\n+\n111111111111111111111111111111111111111111111111111\n";
        const std::string readData_trim_overlap_R2 = "@R2_\nGGTTTATGTTATGGTAATATAGTATAGAGTATAGTTGCGTC\n+\n11111111111111111111111111111111111111111\n";

        const std::string readData_engulf_r2_R1 = "@R1_\nGGTAAACCATTAACCAAATTGGACATATCTCATCTATCCTTATACTTA\n+\n111111111111111111111111111111111111111111111111\n";
        const std::string readData_engulf_r2_R2 = "@R2_\nAGATGAGATATGTCCAATTTGGTTAATGGTT\n+\n1111111111111111111111111111111";

        const std::string readData_engulf_r1_R1 = "@R1_\nCATTAACCAAATTGGACATATCTCATCT\n+\n1111111111111111111111111111\n";
        const std::string readData_engulf_r1_R2 = "@R2_\nTAAGTATAAGGATAGATGAGATATGTCCAATTTGGTTAATGGTTTACC\n+\n111111111111111111111111111111111111111111111111\n";

        // @J00113:301:HFYJ7BBXX:6:1101:26920:4462 1:N:0:GCCAAGAC+GAAGGAAG
        const std::string failed_R1 = "@R1_\nCTCTGTCTCAAAAAAAAAAAAAAAAAAAGTTCAATGAAAGTAGCAGCCAATACCTCCCACACCCACCCCCACCCCGCCCTGCACCTTGCCAGGGACCTTA\n+\nA<A<FAJFFJJJJFFFJJJJJJJJJJJF<7A--<A------7----77----<--<-77A-A-F<A77FJ<7<J-77FFAA7<A-AF-7--77--77---\n";
        const std::string failed_R2 = "@R2_\nAAGGACCCTGGCAATGTGCAGGGCGGGGTGGGGGTGTGTGTGAGATGTATTAGCTGCTACTTTCATCTAACTTTTTTTTTTTTTTTTTTTTAGACAGAAT\n+\nAAFFFJJJJFJJJJJJJJJJJJJAJJJJFJJJJJAJJJJJJJJJJJJJJJJJJJJJJJJJJFJJJJJJJJJJJJJJJFJJJJJJJJJJJJA7--------\n";

        const std::string issue_126_read_1 = "@R1\nNTTTGTAGGGCAGAAAAATTCTGCCTGACCTTTCCCATCTTCCTGATGATGTTAACAAACGACAGCCACTGAGGGCACAGGAGGGGGCTCTACACAGGCCAAGAGCTACAGGCTACAGTCTGAGCACAGAAAACCCTCCAGCTCCAGCAG\n+\n#AAAAJJJJJJJJJJJJJJJJJJJJFJFJJJJJJJJJJJFJJJFJAJFJJFFJFJJJJJJJJJJFJJ-FJJJJJFJJJJJ7JFJFJJAA<JJJFJFFJJJJJJJJJJJJJJJJJJJJJFAAF-F7FFAJJJJJFJJFFJJFJ7FJJFFJJ\n";
        const std::string issue_126_read_2 = "@R2\nNCTGGAGCTGGAGGGTTTTCTGTGCTCAGACTGTAGCCTGTAGCTCTTGGCCTGTGTAGAGCCCCCTCCTGTGCCCTCAGTGGCTGTCGTTTGTTAACATCATCAGGAAGATGGGAAAGGTCAGGCAGAATTTTTCTGCCCTCCAAAGCG\n+\n#AAFFJJFJJJFJJJAAJJJJJJJJFJFJJJJJJJJJJJJJJJJ7-FAJJJJJJJJFAFFFFJFJJJJJJJJ7FJJ7JFJFJJJFJFJAAAJFJJJJFFFJJ<J<FJ-JJ--7A7<JAFJJJ<FJJJJJ<AF<AAF<-<)<))77-<--7\n";


};

TEST_F(Adapter, issue_126) {
    std::istringstream in1(issue_126_read_1);
    std::istringstream in2(issue_126_read_2);

    InputReader<PairedEndRead, PairedEndReadFastqImpl> ifp(in1, in2);
    while(ifp.has_next()) {
        auto i = ifp.next();
        std::shared_ptr<PairedEndRead> per = std::make_shared<PairedEndRead>(std::move(*dynamic_cast<PairedEndRead*>(i.get())));
        check_read_pe(per, misDensity, 100, 8, 20, kmer, kmerOffset, false);
        const std::string s1 =  "CTTTGTAGGGCAGAAAAATTCTGCCTGACCTTTCCCATCTTCCTGATGATGTTAACAAACGACAGCCACTGAGGGCACAGGAGGGGGCTCTACACAGGCCAAGAGCTACAGGCTACAGTCTGAGCACAGAAAACCCTCCAGCTCCAGC";
        const std::string s2 =  "GCTGGAGCTGGAGGGTTTTCTGTGCTCAGACTGTAGCCTGTAGCTCTTGGCCTGTGTAGAGCCCCCTCCTGTGCCCTCAGTGGCTGTCGTTTGTTAACATCATCAGGAAGATGGGAAAGGTCAGGCAGAATTTTTCTGCCCTACAAAG";
        ASSERT_EQ( (per->non_const_read_one()).get_sub_seq(), s1);
        ASSERT_EQ( (per->non_const_read_two()).get_sub_seq(), s2);
    }
}

TEST_F(Adapter, failed_read) {
    std::istringstream in1(failed_R1);
    std::istringstream in2(failed_R2);

    InputReader<PairedEndRead, PairedEndReadFastqImpl> ifp(in1, in2);

    while(ifp.has_next()) {
        auto i = ifp.next();
        std::shared_ptr<PairedEndRead> per = std::make_shared<PairedEndRead>(*dynamic_cast<PairedEndRead*>(i.get()));
        ASSERT_EQ( (per->non_const_read_two()).getRTrim(), 0u);
        check_read_pe(per, misDensity, 100, 8, 20, kmer, kmerOffset, true);
        const std::string s1 =  "CTCTGTCTCAAAAAAAAAAAAAAAAAAAGTTCAATGAAAGTAGCAGCCAATACCTCCCACACCCACCCCCACCCCGCCCTGCACCTTGCCAGGGACCTT";
        const std::string s2 =  "AAGGACCCTGGCAATGTGCAGGGCGGGGTGGGGGTGTGTGTGAGATGTATTAGCTGCTACTTTCATCTAACTTTTTTTTTTTTTTTTTTTTAGACAGAA";
        ASSERT_EQ( (per->non_const_read_two()).getRTrim(), 1u);
        ASSERT_EQ( (per->non_const_read_one()).get_sub_seq(), s1);
        ASSERT_EQ( (per->non_const_read_two()).get_sub_seq(), s2);
    }
};

TEST_F(Adapter, engulfR1) {
    std::istringstream in1(readData_engulf_r1_R1);
    std::istringstream in2(readData_engulf_r1_R2);

    InputReader<PairedEndRead, PairedEndReadFastqImpl> ifp(in1, in2);

    while(ifp.has_next()) {
        auto i = ifp.next();
        std::shared_ptr<PairedEndRead> per = std::make_shared<PairedEndRead>(*dynamic_cast<PairedEndRead*>(i.get()));
        check_read_pe(per, misDensity, 100, 10, 10, kmer, kmerOffset, 10);
        const std::string s1 =  "CATTAACCAAATTGGACATATCTCATCT";
        const std::string s2 =  "TAAGTATAAGGATAGATGAGATATGTCCAATTTGGTTAATG";
        ASSERT_EQ( (per->non_const_read_one()).get_sub_seq(), s1);
        ASSERT_EQ( (per->non_const_read_two()).get_sub_seq(), s2);
    }
};

TEST_F(Adapter, engulfR2) {
    std::istringstream in1(readData_engulf_r2_R1);
    std::istringstream in2(readData_engulf_r2_R2);

    InputReader<PairedEndRead, PairedEndReadFastqImpl> ifp(in1, in2);

    while(ifp.has_next()) {
        auto i = ifp.next();
        std::shared_ptr<PairedEndRead> per = std::make_shared<PairedEndRead>(*dynamic_cast<PairedEndRead*>(i.get()));
        check_read_pe(per, misDensity, 100, 10, 10, kmer, kmerOffset, 10);
        const std::string s1 =  "GGTAAACCATTAACCAAATTGGACATATCTCATCT";
        const std::string s2 =  "AGATGAGATATGTCCAATTTGGTTAATGGTT";
        ASSERT_EQ( (per->non_const_read_one()).get_sub_seq(), s1);
        ASSERT_EQ( (per->non_const_read_two()).get_sub_seq(), s2);

    }
};

TEST_F(Adapter, trim) {
    std::istringstream in1(readData_trim_overlap_R1);
    std::istringstream in2(readData_trim_overlap_R2);

    InputReader<PairedEndRead, PairedEndReadFastqImpl> ifp(in1, in2);

    while(ifp.has_next()) {
        auto i = ifp.next();
        std::shared_ptr<PairedEndRead> per = std::make_shared<PairedEndRead>(*dynamic_cast<PairedEndRead*>(i.get()));
        check_read_pe(per, misDensity, 100, 10, 10, kmer, kmerOffset, 10);
        ASSERT_EQ( (per->non_const_read_one()).get_sub_seq(), "ACCATAACATAAACC");
        ASSERT_EQ( (per->non_const_read_two()).get_sub_seq(), "GGTTTATGTTATGGT");
    }
};

TEST_F(Adapter, trim_short_R1) {
    std::istringstream in1(readData_trim_overlap_R1_short_R1);
    std::istringstream in2(readData_trim_overlap_R2);

    InputReader<PairedEndRead, PairedEndReadFastqImpl> ifp(in1, in2);

    while(ifp.has_next()) {
        auto i = ifp.next();
        std::shared_ptr<PairedEndRead> per = std::make_shared<PairedEndRead>(*dynamic_cast<PairedEndRead*>(i.get()));
        check_read_pe(per, misDensity, 100, 10, 10, kmer, kmerOffset, 10);
        ASSERT_EQ( (per->non_const_read_one()).get_sub_seq(), "ACCATAACATAAACC");
        ASSERT_EQ( (per->non_const_read_two()).get_sub_seq(), "GGTTTATGTTATGGT");
    }
};

TEST_F(Adapter, normal) {
    std::istringstream in1(readData_normal_overlap_R1);
    std::istringstream in2(readData_normal_overlap_R2);

    InputReader<PairedEndRead, PairedEndReadFastqImpl> ifp(in1, in2);

    while(ifp.has_next()) {
        auto i = ifp.next();
        std::shared_ptr<PairedEndRead> per = std::make_shared<PairedEndRead>(*dynamic_cast<PairedEndRead*>(i.get()));
        check_read_pe(per, misDensity, 100, 10, 10, kmer, kmerOffset, 10);
        ASSERT_EQ( (per->non_const_read_one()).get_sub_seq(), "GGTAAACCATTAACCAAATTGGACATATCTCATCTATCCTTATACTTA");
        ASSERT_EQ( (per->non_const_read_two()).get_sub_seq(), "CTCTCGGTCTCCTCTCGTTTCTCGTTCGCGCTAAGTATAAGGATAGA");
    }
};

TEST_F(Adapter, normal_short_R1) {
    std::istringstream in1(readData_normal_overlap_R1_short_R1);
    std::istringstream in2(readData_normal_overlap_R2);

    InputReader<PairedEndRead, PairedEndReadFastqImpl> ifp(in1, in2);

    while(ifp.has_next()) {
        auto i = ifp.next();
        std::shared_ptr<PairedEndRead> per = std::make_shared<PairedEndRead>(*dynamic_cast<PairedEndRead*>(i.get()));
        check_read_pe(per, misDensity, 100, 10, 10, kmer, kmerOffset, 10);
        ASSERT_EQ( (per->non_const_read_one()).get_sub_seq(), "GGTAAACCATTAACCAAATTGGACATATCTCATCTATCCTTATA");
        ASSERT_EQ( (per->non_const_read_two()).get_sub_seq(), "CTCTCGGTCTCCTCTCGTTTCTCGTTCGCGCTAAGTATAAGGATAGA");
    }
};

TEST_F(Adapter, perfectOverlap) {
    std::istringstream in1(readData_perfect_overlap_R1);
    std::istringstream in2(readData_perfect_overlap_R2);

    InputReader<PairedEndRead, PairedEndReadFastqImpl> ifp(in1, in2);

    while(ifp.has_next()) {
        auto i = ifp.next();
        std::shared_ptr<PairedEndRead> per = std::make_shared<PairedEndRead>(*dynamic_cast<PairedEndRead*>(i.get()));
        check_read_pe(per, misDensity, 100, 10, 10, kmer, kmerOffset, 10);
        ASSERT_EQ( (per->non_const_read_one()).get_sub_seq(), "ACTTGACATTAAGCAAGTACCAGTACCGATACCATAGGACCCAAGGTA");
        ASSERT_EQ( (per->non_const_read_two()).get_sub_seq(), "TACCTTGGGTCCTATGGTATCGGTACTGGTACTTGCTTAATGTCAAGT");
    }
};

TEST_F(Adapter, perfectOverlap_short_R1) {
    std::istringstream in1(readData_perfect_overlap_R1_short_R1);
    std::istringstream in2(readData_perfect_overlap_R2);

    InputReader<PairedEndRead, PairedEndReadFastqImpl> ifp(in1, in2);

    while(ifp.has_next()) {
        auto i = ifp.next();
        std::shared_ptr<PairedEndRead> per = std::make_shared<PairedEndRead>(*dynamic_cast<PairedEndRead*>(i.get()));
        check_read_pe(per, misDensity, 100, 10, 10, kmer, kmerOffset, 10);
        ASSERT_EQ( (per->non_const_read_one()).get_sub_seq(), "ACTTGACATTAAGCAAGTACCAGTACCGATACCATAGGACCCAA");
        ASSERT_EQ( (per->non_const_read_two()).get_sub_seq(), "TACCTTGGGTCCTATGGTATCGGTACTGGTACTTGCTTAATGTCAAGT");
    }
};

TEST_F(Adapter, noOverlap) {
    std::istringstream in1(readData_normal_overlap_R1);
    std::istringstream in2(readData_perfect_overlap_R2);

    InputReader<PairedEndRead, PairedEndReadFastqImpl> ifp(in1, in2);

    while(ifp.has_next()) {
        auto i = ifp.next();
        std::shared_ptr<PairedEndRead> per = std::make_shared<PairedEndRead>(*dynamic_cast<PairedEndRead*>(i.get()));
        check_read_pe(per, misDensity, 100, 10, 10, kmer, kmerOffset, 10);
        const std::string s1 =  "GGTAAACCATTAACCAAATTGGACATATCTCATCTATCCTTATACTTA";
        const std::string s2 =  "TACCTTGGGTCCTATGGTATCGGTACTGGTACTTGCTTAATGTCAAGT";
        ASSERT_EQ( (per->non_const_read_one()).get_sub_seq(), s1);
        ASSERT_EQ( (per->non_const_read_two()).get_sub_seq(), s2);
    }

};

TEST_F(Adapter, noOverlap_short_R1) {
    std::istringstream in1(readData_normal_overlap_R1_short_R1);
    std::istringstream in2(readData_perfect_overlap_R2);

    InputReader<PairedEndRead, PairedEndReadFastqImpl> ifp(in1, in2);

    while(ifp.has_next()) {
        auto i = ifp.next();
        std::shared_ptr<PairedEndRead> per = std::make_shared<PairedEndRead>(*dynamic_cast<PairedEndRead*>(i.get()));
        check_read_pe(per, misDensity, 100, 10, 10, kmer, kmerOffset, 10);
        const std::string s1 =  "GGTAAACCATTAACCAAATTGGACATATCTCATCTATCCTTATA";
        const std::string s2 =  "TACCTTGGGTCCTATGGTATCGGTACTGGTACTTGCTTAATGTCAAGT";
        ASSERT_EQ( (per->non_const_read_one()).get_sub_seq(), s1);
        ASSERT_EQ( (per->non_const_read_two()).get_sub_seq(), s2);
    }

};

TEST_F(Adapter, se_noAdapterR1) {
    std::istringstream in1(readData_perfect_overlap_R1);

    InputReader<SingleEndRead, SingleEndReadFastqImpl> ifs(in1);

    while(ifs.has_next()) {
        auto i = ifs.next();
        SingleEndReadPtr ser = std::make_shared<SingleEndRead>(*dynamic_cast<SingleEndRead*>(i.get()));
        check_read_se(ser, misDensity, 100, 8, 20, kmer, kmerOffset, "THEREISNOMATCH");
        const std::string s1 =  "ACTTGACATTAAGCAAGTACCAGTACCGATACCATAGGACCCAAGGTA";
        ASSERT_EQ( (ser->non_const_read_one()).get_sub_seq(), s1);
    }
};

TEST_F(Adapter, se_AdapaterR1) {
    std::istringstream in1(readData_perfect_overlap_R1);

    InputReader<SingleEndRead, SingleEndReadFastqImpl> ifs(in1);

    while(ifs.has_next()) {
        auto i = ifs.next();
        SingleEndReadPtr ser = std::make_shared<SingleEndRead>(*dynamic_cast<SingleEndRead*>(i.get()));
        check_read_se(ser, misDensity, 100, 8, 20, kmer, kmerOffset, "CCAGTACCGATACCATA");
        const std::string s1 =  "ACTTGACATTAAGCAAGTA";
        ASSERT_EQ( (ser->non_const_read_one()).get_sub_seq(), s1);
    }
};

TEST_F(Adapter, se_errorAdapterR1) {
    std::istringstream in1(readData_perfect_overlap_R1);

    InputReader<SingleEndRead, SingleEndReadFastqImpl> ifs(in1);

    while(ifs.has_next()) {
        auto i = ifs.next();
        SingleEndReadPtr ser = std::make_shared<SingleEndRead>(*dynamic_cast<SingleEndRead*>(i.get()));
        check_read_se(ser, misDensity, 100, 8, 20, kmer, kmerOffset, "ATAGTACCGATACCATA");
        const std::string s1 =  "ACTTGACATTAAGCAAGTA";
        ASSERT_EQ( (ser->non_const_read_one()).get_sub_seq(), s1);
    }
};
