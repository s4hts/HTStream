#include <gtest/gtest.h>
#include <boost/dynamic_bitset.hpp>
#include <boost/optional/optional_io.hpp>
#include "read.h"

TEST(EmptyTest, firstTest) {
    ASSERT_EQ(0, 0.0);
}

TEST(CreatePERead, createPEReadWorks){
    Read r1("CCACCCTCATTTCATTCTCAGAAGCATGTATGAAGTTGTAATAGCCCTGACGTATGGTTTACCTACTAAGATACCCTCAGGAGTTCTCATCTAGCAAGTG",

            "88@BCFDE<CEFFCEFEEEA<99,C,C,C9E9,,@CE9E9<CC,C6@E,,C,B8E8,CEE8CEE9CE9,,,C<CCCED<,,,,CCEEF9EFE9,C<,,C,", "Read1");
    Read r2("CTTTCTGGAACTTGAGCAGGAGTTCTGCTCTGTCATCTCTGTTCTCCTGTTCCTTCCACACCTGTTTTTTTCTCACCGTGCCATTTTTCCCTTCATTCTC",
            "-8A@@###############################################################################################", "Read2");

    std::shared_ptr<ReadBase> pe1 = std::make_shared<PairedEndRead>(r1, r2);
    std::shared_ptr<ReadBase> se1 = std::make_shared<SingleEndRead>(r1);
    ASSERT_EQ(pe1->get_key(2, 5), pe1->str_to_bit("TTCTGACCCT"));
    //Because we multiple key by two
    //Also needs a 'C' so SE don't map to the same deal
    ASSERT_EQ(se1->get_key(2, 5), se1->str_to_bit("CACCCTCATTT"));
    std::cout << se1->get_key(2,5) << '\n'; 

}

TEST(ConvertToTwobit, ConvertToTwobitWorks){
  boost::dynamic_bitset<> x(8);
  boost::dynamic_bitset<> bs1;
  //ACTG = 00 01 11 10
  x[7] = 0; x[6] = 0; x[5] = 0; x[4] = 1; x[3] = 1; x[2] = 1; x[1] = 1; x[0] = 0;
  bs1 = *ReadBase::str_to_bit("ACTG");
  ASSERT_EQ(x, bs1);
}

TEST(qual, avgQualScore) {
    SingleEndRead r1(Read("aaa",
                       "aaa",
                       "foo"));
    ASSERT_EQ(r1.avg_q_score(), 97.0);

}

TEST(read, BitToStrWorks) {
    auto bs1 = *ReadBase::str_to_bit("ACTG");
    auto s = ReadBase::bit_to_str(bs1);
    ASSERT_EQ(s, "ACTG");
}

TEST(read, ReverseComplementWorks) {
    ASSERT_EQ(ReadBase::bit_to_str(*ReadBase::reverse_complement("ACTG", 0, 4)), "CAGT");

    ASSERT_EQ(ReadBase::bit_to_str(*ReadBase::reverse_complement("ACTAACTGTA", 2, 4)), "CAGT");
    
}
