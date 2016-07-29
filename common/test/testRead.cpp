#include <gtest/gtest.h>
#include <boost/dynamic_bitset.hpp>
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

    ASSERT_EQ(pe1->getKey(2, 5), pe1->strToBit("ACCCTTTCTG"));
    ASSERT_EQ(se1->getKey(2, 5), se1->strToBit("ACCCTCATTT"));

}

TEST(ConvertToTwobit, ConvertToTwobitWorks){
  boost::dynamic_bitset<> x(8);
  boost::dynamic_bitset<> bs1;
  //ACTG = 00 01 10 11
  x[7] = 0; x[6] = 0; x[5] = 0; x[4] = 1; x[3] = 1; x[2] = 0; x[1] = 1; x[0] = 1;
  bs1 = ReadBase::strToBit("ACTG");
  ASSERT_EQ(x, bs1);
}

TEST(qual, avgQualScore) {
    SingleEndRead r1(Read("aaa",
                       "aaa",
                       "foo"));
    ASSERT_EQ(r1.avg_q_score(), 97.0);

}
