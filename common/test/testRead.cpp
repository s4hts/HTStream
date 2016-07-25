#include <gtest/gtest.h>
#include <boost/dynamic_bitset.hpp>
#include "read.h"

TEST(EmptyTest, firstTest) {
    ASSERT_EQ(0, 0.0);
}

TEST(CreatePERead, createPEReadWorks){
    Read r1("CCACCCTCATTTCATTCTCAGAAGCATGTATGAAGTTGTAATAGCCCTGACGTATGGTTTACCTACTAAGATACCCTCAGGAGTTCTCATCTAGCAAGTG",
          "88@BCFDE<CEFFCEFEEEA<99,C,C,C9E9,,@CE9E9<CC,C6@E,,C,B8E8,CEE8CEE9CE9,,,C<CCCED<,,,,CCEEF9EFE9,C<,,C,");
    Read r2("CTTTCTGGAACTTGAGCAGGAGTTCTGCTCTGTCATCTCTGTTCTCCTGTTCCTTCCACACCTGTTTTTTTCTCACCGTGCCATTTTTCCCTTCATTCTC",
          "-8A@@###############################################################################################");

    std::shared_ptr<ReadBase> pe1 = std::make_shared<PairedEndRead>(r1, r2, "PERead");
    std::shared_ptr<ReadBase> se1 = std::make_shared<SingleEndRead>(r1, "SERead1");

    ASSERT_EQ(se1->getId(), "SERead1");
    ASSERT_EQ(pe1->getId(), "PERead");

    ASSERT_EQ(pe1->getStrKey(2, 5), "ACCCTTTCTG");
    ASSERT_EQ(se1->getStrKey(2, 5), "ACCCTCATTT");

  }

TEST(ConvertToTwobit, ConvertToTwobitWorks){
  Read r1("CCACCCTCATTTCATTCTCAGAAGCATGTATGAAGTTGTAATAGCCCTGACGTATGGTTTACCTACTAAGATACCCTCAGGAGTTCTCATCTAGCAAGTG",
        "88@BCFDE<CEFFCEFEEEA<99,C,C,C9E9,,@CE9E9<CC,C6@E,,C,B8E8,CEE8CEE9CE9,,,C<CCCED<,,,,CCEEF9EFE9,C<,,C,");
  std::shared_ptr<ReadBase> se1 = std::make_shared<SingleEndRead>(r1, "SERead1");

  boost::dynamic_bitset<> x(8);
  boost::dynamic_bitset<> bs1;
  //ACTG = 00011011
  x[7] = 0; x[6] = 0; x[5] = 0; x[4] = 1; x[3] = 1; x[2] = 0; x[1] = 1; x[0] = 1;
  bs1 = se1->strToBit("ACTG");
  ASSERT_EQ(x, bs1);
}
