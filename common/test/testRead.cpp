#include <gtest/gtest.h>
#include "read.h"

TEST(EmptyTest, firstTest) {
    ASSERT_EQ(0, 0.0);
}

TEST(CreatePERead, createPEReadWorks){
  Read r1("CCACCCTCATTTCATTCTCAGAAGCATGTATGAAGTTGTAATAGCCCTGACGTATGGTTTACCTACTAAGATACCCTCAGGAGTTCTCATCTAGCAAGTG",
          "88@BCFDE<CEFFCEFEEEA<99,C,C,C9E9,,@CE9E9<CC,C6@E,,C,B8E8,CEE8CEE9CE9,,,C<CCCED<,,,,CCEEF9EFE9,C<,,C,");
  Read r2("CTTTCTGGAACTTGAGCAGGAGTTCTGCTCTGTCATCTCTGTTCTCCTGTTCCTTCCACACCTGTTTTTTTCTCACCGTGCCATTTTTCCCTTCATTCTC",
          "-8A@@###############################################################################################");

  PairedEndRead pe1(r1, r2, "PERead");
    ASSERT_EQ(0, 0.0);
}
