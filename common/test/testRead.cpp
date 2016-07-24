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

  PairedEndRead pe1(r1, r2, "PE1")

}

class BinarySearchTest: public ::testing::Test {
public:
    BinarySearchTest() : bst(4, 2) {}

protected:

    virtual void SetUp() {

        // todo:  make sensible strings here
        auto r1 = std::make_shared<readInfo>("aattttaaa", "aattttaaa", "aattttaaa");
        auto r2 = std::make_shared<readInfo>("aaggggaaa", "aaggggaaa", "aaggggaaa");

        bst.AddNode(r1, r2);

    }

    virtual void TearDown() {
    }

    BinarySearchTree bst;

};




TEST_F(BinarySearchTest, getIDWorks) {
    auto r1 = std::make_shared<readInfo>("aattttaaa", "aattttaaa", "aattttaaa");
    auto r2 = std::make_shared<readInfo>("aaggggaaa", "aaggggaaa", "aaggggaaa");

    BinarySearchTree::idptr id = 0;
    ASSERT_TRUE(bst.getID(r1, r2, id)) << "reads expected to be found, already added";

}
