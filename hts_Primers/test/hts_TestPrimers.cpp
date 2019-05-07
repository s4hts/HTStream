#include "gtest/gtest.h"
#include <sstream>
#include <iostream>
#include "hts_Primers.h"

class Primer : public ::testing::Test {
    public:
        const std::string primer_1 =        "TTCATTAAAAATTGAATTGACATTAACCT";
        const std::string read_seq_1 = "CGGGTTTCATTAAAAATTGAATTGACATTAACCTATAAAAATAGGCGTCGAGGCCCTTTCGTCTTCTATCGGAGCTCCAAGACCGCCTCGGCGTGAAGGTGGTGATAGCGCCCGGAAGAGAGTCAATTCAGGGTGGTGAATACTCTAGATC";
        const std::string readData_1 = "@Read1\nTTTTTNGAAAAAAAAAGNTTTTT\n+\n#######################\n";
};

TEST_F(Primer, bounded_edit_distance) {
  ALIGNPOS r;

  r = bounded_edit_distance(primer_1, read_seq_1, 6, 4, 4);

  ASSERT_EQ(0, r.dist);
  ASSERT_EQ(5, r.spos);
  ASSERT_EQ(34, r.epos);
};
