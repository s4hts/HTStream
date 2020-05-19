#include "gtest/gtest.h"
#include <sstream>
#include <iostream>
#include "hts_PolyATTrim.h"

class PolyATTail : public ::testing::Test {
public:
    const std::string readData_1 = "@Read1\nTTTTTTGGAAAAAAAAAGTTTTTTTTG\n+\n###########################\n";
    const std::string readData_2 = "@Read1\nAAACAAAAAAGGAAAAAAATAAA\n+\n#######################\n";
    const std::string readData_3 = "@Read1\nAAAAAAAAAAAAAAAAAAAAAAAA\n+\n########################\n";
    const std::string readData_4 = "@Read1\nGNTTTTTTCATTGGATGCATTAATAACCCATGTTTTACCTTTTGAAAAAATAAATGAAGGATTTGACCTGCTTCACTCTGGGAAAAGGTAGATTTTTTAG\n+\nA#AFFJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ\n";
    size_t min_trim = 5;
    size_t max_trim = 30;
    size_t window_size = 6;
    size_t perfect_windows = 1;
    size_t min_length = 5;
    double max_mismatch_errorDensity = 0.3;
    PolyATTrim pa;
};

TEST_F(PolyATTail, BasicTrim) {
    std::istringstream in1(readData_1);
    std::istringstream in2(readData_2);

    InputReader<PairedEndRead, PairedEndReadFastqImpl> ifp(in1, in2);

    while(ifp.has_next()) {
        auto i = ifp.next();
        PairedEndRead *per = dynamic_cast<PairedEndRead*>(i.get());
        pa.trim_left(per->non_const_read_one(), 'T', min_trim, max_trim, window_size, max_mismatch_errorDensity, perfect_windows);
        pa.trim_right(per->non_const_read_one(), 'T', min_trim, max_trim, window_size, max_mismatch_errorDensity, perfect_windows);
        ASSERT_EQ("GGAAAAAAAAAG", (per->non_const_read_one()).get_sub_seq());
    }
};

TEST_F(PolyATTail, perfectMatch2) {
    std::istringstream in1(readData_1);
    std::istringstream in2(readData_2);

    InputReader<PairedEndRead, PairedEndReadFastqImpl> ifp(in1, in2);

    while(ifp.has_next()) {
        auto i = ifp.next();
        PairedEndRead *per = dynamic_cast<PairedEndRead*>(i.get());
        pa.trim_left(per->non_const_read_two(), 'A', min_trim, max_trim, window_size, max_mismatch_errorDensity, 2);
        pa.trim_right(per->non_const_read_two(), 'A', min_trim, max_trim, window_size, max_mismatch_errorDensity, 2);
        ASSERT_EQ("AAACAAAAAAGG", (per->non_const_read_two()).get_sub_seq());
    }
};

TEST_F(PolyATTail, mostTrim2) {
    std::istringstream in1(readData_1);
    std::istringstream in2(readData_2);

    InputReader<PairedEndRead, PairedEndReadFastqImpl> ifp(in1, in2);

    while(ifp.has_next()) {
        auto i = ifp.next();
        PairedEndRead *per = dynamic_cast<PairedEndRead*>(i.get());
        pa.trim_left(per->non_const_read_two(), 'A', min_trim, max_trim, window_size, max_mismatch_errorDensity, perfect_windows);
        pa.trim_right(per->non_const_read_two(), 'A', min_trim, max_trim, window_size, max_mismatch_errorDensity, perfect_windows);
        ASSERT_EQ("GG", (per->non_const_read_two()).get_sub_seq());
    }
};

TEST_F(PolyATTail, AllTrim) {
    std::istringstream in1(readData_1);
    std::istringstream in2(readData_3);

    InputReader<PairedEndRead, PairedEndReadFastqImpl> ifp(in1, in2);

    while(ifp.has_next()) {
        auto i = ifp.next();
        PairedEndRead *per = dynamic_cast<PairedEndRead*>(i.get());
        pa.trim_left(per->non_const_read_two(), 'A', min_trim, max_trim, window_size, max_mismatch_errorDensity, perfect_windows);
        pa.trim_right(per->non_const_read_two(), 'A', min_trim, max_trim, window_size, max_mismatch_errorDensity, perfect_windows);
        ASSERT_EQ("N", (per->non_const_read_two()).get_sub_seq());
    }
};

TEST_F(PolyATTail, NoTrim) {
    std::istringstream in1(readData_1);
    std::istringstream in2(readData_3);

    InputReader<PairedEndRead, PairedEndReadFastqImpl> ifp(in1, in2);

    while(ifp.has_next()) {
        auto i = ifp.next();
        PairedEndRead *per = dynamic_cast<PairedEndRead*>(i.get());
        pa.trim_left(per->non_const_read_two(), 'T', min_trim, max_trim, window_size, max_mismatch_errorDensity, perfect_windows);
        pa.trim_right(per->non_const_read_two(), 'T', min_trim, max_trim, window_size, max_mismatch_errorDensity, perfect_windows);
        ASSERT_EQ("AAAAAAAAAAAAAAAAAAAAAAAA", (per->non_const_read_two()).get_sub_seq());
    }
};

TEST_F(PolyATTail, closebutno) {
    std::istringstream in1(readData_4);

    InputReader<SingleEndRead, SingleEndReadFastqImpl> ifs(in1);
    std::shared_ptr<std::ostringstream> out1(new std::ostringstream);
    {
        std::shared_ptr<HtsOfstream> hts_of(new HtsOfstream(out1));
        std::shared_ptr<OutputWriter> tab(new ReadBaseOutTab(hts_of));
        WriterHelper writer(tab, tab, true);
        while(ifs.has_next()) {
            auto i = ifs.next();
            SingleEndRead *ser = dynamic_cast<SingleEndRead*>(i.get());
            pa.trim_left(ser->non_const_read_one(), 'T', min_trim, max_trim, window_size, max_mismatch_errorDensity, perfect_windows);
            pa.trim_right(ser->non_const_read_one(), 'T', min_trim, max_trim, window_size, max_mismatch_errorDensity, perfect_windows);
            writer(*ser);
        }
    }
    ASSERT_EQ("Read1\tGNTTTTTTCATTGGATGCATTAATAACCCATGTTTTACCTTTTGAAAAAATAAATGAAGGATTTGACCTGCTTCACTCTGGGAAAAGGTAGATTTTTTAG\tA#AFFJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ\n", out1->str());
};
