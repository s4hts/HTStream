#include "gtest/gtest.h"
#include <sstream>
#include <iostream>
#include "hts_ExtractUMI.h"

class ExtractUMITest : public ::testing::Test {
public:

    const std::string readData_1 = "@Read1\nNAAAAAGACATTAAGCAA\n+\n!!!!!!############\n";
    const std::string readData_2 = "@Read2\nTTTTTTGACATTAAGCAA\n+\n!!!!!!############\n";
    const std::string dragenData_1 = "@A00887_1:1:2:3:4:5:6\nNAAAAAGACATTAAGCAA\n+\n!!!!!!############\n";
    const std::string dragenData_2 = "@A00887_2:1:2:3:4:5:6\nTTTTTTGACATTAAGCAA\n+\n!!!!!!############\n";

    ExtractUMI eu;

    // init struct for testing
    ExtractUMI::UMI umi = {
        "",         // umi R1 sequence
        "",         // qual sequence
        6,          // umi length
        false,      // discard status (init)
        0,          // quality threshold
        0,          // avg quality threshold
        33,         // quality offset
        false,      // homopolymer discard
        false,       // discard N containing UMIs
        false      // add as tag
    };
};

TEST_F(ExtractUMITest, BasicExtract) { // SE extract test
    std::istringstream in1(readData_1);

    InputReader<SingleEndRead, SingleEndReadFastqImpl> ifp(in1);

    auto i = ifp.next();
    SingleEndRead *ser = dynamic_cast<SingleEndRead*>(i.get());
    eu.extract_umi(ser->non_const_read_one(), umi, ':', '-');
    ASSERT_EQ("Read1:NAAAAA", (ser->non_const_read_one()).get_id_first());
};


TEST_F(ExtractUMITest, ReadOneExtract) { // R1 extract test
    std::istringstream in1(readData_1);
    std::istringstream in2(readData_2);

    InputReader<PairedEndRead, PairedEndReadFastqImpl> ifp(in1, in2);
    std::shared_ptr<std::ostringstream> out1(new std::ostringstream);
    std::shared_ptr<HtsOfstream> hts_of(new HtsOfstream(out1));
    std::shared_ptr<OutputWriter> tab(new ReadBaseOutTab(hts_of));

    // reset umi struct
    std::tie(umi.seq, umi.qual) = std::make_tuple("", "");

    WriterHelper writer(tab, tab);
    while (ifp.has_next()) {
        auto i = ifp.next();
        PairedEndRead *per = dynamic_cast<PairedEndRead*>(i.get());
        eu.extract_umi(per->non_const_read_one(), umi, ':', '-');
        eu.extract_umi(per->non_const_read_two(), umi, ':', '-');
        writer(*per);
    }
    ASSERT_EQ("Read1:NAAAAA\tGACATTAAGCAA\t############\tRead2:NAAAAA\tTTTTTTGACATTAAGCAA\t!!!!!!############\n", out1->str());
};


TEST_F(ExtractUMITest, ReadTwoExtract) { // R2 extract test
    std::istringstream in1(readData_1);
    std::istringstream in2(readData_2);

    InputReader<PairedEndRead, PairedEndReadFastqImpl> ifp(in1, in2);
    std::shared_ptr<std::ostringstream> out1(new std::ostringstream);
    std::shared_ptr<HtsOfstream> hts_of(new HtsOfstream(out1));
    std::shared_ptr<OutputWriter> tab(new ReadBaseOutTab(hts_of));

    // reset umi struct
    std::tie(umi.seq, umi.qual) = std::make_tuple("", "");

    WriterHelper writer(tab, tab);
    while (ifp.has_next()) {
        auto i = ifp.next();
        PairedEndRead *per = dynamic_cast<PairedEndRead*>(i.get());
        eu.extract_umi(per->non_const_read_two(), umi, ':', '-');
        eu.extract_umi(per->non_const_read_one(), umi, ':', '-');
        writer(*per);
    }
    ASSERT_EQ("Read1:TTTTTT\tNAAAAAGACATTAAGCAA\t!!!!!!############\tRead2:TTTTTT\tGACATTAAGCAA\t############\n", out1->str());
};

TEST_F(ExtractUMITest, BothExtract) { // Both extract test
    std::istringstream in1(readData_1);
    std::istringstream in2(readData_2);

    InputReader<PairedEndRead, PairedEndReadFastqImpl> ifp(in1, in2);
    std::shared_ptr<std::ostringstream> out1(new std::ostringstream);
    std::shared_ptr<HtsOfstream> hts_of(new HtsOfstream(out1));
    std::shared_ptr<OutputWriter> tab(new ReadBaseOutTab(hts_of));

    // reset umi struct
    std::tie(umi.seq, umi.qual) = std::make_tuple("", "");

    WriterHelper writer(tab, tab);
    while (ifp.has_next()) {
        auto i = ifp.next();
        PairedEndRead *per = dynamic_cast<PairedEndRead*>(i.get());
        eu.extract_umi(per->non_const_read_one(), umi, ':', '-', false, false);
        std::tie(umi.seq) = std::make_tuple(""); // reset umi struct
        eu.extract_umi(per->non_const_read_two(), umi, ':', '-', false, false);
        writer(*per);
    }
    ASSERT_EQ("Read1:NAAAAA\tGACATTAAGCAA\t############\tRead2:TTTTTT\tGACATTAAGCAA\t############\n", out1->str());
};


TEST_F(ExtractUMITest, QualFilt) { // Hard Filter extract test
    std::istringstream in1(readData_1);

    InputReader<SingleEndRead, SingleEndReadFastqImpl> ifp(in1);

    // reset umi struct
    std::tie(umi.seq, umi.qual, umi.qual_threshold) = std::make_tuple("", "", 20);

    auto i = ifp.next();
    SingleEndRead *ser = dynamic_cast<SingleEndRead*>(i.get());
    eu.extract_umi(ser->non_const_read_one(), umi, ':', '-');
    ASSERT_EQ(true, (ser->non_const_read_one()).getDiscard());
};


TEST_F(ExtractUMITest, AvgQualFilt) { // Average Filter extract test
    std::istringstream in1(readData_1);

    InputReader<SingleEndRead, SingleEndReadFastqImpl> ifp(in1);

    // reset umi struct
    std::tie(umi.seq, umi.qual, umi.discard, umi.qual_threshold, umi.avg_qual_threshold) = std::make_tuple("", "", false, 0, 20);

    auto i = ifp.next();
    SingleEndRead *ser = dynamic_cast<SingleEndRead*>(i.get());
    eu.extract_umi(ser->non_const_read_one(), umi, ':', '-');
    ASSERT_EQ(true, (ser->non_const_read_one()).getDiscard());
};


TEST_F(ExtractUMITest, HomopolymerFilt) { // Homopolymer filter test
    std::istringstream in2(readData_2);

    InputReader<SingleEndRead, SingleEndReadFastqImpl> ifp(in2);

    // reset umi struct
    std::tie(umi.seq, umi.qual, umi.discard, umi.avg_qual_threshold, umi.homopolymer) = std::make_tuple("", "", false, 0, true);

    auto i = ifp.next();
    SingleEndRead *ser = dynamic_cast<SingleEndRead*>(i.get());
    eu.extract_umi(ser->non_const_read_one(), umi, ':', '-');
    ASSERT_EQ(true, (ser->non_const_read_one()).getDiscard());
};


TEST_F(ExtractUMITest, NFilt) { // N filter test
    std::istringstream in1(readData_1);

    InputReader<SingleEndRead, SingleEndReadFastqImpl> ifp(in1);

    // reset umi struct
    std::tie(umi.seq, umi.qual, umi.discard, umi.homopolymer, umi.discard_n) = std::make_tuple("", "", false, false, true);

    auto i = ifp.next();
    SingleEndRead *ser = dynamic_cast<SingleEndRead*>(i.get());
    eu.extract_umi(ser->non_const_read_one(), umi, ':', '-');
    ASSERT_EQ(true, (ser->non_const_read_one()).getDiscard());
};


TEST_F(ExtractUMITest, DRAGEN) { // DRAGEN test
    std::istringstream in1(dragenData_1);
    std::istringstream in2(dragenData_2);

    InputReader<PairedEndRead, PairedEndReadFastqImpl> ifp(in1, in2);
    std::shared_ptr<std::ostringstream> out1(new std::ostringstream);
    std::shared_ptr<HtsOfstream> hts_of(new HtsOfstream(out1));
    std::shared_ptr<OutputWriter> tab(new ReadBaseOutTab(hts_of));

     // reset umi struct
    std::tie(umi.seq, umi.qual, umi.discard, umi.homopolymer, umi.discard_n) = std::make_tuple("", "", false, false, false);

    WriterHelper writer(tab, tab);
    while (ifp.has_next()) {
        auto i = ifp.next();
        PairedEndRead *per = dynamic_cast<PairedEndRead*>(i.get());
        eu.extract_umi(per->non_const_read_one(), umi, ':', '+', true, false);
        eu.extract_umi(per->non_const_read_two(), umi, ':', '+', true, false);
        eu.set_dragen(per->non_const_read_one(), umi);
        eu.set_dragen(per->non_const_read_two(), umi);
        writer(*per);
    }
    ASSERT_EQ("A00887_1:1:2:3:4:5:6:NAAAAA+TTTTTT\tGACATTAAGCAA\t############\tA00887_2:1:2:3:4:5:6:NAAAAA+TTTTTT\tGACATTAAGCAA\t############\n", out1->str());
};

TEST_F(ExtractUMITest, TagExtract) { // Tag extract test
    std::istringstream in1(readData_1);
    std::istringstream in2(readData_2);

    InputReader<PairedEndRead, PairedEndReadFastqImpl> ifp(in1, in2);
    std::shared_ptr<std::ostringstream> out1(new std::ostringstream);
    std::shared_ptr<HtsOfstream> hts_of(new HtsOfstream(out1));
    std::shared_ptr<OutputWriter> tab(new ReadBaseOutTab(hts_of));

    // reset umi struct
    std::tie(umi.seq, umi.qual, umi.add_as_tag) = std::make_tuple("", "", true);

    WriterHelper writer(tab, tab);
    while (ifp.has_next()) {
        auto i = ifp.next();
        PairedEndRead *per = dynamic_cast<PairedEndRead*>(i.get());
        eu.extract_umi( per->non_const_read_one(), umi, ':', '-', false, true );
        eu.set_tag( per->non_const_read_one(), umi );
        std::tie(umi.seq) = std::make_tuple(""); // reset umi struct
        eu.extract_umi( per->non_const_read_two(), umi, ':', '-', false, true );
        eu.set_tag( per->non_const_read_two(), umi );
        writer(*per);
    }
    ASSERT_EQ("Read1\tGACATTAAGCAA\t############\tRead2\tGACATTAAGCAA\t############\tRX:Z:NAAAAA\tRX:Z:TTTTTT\n", out1->str());
};