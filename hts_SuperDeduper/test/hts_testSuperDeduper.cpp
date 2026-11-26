#include "gtest/gtest.h"
#include <sstream>
#include <iostream>
#include <boost/program_options.hpp>
#include "hts_SuperDeduper.h"

namespace po = boost::program_options;

class SDTest : public ::testing::Test {
public:
    po::variables_map vm;
    const std::string readData = "@M03610:15:000000000-ALE1U:1:1101:17702:1965 1:N:0:1\nCCTTTTGTTTGTATTGTTTCTCTTAAGCATATTTCTTAATTACTTAACTAACTACCTCAAATGAAATTTGTAAAAAATCCATCAATTTTAACCGCTTGTTTGCTTGCCAGTGTCCTCGGAGGCATTGGGACACTCGCCTACTCCCCTTTTG\n+\nAAAAAF5CFFFF6FGFGGGGGGHHGHHFBHFFHHHHHBEGHHHHHEGGHHHHHFGFGGFEFGHBFBHHHHHEEGFFEGHHFHHEBGHHHEGHHGEGGEHHHGHHHEHHFFG3GHHHHGCDEGGGHHDGEHHGHHHFEEGGFHHHHGHHHH1\n@M03610:15:000000000-ALE1U:1:1101:17702:1965 1:N:0:1\nCCTTTTGTTTGTATTGTTTCTCTTAAGCATATTTCTTAATTACTTAACTAACTACCTCAAATGAAATTTGTAAAAAATCCATCAATTTTAACCGCTTGTTTGCTTGCCAGTGTCCTCGGAGGCATTGGGACACTCGCCTACTCCCCTTTTG\n+\nBAAAAF5CFFFF6FGFGGGGGGHHGHHFBHFFHHHHHBEGHHHHHEGGHHHHHFGFGGFEFGHBFBHHHHHEEGFFEGHHFHHEBGHHHEGHHGEGGEHHHGHHHEHHFFG3GHHHHGCDEGGGHHDGEHHGHHHFEEGGFHHHHGHHHH1\n@M03610:15:000000000-ALE1U:1:1101:17702:1965 1:N:0:1\nCCTTTTGTTTGTATTGTTTCTCTTAAGCATATTTCTTAATTACTTAACTAACTACCTCAAATGAAATTTGTAAAAAATCCATCAATTTTAACCGCTTGTTTGCTTGCCAGTGTCCTCGGAGGCATTGGGACACTCGCCTACTCCCCTTTTG\n+\nCAAAAF5CFFFF6FGFGGGGGGHHGHHFBHFFHHHHHBEGHHHHHEGGHHHHHFGFGGFEFGHBFBHHHHHEEGFFEGHHFHHEBGHHHEGHHGEGGEHHHGHHHEHHFFG3GHHHHGCDEGGGHHDGEHHGHHHFEEGGFHHHHGHHHH1\n@M03610:15:000000000-ALE1U:1:1101:17702:1965 1:N:0:1\nCNNNNNNNNNNNNNNNNNNNNCTTAAGCATATTTCTTAATTACTTAACTAACTACCTCAAATGAAATTTGTAAAAAATCCATCAATTTTAACCGCTTGTTTGCTTGCCAGTGTCCTCGGAGGCATTGGGACACTCGCCTACTCCCCTTTTG\n+\nAAAAAF5CFFFF6FGFGGGGGGHHGHHFBHFFHHHHHBEGHHHHHEGGHHHHHFGFGGFEFGHBFBHHHHHEEGFFEGHHFHHEBGHHHEGHHGEGGEHHHGHHHEHHFFG3GHHHHGCDEGGGHHDGEHHGHHHFEEGGFHHHHGHHHH1\n";
    const std::string umiData1 = "@M03610:15:000000000-ALE1U:1:1101:17702:1965:CCTTTT 1:N:0:1\nGTTTGTATTGTTTCTCTTAAGCATATTTCTTAATTACTTAACTAACTACCTCAAATGAAATTTGTAAAAAATCCATCAATTTTAACCGCTTGTTTGCTTGCCAGTGTCCTCGGAGGCATTGGGACACTCGCCTACTCCCCTTTTG\n+\n5CFFFF6FGFGGGGGGHHGHHFBHFFHHHHHBEGHHHHHEGGHHHHHFGFGGFEFGHBFBHHHHHEEGFFEGHHFHHEBGHHHEGHHGEGGEHHHGHHHEHHFFG3GHHHHGCDEGGGHHDGEHHGHHHFEEGGFHHHHGHHHH1\n@M03610:15:000000000-ALE1U:1:1101:17702:1965:CCTTTT 1:N:0:1\nGTTTGTATTGTTTCTCTTAAGCATATTTCTTAATTACTTAACTAACTACCTCAAATGAAATTTGTAAAAAATCCATCAATTTTAACCGCTTGTTTGCTTGCCAGTGTCCTCGGAGGCATTGGGACACTCGCCTACTCCCCTTTTG\n+\n5CFFFF6FGFGGGGGGHHGHHFBHFFHHHHHBEGHHHHHEGGHHHHHFGFGGFEFGHBFBHHHHHEEGFFEGHHFHHEBGHHHEGHHGEGGEHHHGHHHEHHFFG3GHHHHGCDEGGGHHDGEHHGHHHFEEGGFHHHHGHHHH1\n@M03610:15:000000000-ALE1U:1:1101:17702:1965:CCTTTT 1:N:0:1\nGTTTGTATTGTTTCTCTTAAGCATATTTCTTAATTACTTAACTAACTACCTCAAATGAAATTTGTAAAAAATCCATCAATTTTAACCGCTTGTTTGCTTGCCAGTGTCCTCGGAGGCATTGGGACACTCGCCTACTCCCCTTTTG\n+\n5CFFFF6FGFGGGGGGHHGHHFBHFFHHHHHBEGHHHHHEGGHHHHHFGFGGFEFGHBFBHHHHHEEGFFEGHHFHHEBGHHHEGHHGEGGEHHHGHHHEHHFFG3GHHHHGCDEGGGHHDGEHHGHHHFEEGGFHHHHGHHHH1\n@M03610:15:000000000-ALE1U:1:1101:17702:1965:CNNNNN 1:N:0:1\nNNNNNNNNNNNNNNNCTTAAGCATATTTCTTAATTACTTAACTAACTACCTCAAATGAAATTTGTAAAAAATCCATCAATTTTAACCGCTTGTTTGCTTGCCAGTGTCCTCGGAGGCATTGGGACACTCGCCTACTCCCCTTTTG\n+\n5CFFFF6FGFGGGGGGHHGHHFBHFFHHHHHBEGHHHHHEGGHHHHHFGFGGFEFGHBFBHHHHHEEGFFEGHHFHHEBGHHHEGHHGEGGEHHHGHHHEHHFFG3GHHHHGCDEGGGHHDGEHHGHHHFEEGGFHHHHGHHHH1\n";
    const std::string umiData2 = "@M03610:15:000000000-ALE1U:1:1101:17702:1965:CCTTTT+TTTCCC 1:N:0:1\nGTTTGTATTGTTTCTCTTAAGCATATTTCTTAATTACTTAACTAACTACCTCAAATGAAATTTGTAAAAAATCCATCAATTTTAACCGCTTGTTTGCTTGCCAGTGTCCTCGGAGGCATTGGGACACTCGCCTACTCCCCTTTTG\n+\n5CFFFF6FGFGGGGGGHHGHHFBHFFHHHHHBEGHHHHHEGGHHHHHFGFGGFEFGHBFBHHHHHEEGFFEGHHFHHEBGHHHEGHHGEGGEHHHGHHHEHHFFG3GHHHHGCDEGGGHHDGEHHGHHHFEEGGFHHHHGHHHH1\n@M03610:15:000000000-ALE1U:1:1101:17702:1965:CCTTTT+TTTCCC 1:N:0:1\nGTTTGTATTGTTTCTCTTAAGCATATTTCTTAATTACTTAACTAACTACCTCAAATGAAATTTGTAAAAAATCCATCAATTTTAACCGCTTGTTTGCTTGCCAGTGTCCTCGGAGGCATTGGGACACTCGCCTACTCCCCTTTTG\n+\n5CFFFF6FGFGGGGGGHHGHHFBHFFHHHHHBEGHHHHHEGGHHHHHFGFGGFEFGHBFBHHHHHEEGFFEGHHFHHEBGHHHEGHHGEGGEHHHGHHHEHHFFG3GHHHHGCDEGGGHHDGEHHGHHHFEEGGFHHHHGHHHH1\n@M03610:15:000000000-ALE1U:1:1101:17702:1965:CCTTTT+TTTCCC 1:N:0:1\nGTTTGTATTGTTTCTCTTAAGCATATTTCTTAATTACTTAACTAACTACCTCAAATGAAATTTGTAAAAAATCCATCAATTTTAACCGCTTGTTTGCTTGCCAGTGTCCTCGGAGGCATTGGGACACTCGCCTACTCCCCTTTTG\n+\n5CFFFF6FGFGGGGGGHHGHHFBHFFHHHHHBEGHHHHHEGGHHHHHFGFGGFEFGHBFBHHHHHEEGFFEGHHFHHEBGHHHEGHHGEGGEHHHGHHHEHHFFG3GHHHHGCDEGGGHHDGEHHGHHHFEEGGFHHHHGHHHH1\n@M03610:15:000000000-ALE1U:1:1101:17702:1965:CNNNNN+TTTCCC 1:N:0:1\nNNNNNNNNNNNNNNNCTTAAGCATATTTCTTAATTACTTAACTAACTACCTCAAATGAAATTTGTAAAAAATCCATCAATTTTAACCGCTTGTTTGCTTGCCAGTGTCCTCGGAGGCATTGGGACACTCGCCTACTCCCCTTTTG\n+\n5CFFFF6FGFGGGGGGHHGHHFBHFFHHHHHBEGHHHHHEGGHHHHHFGFGGFEFGHBFBHHHHHEEGFFEGHHFHHEBGHHHEGHHGEGGEHHHGHHHEHHFFG3GHHHHGCDEGGGHHDGEHHGHHHFEEGGFHHHHGHHHH1\n";
    const std::string recordData = "@M03610:15:000000000-ALE1U:1:1101:17702:1965 1:N:0:1\nACGTACGTAGCTAGCTAGGATCGATCGATC\n+\n++++++++++++++++++++++++++++++\n@M03610:15:000000000-ALE1U:1:1101:17702:1965 1:N:0:1\nACGTACGTAGCTAGCTAGGATCGATCGATC\n+\n555555555555555555555555555555\n@M03610:15:000000000-ALE1U:1:1101:17702:1965 1:N:0:1\nACGTACGTAGCTAGCTAGGATCGATCGATC\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n@M03610:15:000000000-ALE1U:1:1101:17702:1965 1:N:0:1\nACGTACGTAGCTAGCTAGGATCGATCGATC\n+\n000000000000000000000000000000\n@M03610:15:000000000-ALE1U:1:1101:17702:1965 1:N:0:1\nACGTANNNNNCTAGCTAGGATCGATCGATC\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n";

    SuperDeduper sd;
};

TEST_F(SDTest, HashMapLoadTest) {
    std::istringstream in1(readData);
    std::istringstream in2(readData);
    size_t start = 5;
    size_t length = 5;
    double avg_auto_write = 100;
    BitMap read_map;
    SuperDeduperCounters counter("hts_SuperDeduper", vm);

    InputReader<PairedEndRead, PairedEndReadFastqImpl> ifp(in1, in2);
    std::shared_ptr<std::ostringstream> out1(new std::ostringstream);
    std::shared_ptr<HtsOfstream> hts_of(new HtsOfstream(out1));
    std::shared_ptr<OutputWriter> tab(new ReadBaseOutTab(hts_of));

    sd.load_map(ifp, counter, tab, tab, avg_auto_write, 3, start, length, 100, false);
    std::cout << read_map.size() << '\n';
    ASSERT_EQ(sd.read_map.size(), 1u);
    ASSERT_EQ(counter.TotalFragmentsInput, 4u);
    ASSERT_EQ(counter.Duplicate, 2u);
    ASSERT_EQ(counter.Ignored, 1u);
};


TEST_F(SDTest, RecordDataDedupTest) {
    // Baseline: same data without record, only the best read (rec_c) is kept.
    std::istringstream in1(recordData);
    std::istringstream in2(recordData);
    size_t start = 5;
    size_t length = 5;
    double avg_auto_write = 30;
    SuperDeduperCounters counter("hts_SuperDeduper", vm);

    InputReader<PairedEndRead, PairedEndReadFastqImpl> ifp(in1, in2);
    std::shared_ptr<std::ostringstream> out1(new std::ostringstream);
    std::shared_ptr<HtsOfstream> hts_of(new HtsOfstream(out1));
    std::shared_ptr<OutputWriter> tab(new ReadBaseOutTab(hts_of));

    sd.load_map(ifp, counter, tab, tab, avg_auto_write, 3, start, length, 100, false);

    ASSERT_EQ(sd.read_map.size(), 1u);
    ASSERT_EQ(sd.read_map.begin()->second, nullptr);
    ASSERT_EQ(counter.TotalFragmentsInput, 5u);
    ASSERT_EQ(counter.TotalFragmentsOutput, 1u);
    ASSERT_EQ(counter.Duplicate, 3u);
    ASSERT_EQ(counter.Ignored, 1u);
};

TEST_F(SDTest, RecordTest) {
    // Record mode: every fragment is written exactly once, dup/ignored counts match dedup mode.
    std::istringstream in1(recordData);
    std::istringstream in2(recordData);
    size_t start = 5;
    size_t length = 5;
    double avg_auto_write = 30;
    SuperDeduperCounters counter("hts_SuperDeduper", vm);

    InputReader<PairedEndRead, PairedEndReadFastqImpl> ifp(in1, in2);
    std::shared_ptr<std::ostringstream> out1(new std::ostringstream);
    std::shared_ptr<HtsOfstream> hts_of(new HtsOfstream(out1));
    std::shared_ptr<OutputWriter> tab(new ReadBaseOutTab(hts_of));

    sd.load_map(ifp, counter, tab, tab, avg_auto_write, 3, start, length, 100, true);

    // rec_b, rec_c, rec_d, rec_e written in-loop; rec_a still held in the map
    ASSERT_EQ(counter.TotalFragmentsOutput, 4u);
    ASSERT_EQ(sd.read_map.size(), 1u);
    ASSERT_NE(sd.read_map.begin()->second, nullptr);

    // Flush readmap, as do_app does
    for (auto const &i : sd.read_map) {
        if (i.second.get() != nullptr) {
            counter.output(*i.second.get());
        }
    }

    ASSERT_EQ(counter.TotalFragmentsInput, 5u);
    ASSERT_EQ(counter.TotalFragmentsOutput, 5u);
    ASSERT_EQ(counter.Duplicate, 3u);
    ASSERT_EQ(counter.Ignored, 1u);
};

TEST_F(SDTest, UMITest_Both) {
    std::istringstream in1(umiData1);
    size_t start = 5;
    size_t length = 5;
    double avg_auto_write = 100;
    BitMap read_map;
    SuperDeduperCounters counter("hts_SuperDeduper", vm);

    InputReader<SingleEndRead, SingleEndReadFastqImpl> ifp(in1);
    std::shared_ptr<std::ostringstream> out1(new std::ostringstream);
    std::shared_ptr<HtsOfstream> hts_of(new HtsOfstream(out1));
    std::shared_ptr<OutputWriter> tab(new ReadBaseOutTab(hts_of));

    sd.load_map(ifp, counter, tab, tab, avg_auto_write, 3, start, length, 100, false, 30, ':');
    std::cout << read_map.size() << '\n';
    ASSERT_EQ(sd.read_map.size(), 1u);
    ASSERT_EQ(counter.TotalFragmentsInput, 4u);
    ASSERT_EQ(counter.Duplicate, 2u);
    ASSERT_EQ(counter.Ignored, 1u);


};

TEST_F(SDTest, UMITest_Paired) {
    std::istringstream in1(umiData2);
    size_t start = 5;
    size_t length = 5;
    double avg_auto_write = 100;
    BitMap read_map;
    SuperDeduperCounters counter("hts_SuperDeduper", vm);

    InputReader<SingleEndRead, SingleEndReadFastqImpl> ifp(in1);
    std::shared_ptr<std::ostringstream> out1(new std::ostringstream);
    std::shared_ptr<HtsOfstream> hts_of(new HtsOfstream(out1));
    std::shared_ptr<OutputWriter> tab(new ReadBaseOutTab(hts_of));

    sd.load_map(ifp, counter, tab, tab, avg_auto_write, 3, start, length, 100, false, 30, ':');
    std::cout << read_map.size() << '\n';
    ASSERT_EQ(sd.read_map.size(), 1u);
    ASSERT_EQ(counter.TotalFragmentsInput, 4u);
    ASSERT_EQ(counter.Duplicate, 2u);
    ASSERT_EQ(counter.Ignored, 1u);

};

TEST_F(SDTest, UMITest_Tag) {
    std::istringstream in1(readData);
    size_t start = 5;
    size_t length = 5;
    double avg_auto_write = 100;
    BitMap read_map;
    SuperDeduperCounters counter("hts_SuperDeduper", vm);

    InputReader<SingleEndRead, SingleEndReadFastqImpl> ifp(in1);
    std::shared_ptr<std::ostringstream> out1(new std::ostringstream);
    std::shared_ptr<HtsOfstream> hts_of(new HtsOfstream(out1));
    std::shared_ptr<OutputWriter> tab(new ReadBaseOutTab(hts_of));

    WriterHelper writer(tab, tab);
    while (ifp.has_next()) {
        auto i = ifp.next();
        SingleEndRead *ser = dynamic_cast<SingleEndRead*>(i.get());
        ser->non_const_read_one().add_comment("RX:Z:CCTTTT+TTTCCC");
        writer(*ser);
    }

    std::string umiData3 = out1->str();
    std::istringstream in2(umiData3);
    InputReader<ReadBase, TabReadImpl> ift(in2);

    sd.load_map(ift, counter, tab, tab, avg_auto_write, 3, start, length, 100, false, 30, ' ', 0, true);

    std::cout << read_map.size() << '\n';
    ASSERT_EQ(sd.read_map.size(), 1u);
    ASSERT_EQ(counter.TotalFragmentsInput, 4u);
    ASSERT_EQ(counter.Duplicate, 2u);
    ASSERT_EQ(counter.Ignored, 1u);

};