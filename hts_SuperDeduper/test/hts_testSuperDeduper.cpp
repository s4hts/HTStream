#include "gtest/gtest.h"
#include <sstream>
#include <iostream>
#include "hts_SuperDeduper.h"

class SDTest : public ::testing::Test {
    public:
        const std::string readData = "@M03610:15:000000000-ALE1U:1:1101:17702:1965 1:N:0:1\nCCTTTTGTTTGTATTGTTTCTCTTAAGCATATTTCTTAATTACTTAACTAACTACCTCAAATGAAATTTGTAAAAAATCCATCAATTTTAACCGCTTGTTTGCTTGCCAGTGTCCTCGGAGGCATTGGGACACTCGCCTACTCCCCTTTTG\n+\nAAAAAF5CFFFF6FGFGGGGGGHHGHHFBHFFHHHHHBEGHHHHHEGGHHHHHFGFGGFEFGHBFBHHHHHEEGFFEGHHFHHEBGHHHEGHHGEGGEHHHGHHHEHHFFG3GHHHHGCDEGGGHHDGEHHGHHHFEEGGFHHHHGHHHH1\n@M03610:15:000000000-ALE1U:1:1101:17702:1965 1:N:0:1\nCCTTTTGTTTGTATTGTTTCTCTTAAGCATATTTCTTAATTACTTAACTAACTACCTCAAATGAAATTTGTAAAAAATCCATCAATTTTAACCGCTTGTTTGCTTGCCAGTGTCCTCGGAGGCATTGGGACACTCGCCTACTCCCCTTTTG\n+\nBAAAAF5CFFFF6FGFGGGGGGHHGHHFBHFFHHHHHBEGHHHHHEGGHHHHHFGFGGFEFGHBFBHHHHHEEGFFEGHHFHHEBGHHHEGHHGEGGEHHHGHHHEHHFFG3GHHHHGCDEGGGHHDGEHHGHHHFEEGGFHHHHGHHHH1\n@M03610:15:000000000-ALE1U:1:1101:17702:1965 1:N:0:1\nCCTTTTGTTTGTATTGTTTCTCTTAAGCATATTTCTTAATTACTTAACTAACTACCTCAAATGAAATTTGTAAAAAATCCATCAATTTTAACCGCTTGTTTGCTTGCCAGTGTCCTCGGAGGCATTGGGACACTCGCCTACTCCCCTTTTG\n+\nCAAAAF5CFFFF6FGFGGGGGGHHGHHFBHFFHHHHHBEGHHHHHEGGHHHHHFGFGGFEFGHBFBHHHHHEEGFFEGHHFHHEBGHHHEGHHGEGGEHHHGHHHEHHFFG3GHHHHGCDEGGGHHDGEHHGHHHFEEGGFHHHHGHHHH1\n@M03610:15:000000000-ALE1U:1:1101:17702:1965 1:N:0:1\nCNNNNNNNNNNNNNNNNNNNNCTTAAGCATATTTCTTAATTACTTAACTAACTACCTCAAATGAAATTTGTAAAAAATCCATCAATTTTAACCGCTTGTTTGCTTGCCAGTGTCCTCGGAGGCATTGGGACACTCGCCTACTCCCCTTTTG\n+\nAAAAAF5CFFFF6FGFGGGGGGHHGHHFBHFFHHHHHBEGHHHHHEGGHHHHHFGFGGFEFGHBFBHHHHHEEGFFEGHHFHHEBGHHHEGHHGEGGEHHHGHHHEHHFFG3GHHHHGCDEGGGHHDGEHHGHHHFEEGGFHHHHGHHHH1\n";

};

TEST_F(SDTest, HashMapLoadTest) {

    std::istringstream in1(readData);
    std::istringstream in2(readData);
    size_t start = 5;
    size_t length = 5;
    double avg_auto_write = 100;
    BitMap read_map;
    SuperDeduperCounters counter("/dev/null", false, false, "hts_SuperDeduper", "");

    InputReader<PairedEndRead, PairedEndReadFastqImpl> ifp(in1, in2);
    std::shared_ptr<std::ostringstream> out1(new std::ostringstream);
    std::shared_ptr<HtsOfstream> hts_of(new HtsOfstream(out1));
    std::shared_ptr<OutputWriter> tab(new ReadBaseOutTab(hts_of));


    load_map(ifp, counter, read_map, tab, tab, avg_auto_write, 3, start, length, 100);
    std::cout << read_map.size() << '\n';
    ASSERT_EQ(read_map.size(), 1);
    ASSERT_EQ(counter.TotalFragmentsInput, 4);
    ASSERT_EQ(counter.Duplicate, 2);
    ASSERT_EQ(counter.Ignored, 1);


};
