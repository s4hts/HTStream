#include <gtest/gtest.h>
#include "ioHandler.h"
#include <sstream>
#include <iostream>

std::string readData("@M03610:15:000000000-ALE1U:1:1101:17702:1965 1:N:0:1\nCCTTTTGTTTGTATTGTTTCTCTTAAGCATATTTCTTAATTACTTAACTAACTACCTCAAATGAAATTTGTAAAAAATCCATCAATTTTAACCGCTTGTTTGCTTGCCAGTGTCCTCGGAGGCATTGGGACACTCGCCTACTCCCCTTTTG\n+\nAAAAAF5CFFFF6FGFGGGGGGHHGHHFBHFFHHHHHBEGHHHHHEGGHHHHHFGFGGFEFGHBFBHHHHHEEGFFEGHHFHHEBGHHHEGHHGEGGEHHHGHHHEHHFFG3GHHHHGCDEGGGHHDGEHHGHHHFEEGGFHHHHGHHHH1\n@M03610:15:000000000-ALE1U:1:1101:13392:1966 1:N:0:1\nGTCGCCGGATCTAATGCACTTGTCGCTTCATCGCATAATAAAACGTGCGGATCACTTGCCAACGCACGGGCAATCGCTACACGTTGTTTCTGACCACCGGATAAATTGCTCGGATAAACGTCTTTACGTTCAGTTAAGCCGACCAATTCAA\n+\n>AA3>A>22ADFC5FBGGGGGCGHGHGHGGHHHGGGGFHFFHFGGHHHGGGGEHHHHFHHGHGGGGGGGCGGFHHGGHGGHHGHHHGHHHHHHHHGHGGGGGGFHHHHHHHGGGGGHHHGHHHGHHHGHHHHHHHGGEGFC?DGGHHHHG0\n@M03610:15:000000000-ALE1U:1:1101:19181:1996 1:N:0:1\nATTTGTGGCAATGGCGAGTTTGGGGCGATACCAATATCCACACAATAGCCAGTAAACTTAACGTTAGCGGGCAAGGAAATTGGCACCAAAATTTGGTACAACTCAACAAATTAAGTGGTGAATTAGGCAAAATTCACCAGAAAGGCGTGTA\n+\n>AAA3DCAAFFFB4EAEFGGGE2AFGGGGFGHBGHFHHHFHHGGGHHHHHHHHHHHHHHHHHGHHEGHCGGGGGGGHHHHHFHHGHGHHHHHHHHGHHHHGHHHGHHEHGHHHHHGHGHGHHHHHHHHHHGHHHHHHHHH0GHGHHGGGGF\n@M03610:15:000000000-ALE1U:1:1101:18865:1996 1:N:0:1\nGCTACAAAGATAATGATAGCGCTTATAGTTAGCTCGGTTTGCTAATTCGGTTTGACTAATATTATCCTTTAAATAATCTGTTTTCGCTTTTTTAGTCGGCGAAGGGGAACAAGCAAATAAAGTTGTCAGTGTTGAAATTAAAAAAATAAGT\n+\n>3A3>4@4C4FFBGFGGGGGGGGGGHGHHHHGGHHGGFGGGHHGHHHHGHGGGGEHHGGHHHHHHFHHHHFFHHHFHHHHHHHHHGGGHGHHGHGHHDCEGCEG?EGGGGHHHHHHHHHHHGHHHHGHHHHHHHHFFGHHEHHHGGHHHHH\n@M03610:15:000000000-ALE1U:1:1101:18806:2006 1:N:0:1\nGCGGATGGCTCGCCAAAAGTTGGCTCTTTACGAAAATTATTTTGTGACATTTAACTTTCCAATTTTAGGTGAAAATATAATGGGAAATTTCCCGCTGAATAAATTATCCCCATTTTAAAGTATTTCTCTTTAATCAGCACACGATTTTTCA\n+\n>>111>1AFFFDCEECGGGGGGHHHHHHHHHGFGC?GHFHHHHGGHHHHHHHHHHHHHHHHGHHHHHHGHHFHHHHHHHHGHGHHHHHHHBHHGGGGGHHHHHHHHHGHHGHHHHHHHHHEGHGHHHHGHHHHHHHFGGHGHGGHHGHHF2\n");

TEST(testIoHandle, parseSingleReadFastq) {
    std::istringstream in(readData);
    inputFastqSingle ifs(in);
    size_t read_count = 0;
    for(auto i = ifs.begin(); i != ifs.end(); ++i) {
        std::cout << i->get_read().get_qual() << std::endl;
        read_count++;
    }
    ASSERT_EQ(read_count, 5);
}    
