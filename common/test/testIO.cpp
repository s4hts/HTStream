#include <gtest/gtest.h>
#include "ioHandler.h"
#include <sstream>
#include <iostream>

extern template class InputReader<SingleEndRead, SingleEndReadFastqImpl>;
extern template class InputReader<SingleEndRead, FastaReadImpl>;
extern template class InputReader<PairedEndRead, PairedEndReadFastqImpl>;
extern template class InputReader<PairedEndRead, InterReadImpl>;
extern template class InputReader<ReadBase, TabReadImpl>;

class ReadsTest : public ::testing::Test {
public:
    const std::string readData = "@M03610:15:000000000-ALE1U:1:1101:17702:1965 1:N:0:1\nCCTTTTGTTTGTATTGTTTCTCTTAAGCATATTTCTTAATTACTTAACTAACTACCTCAAATGAAATTTGTAAAAAATCCATCAATTTTAACCGCTTGTTTGCTTGCCAGTGTCCTCGGAGGCATTGGGACACTCGCCTACTCCCCTTTTG\n+\nAAAAAF5CFFFF6FGFGGGGGGHHGHHFBHFFHHHHHBEGHHHHHEGGHHHHHFGFGGFEFGHBFBHHHHHEEGFFEGHHFHHEBGHHHEGHHGEGGEHHHGHHHEHHFFG3GHHHHGCDEGGGHHDGEHHGHHHFEEGGFHHHHGHHHH1\n@M03610:15:000000000-ALE1U:1:1101:13392:1966 1:N:0:1\nGTCGCCGGATCTAATGCACTTGTCGCTTCATCGCATAATAAAACGTGCGGATCACTTGCCAACGCACGGGCAATCGCTACACGTTGTTTCTGACCACCGGATAAATTGCTCGGATAAACGTCTTTACGTTCAGTTAAGCCGACCAATTCAA\n+\n>AA3>A>22ADFC5FBGGGGGCGHGHGHGGHHHGGGGFHFFHFGGHHHGGGGEHHHHFHHGHGGGGGGGCGGFHHGGHGGHHGHHHGHHHHHHHHGHGGGGGGFHHHHHHHGGGGGHHHGHHHGHHHGHHHHHHHGGEGFC?DGGHHHHG0\n@M03610:15:000000000-ALE1U:1:1101:19181:1996 1:N:0:1\nATTTGTGGCAATGGCGAGTTTGGGGCGATACCAATATCCACACAATAGCCAGTAAACTTAACGTTAGCGGGCAAGGAAATTGGCACCAAAATTTGGTACAACTCAACAAATTAAGTGGTGAATTAGGCAAAATTCACCAGAAAGGCGTGTA\n+\n>AAA3DCAAFFFB4EAEFGGGE2AFGGGGFGHBGHFHHHFHHGGGHHHHHHHHHHHHHHHHHGHHEGHCGGGGGGGHHHHHFHHGHGHHHHHHHHGHHHHGHHHGHHEHGHHHHHGHGHGHHHHHHHHHHGHHHHHHHHH0GHGHHGGGGF\n@M03610:15:000000000-ALE1U:1:1101:18865:1996 1:N:0:1\nGCTACAAAGATAATGATAGCGCTTATAGTTAGCTCGGTTTGCTAATTCGGTTTGACTAATATTATCCTTTAAATAATCTGTTTTCGCTTTTTTAGTCGGCGAAGGGGAACAAGCAAATAAAGTTGTCAGTGTTGAAATTAAAAAAATAAGT\n+\n>3A3>4@4C4FFBGFGGGGGGGGGGHGHHHHGGHHGGFGGGHHGHHHHGHGGGGEHHGGHHHHHHFHHHHFFHHHFHHHHHHHHHGGGHGHHGHGHHDCEGCEG?EGGGGHHHHHHHHHHHGHHHHGHHHHHHHHFFGHHEHHHGGHHHHH\n@M03610:15:000000000-ALE1U:1:1101:18806:2006 1:N:0:1\nGCGGATGGCTCGCCAAAAGTTGGCTCTTTACGAAAATTATTTTGTGACATTTAACTTTCCAATTTTAGGTGAAAATATAATGGGAAATTTCCCGCTGAATAAATTATCCCCATTTTAAAGTATTTCTCTTTAATCAGCACACGATTTTTCA\n+\n>>111>1AFFFDCEECGGGGGGHHHHHHHHHGFGC?GHFHHHHGGHHHHHHHHHHHHHHHHGHHHHHHGHHFHHHHHHHHGHGHHHHHHHBHHGGGGGHHHHHHHHHGHHGHHHHHHHHHEGHGHHHHGHHHHHHHFGGHGHGGHHGHHF2\n";
    const std::string readTabData = "TESTOne\tACTG\t####\tTESTOne\tACTG\t####\nTestTwo\tACTG\t####\n";
    const std::string readDataInter = "@Whatever1\nAAAAAAA\n+\n#######\n@Whatever2\nGGGCCCT\n+\nBBBBBBB\n";
    const std::string fastaFile = ">1\nTACCCACACA\nTTTTACACACAACAC\n>2\nACGATGACA\n";
    const std::string fastaString = "TACCCACACA,TTTTACACACAACAC,ACGATGACA";
};

TEST_F(ReadsTest, fasta_to_read) {
    std::istringstream in1(fastaFile);
    InputReader<SingleEndRead, FastaReadImpl> ifs(in1);
    size_t read_count = 0;
    while(ifs.has_next()) {
        auto r = ifs.next();
        std::string s(r->get_read().get_seq() );
        if (read_count == 0) {
            ASSERT_EQ("TACCCACACATTTTACACACAACAC", s);
        } else {
            ASSERT_EQ("ACGATGACA", s);
        }
        ++read_count;
    }
    ASSERT_EQ(read_count, 2u);
};


TEST_F(ReadsTest, parseTabRead) {
    std::istringstream in1(readTabData);
    InputReader<ReadBase, TabReadImpl> ifs(in1);
    size_t read_count = 0;
    while(ifs.has_next()) {
        auto r = ifs.next();
        ++read_count;
    }
    ASSERT_EQ(read_count, 2u);
}

TEST_F(ReadsTest, testInterRead) {
    std::istringstream in1(readDataInter);
    InputReader<PairedEndRead, InterReadImpl> ifs(in1);
    size_t read_count = 0;
    while(ifs.has_next()) {
        auto r = ifs.next();
        ++read_count;
    }

    ASSERT_EQ(read_count, 1u);
}


TEST_F(ReadsTest, parseSingleReadFastq) {
    std::istringstream in1(readData);
    InputReader<SingleEndRead, SingleEndReadFastqImpl> ifs(in1);
    size_t read_count = 0;
    while(ifs.has_next()) {
        auto r = ifs.next();
        read_count++;
    }
    ASSERT_EQ(read_count, 5u);
}

TEST_F(ReadsTest, parsePairedReadFastq) {
    std::istringstream in1(readData);
    std::istringstream in2(readData);
    size_t read_count = 0;

    InputReader<PairedEndRead, PairedEndReadFastqImpl> ifp(in1, in2);
    while(ifp.has_next()) {
        auto i = ifp.next();
        std::cout << i->get_read_one().get_qual() << std::endl;
        read_count++;
    }
    ASSERT_EQ(read_count, 5u);
}

TEST_F(ReadsTest, testWriteFastqSingle) {

    std::istringstream in1(readData);
    InputReader<SingleEndRead, SingleEndReadFastqImpl> ifs(in1);
    std::shared_ptr<std::ostringstream> out1(new std::ostringstream());

    {
        std::shared_ptr<HtsOfstream> hts_of(new HtsOfstream(out1));
        std::unique_ptr<OutputWriter> se(new SingleEndReadOutFastq(hts_of));

        while(ifs.has_next()) {
            auto r = ifs.next();
            se->write(*r);
        }
    }
    ASSERT_EQ(readData, out1->str());
}

TEST_F(ReadsTest, testTabWrite) {
    std::istringstream in1(readTabData);
    InputReader<ReadBase, TabReadImpl> ifs(in1);
    std::shared_ptr<std::ostringstream> out1(new std::ostringstream());

    {
        std::shared_ptr<HtsOfstream> hts_of(new HtsOfstream(out1));
        std::shared_ptr<OutputWriter> ofs(new ReadBaseOutTab(hts_of));
        WriterHelper writer = WriterHelper(ofs, ofs);
        while(ifs.has_next()) {
            auto r = ifs.next();
            writer(*r);
        }
    }

    ASSERT_EQ(readTabData, out1->str());
}

TEST_F(ReadsTest, testInterWrite) {
    std::istringstream in1(readDataInter);
    InputReader<PairedEndRead, InterReadImpl> ifs(in1);
    std::shared_ptr<std::ostringstream> out1(new std::ostringstream());

    {
        std::shared_ptr<HtsOfstream> hts_of(new HtsOfstream(out1));
        std::unique_ptr<OutputWriter> ofs(new PairedEndReadOutInter(hts_of));
        while(ifs.has_next()) {
            auto r = ifs.next();
            ofs->write(*r);
        }
    }

    ASSERT_EQ(readDataInter, out1->str());
}

TEST_F(ReadsTest, testString2Fasta) {

    std::istringstream fa_to_read(string2fasta(fastaString, "seq"));
    InputReader<SingleEndRead, FastaReadImpl> faReader(fa_to_read);

    auto a = faReader.next();
    auto b = faReader.next();
    Read readSeq = b->get_read();

    ASSERT_EQ("TTTTACACACAACAC", readSeq.get_seq());
    ASSERT_EQ("seq2", readSeq.get_id());
}
