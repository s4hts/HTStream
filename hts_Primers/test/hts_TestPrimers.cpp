#include "gtest/gtest.h"
#include <sstream>
#include <iostream>
#include <boost/program_options.hpp>
#include "hts_Primers.h"

namespace po = boost::program_options;

class Primer : public ::testing::Test {
public:
    po::variables_map vm;
    const std::string primer_1 = "TTCATTAAAAATTGAATTGACATTAACCT";
    const std::string read_seq_1 = "CGGGTTTCATTAAAAATTGAATTGACATTAACCTATAAAAATAGGCGTCGAGGCCCTTTCGTCTTCTATCGGAGCTCCAAGACCGCCTCGGCGTGAAGGTGGTGATAGCGCCCGGAAGAGAGTCAATTCAGGGTGGTGAATACTCTAGATC";

    const std::string primer_5p = "GTGYCAGCMGCCGCGGTAA";
    const std::string primer_3p = "GGACTACNVGGGTWTCTAAT";
    const std::string read_seq_1f = "@M02326:75:000000000-BL25H:1:1103:15116:4927 1:N:0:CTGCCATA+CTGCCATA\nGTAGGGTGTCAGCCGCCGCGGTAATACGAAGGTGGCAAGCGTTGTTCGGATTTACTGGGCGTAAAGGGAGCGTAGGTGGTTAGGTAAGCCCTCCGGGAAATCTTCAGGCTTAACCTGAAAAGGTCGGGGGGGACTGCCTAGCTAGAGGGCGGGAGAGGAGCGCGGAATTCCCGGTGTAGCGGTGAAATGCGTAGAGATCGGGAGGAAGGCCGGTGGCGAAGGCGGCGCTCTGGAACGTTTCTGACACTGAGGCTCGGAAGCGTGGGGAGCAAACAGGATTAGAAACCCCTGTAGTCCGTAG\n+\nCCCCCGGGGGGGFGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGEGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGFGGFGGGGEGGGGGGGGGGGGGGGGGFFGGGGGGGGGGGGGGGGGGGGGGGGGGGEGCGGGGGGGGGGGGGGGGGGGGGGGGGGGGEGDGFGGGGGGGGGGGGGGGGGEGGCGGGGGGGGFGGGGGGFGGFFGFGGCFGGC?FGGGGGFFFFFFFFF:FFD?BFFFFF1?FF::D<FABB:<B4";
    const std::string read_seq_1r = "@M02326:75:000000000-BL25H:1:1103:15116:4927 2:N:0:CTGCCATA+CTGCCATA\nACGGACTACAGGGGTTTCTAATCCTGTTTGCTCCCCACGCTTCCGAGCCTCAGTGTCAGAAACGTTCCAGAGCGCCGCCTTCGCCACCGGCCTTCCTCCCGATCTCTACGCATTTCACCGCTACACCGGGAATTCCGCGCTCCTCTCCCGCCCTCTAGCTAGGCAGTCCCCCCCGACCTTTTCAGGTTAAGCCTGAAGATTTCCCGGAGGGCTTACCTAACCACCTACGCTCCCTTTACGCCCAGTAAATCCGAACAACGCTTGCCACCTTCGTATTACCCCGGCGGCTGACACCCTACAC\n+\nCCCCCGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGDGGGGGEGGGGGGGGGGGGGGGGGGGECFEGGGGDGGGGGGGGGGGGDGGGGGGGFGGGGGDGGGFGGGGGGGGGGGGGGGDFGGGGGCFFDEDEEDGGGGGF7@FFDGDDGFEFFEFGDFGGDFFGGGGEEDEGDGFFFFGFC>ECGGGFG4>?CFGF6AFGDD<5>>FF<?FGBFGFF7?FFFBB>15(8?2:<<<5@1?:?FFF2249?FB;?<(><CF>FC?FF>>328:FA24(32,497>>63*19<51:0,(";
    const std::string read_seq_2f = "@M02326:75:000000000-BL25H:1:1103:24744:3455 1:N:0:CTGCCATA+CTGCCATA\nGTCTTGTGTCAGCCGCCGCGTTAATACAGAGGTCCCAAGCGTTGTTCGGATTCATTGGGCGTAAAGGGTGCGTAGGTGGCGCCGTAAGTGGGGTGTGAAATTTCGGAGCTTAACTCCGAAACTGCATTCCATACTGCGGTGCTTGAGGACTGGAGAGGAGACCGGAATTCATGGTGTAGCAGTGAAATGCGTAGAGATCATGAGGAAGACCAGTGGCGAAGGCGGGTCTCTGGACAGTTCCTGACACTGAGGCACGAAGGCTAGGGGAGCGAACGGGATTAGATACCCTAGTAGTCCAGAA\n+\nC@BCCFDEGGGGGGGGGGGGGEGECFGGGGGGGCCGCEFGEFGGGFFEGD7FEFGGGGGGG@CCGFEEGGGGFEGGGGGDGGFG@C@FBA@FGCCF?<@<<EGGGGGGGFGGGGG<=7FGDGGEAFFFFGGGFFD<F:CFFEGCEGGGGFE:FGDFGGFGGGCFGGEG7CFGGGGGFFGFGGGFFFFFF9FEC<8FF9?F7<BCECGGFCB@FCFGC1C>8*CCCB83*+:CCCC6>?7CCFCD<476?GFFFGGDGC?@C*CFFF7;;@DEFFFFBD3:D3==>))47F1AAFA9<48<(";
    const std::string read_seq_2r = "@M02326:75:000000000-BL25H:1:1103:24744:3455 2:N:0:CTGCCATA+CTGCCATA\nGGACTACTAGGGTATCTAATCCCGTTCGCTCCCCTAGCCTTCGTGCCTCAGTGTCAGGAACTGTCCAGAGACCCGCCTTCGCCACTGGTCTTCCTCATGATCTCTACGCATTTCACTGCTACACCATGAATTCCGGTCTCCTCTCCAGTCCTCAAGCACCGCAGTATGGAATGCAGTTTCGGAGTTAAGCTCCGAAATTTCACACCCCACTTACGGCGCCACCTACGCACCCTTTACGCCCAATGAATCCGAACAACGCTTGGGCCATCTGTATCACCGCGGCGGCTGACACAGACCACAA\n+\nC@CCCFFGGFFGGGGEGGFFGEGEGGGGG,F@FGGGG?F@FEFEFG@,EF<,EFFCEGEFEFEGGGG<FGFGDG:7C,@FFCB8FGG8FF9E?FGGC9EEFCGGCFCBCEC:FGEG9<FAA;FG8?B<F9B,949+6@DGGGFFFFG8@9D8DFG9F,==1>C@EG,3==@BEA?DA9EF5E8?6,@E7,25?:6+C?CEEFF7;9DGF>)87C5/;)8C86(56=@29@@DFD8)66(77518@;;CC=1(,-4(400(0,*))0*(5:4<2-))-,-(,42,21(3(4?FF*(1(.(((";
    const std::string read_seq_3f = "@M02326:75:000000000-BL25H:1:2113:25713:12283 1:N:0:GCATCCTA+GCATCCTA\nTGGACTACTGGGGTATCTAATCCTGTTTGCTCCCCACGCTTTCGCGTCTCAGCGTCAATATCGGTCCAGGTAGCCGCCTTCGCCACCGGTGTTCCTCCTAATATCTACGGATTTCACTCCTACACTAGGAATTCCACTACCCTCTCCCGTATTCAAGTCTCCCAGTATCCAATGCACTTCCTGGGTTGAGCCCAGGGCTTTCACATCAGACTTAGAAAACCGCCTACACGCGCTTTCCGCCACATAAATCCGAACAACGCTTGCCCCCTCCGTATTACCTCGGCGGCTGCCACCGACCGGA\n+\nCCCCCFGGG@<FGGGGGGGEGFFGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGEGGGFGGGGFGGGGGGEGGFEFGGGGGFGGGGGGGGGF9FEGGGFF@ECGDGGFGGGGFGGGGFDGFFGGGGF,EFCGFDDDFDDFGFDFGFEFFGFDF9EC8@ADDGGGGFGGGFGFGGGF7D@C*>:@@,3@EDDFAFGA;@FFBFE>BCF6EBCF?E@A>D5=>-5).65?))678)(/=)).58@>0611086(/1*;+*1(0(49?(,4442)5:1)21(4(((2*/(((-2((";
    const std::string read_seq_3r = "@M02326:75:000000000-BL25H:1:2113:25713:12283 2:N:0:GCATCCTA+GCATCCTA\nGTGTCAGCCGCCGCGGTAATACGGAGGGTGCAAGCGTTGTTCGGATTTATTGGGCGTAAAGCGCGTGTAGGCGGTTTTCTAAGTCTGATGTGAAAGCCCTGGGCTCAACCCAGGAAGTGCATTGGATACTGGGAGACTTGAATACGGGAGAGGGTAGTGGAATTCCTAGTGTAGGAGTGAAATCCGTAGATATTAGGAGGAACACCGGTGGCGAAGGCGGCTACCTGGACCGATATTGACGCTGAGACGCGAAAGCGTGGGGAGCAAACAGGATTAGATACCCCAGTAGTCCAAGATCGGA\n+\nCCCCCGGGGGGDEDGGGGGGGGGG7FGGFGGGFGGGGGGGGGGGGGFGGGGFGGGGGGGGGGGGGGGGGGGGGGGGGGGGGAFFGFGGGGGGGGGGGFGGGGGGGGGGGGGGGGDEGFGGGGGGGFGGGGGGGGGDFGGGGFGGGGGGGF>FGEGGGGGFGEGGFEGGGGGGGGGC7FGGFGGGGDG7EEGGGGGGGEGGGGGFGGGGGGGGG5BBE>EGGGGGGGGFFFGEGC=EGFGFEGGGGGGF3@D=7:EFGGFGD3CC>F4::BGFFFFF?FAABF7<16F<><BAFFF?@ABF3";
    const std::string read_seq_se = "@M02326:75:000000000-BL25H:1:2113:25713:12283 1:N:0:GCATCCTA+GCATCCTA\nGTGTCAGCCGCCGCGGTAATACGGAGGGTGCAAGCGTTGTTCGGATTTATTGGGCGTAAAGCGCGTGTAGGCGGTTTTCTAAGTCTGATGTGAAAGCCCTGGGCTCAACCCAGGAAGTGCATTGGATACTGGGAGACTTGAATACGGGAGAGGGTAGTGGAATTCCTAGTGTAGGAGTGAAATCCGTAGATATTAGGAGGAACACCGGTGGCGAAGGCGGCTACCTGGACCGATATTGACGCTGAGACGCGAAAGCGTGGGGAGCAAACAGGATTAGATACCCCAGTAGTCCA\n+\nIII<IIIIIIIII<IIIIIIIIIIFIII<IIIIIIIIIIIIIIIIIIIII@>IIII?IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII";
    Primers p;
};

TEST_F(Primer, test_bounded_edit_distance) {
    ALIGNPOS r;
    // primer, seq, float_bp, max_error, end_matches
    r = p.bounded_edit_distance(primer_1, read_seq_1, 6, 4, 4);

    ASSERT_EQ(0u, r.dist);
    ASSERT_EQ(5u, r.spos);
    ASSERT_EQ(34u, r.epos);
};

TEST_F(Primer, test_pairs_both_match) {
    po::variables_map vm;

    std::istringstream in1(read_seq_1f);
    std::istringstream in2(read_seq_1r);
    InputReader<PairedEndRead, PairedEndReadFastqImpl> ifp(in1, in2);
    SeqMap primer5p = p.fasta2dict(primer_5p, "seq");
    SeqMap primer3p = p.fasta2dict(primer_3p, "seq");

    PrimerCounters counter("hts_Primers", vm);
    while(ifp.has_next()) {
        auto i = ifp.next();
        PairedEndRead* per = dynamic_cast<PairedEndRead*>(i.get());
        counter.input(*per);
        // check_read_pe(PairedEndRead &pe, PrimerCounters &counter, SeqMap &primer5p, SeqMap &primer3p, const size_t pMismatches, const size_t pEndMismatches, const size_t pfloat, const size_t flip, const size_t keep, const size_t mpmatches
        p.check_read_pe(*per, counter, primer5p, primer3p, 4, 4, 5, true, false, 2);
        counter.output(*per);
        ASSERT_EQ("TACGAAGGTGGCAAGCGTTGTTCGGATTTACTGGGCGTAAAGGGAGCGTAGGTGGTTAGGTAAGCCCTCCGGGAAATCTTCAGGCTTAACCTGAAAAGGTCGGGGGGGACTGCCTAGCTAGAGGGCGGGAGAGGAGCGCGGAATTCCCGGTGTAGCGGTGAAATGCGTAGAGATCGGGAGGAAGGCCGGTGGCGAAGGCGGCGCTCTGGAACGTTTCTGACACTGAGGCTCGGAAGCGTGGGGAGCAAACAGGATTAGAAACCCCTGTAGTCCGTAG", (per->non_const_read_one()).get_sub_seq());
        ASSERT_EQ(24u, (per->non_const_read_one()).getLTrim());
        ASSERT_EQ("CCTGTTTGCTCCCCACGCTTCCGAGCCTCAGTGTCAGAAACGTTCCAGAGCGCCGCCTTCGCCACCGGCCTTCCTCCCGATCTCTACGCATTTCACCGCTACACCGGGAATTCCGCGCTCCTCTCCCGCCCTCTAGCTAGGCAGTCCCCCCCGACCTTTTCAGGTTAAGCCTGAAGATTTCCCGGAGGGCTTACCTAACCACCTACGCTCCCTTTACGCCCAGTAAATCCGAACAACGCTTGCCACCTTCGTATTACCCCGGCGGCTGACACCCTACAC", (per->non_const_read_two()).get_sub_seq());
        ASSERT_EQ(22u, (per->non_const_read_two()).getLTrim());
        ASSERT_EQ(0u, counter.flipped);
    }
};

TEST_F(Primer, test_pairs_flipped) {
    std::istringstream in1(read_seq_3f);
    std::istringstream in2(read_seq_3r);
    InputReader<PairedEndRead, PairedEndReadFastqImpl> ifp(in1, in2);
    SeqMap primer5p = p.fasta2dict(primer_5p, "seq");
    SeqMap primer3p = p.fasta2dict(primer_3p, "seq");

    PrimerCounters counter("hts_Primers", vm);
    while(ifp.has_next()) {
        auto i = ifp.next();
        PairedEndRead* per = dynamic_cast<PairedEndRead*>(i.get());
        counter.input(*per);
        // check_read_pe(PairedEndRead &pe, PrimerCounters &counter, SeqMap &primer5p, SeqMap &primer3p, const size_t pMismatches, const size_t pEndMismatches, const size_t pfloat, const size_t flip, const size_t keep, const size_t mpmatches
        p.check_read_pe(*per, counter, primer5p, primer3p, 4, 4, 5, true, false, 2);
        counter.output(*per);
        ASSERT_EQ("TACGGAGGGTGCAAGCGTTGTTCGGATTTATTGGGCGTAAAGCGCGTGTAGGCGGTTTTCTAAGTCTGATGTGAAAGCCCTGGGCTCAACCCAGGAAGTGCATTGGATACTGGGAGACTTGAATACGGGAGAGGGTAGTGGAATTCCTAGTGTAGGAGTGAAATCCGTAGATATTAGGAGGAACACCGGTGGCGAAGGCGGCTACCTGGACCGATATTGACGCTGAGACGCGAAAGCGTGGGGAGCAAACAGGATTAGATACCCCAGTAGTCCAAGATCGGA", (per->non_const_read_one()).get_sub_seq());
        ASSERT_EQ(19u, (per->non_const_read_one()).getLTrim());
        ASSERT_EQ("CCTGTTTGCTCCCCACGCTTTCGCGTCTCAGCGTCAATATCGGTCCAGGTAGCCGCCTTCGCCACCGGTGTTCCTCCTAATATCTACGGATTTCACTCCTACACTAGGAATTCCACTACCCTCTCCCGTATTCAAGTCTCCCAGTATCCAATGCACTTCCTGGGTTGAGCCCAGGGCTTTCACATCAGACTTAGAAAACCGCCTACACGCGCTTTCCGCCACATAAATCCGAACAACGCTTGCCCCCTCCGTATTACCTCGGCGGCTGCCACCGACCGGA", (per->non_const_read_two()).get_sub_seq());
        ASSERT_EQ(21u, (per->non_const_read_two()).getLTrim());
        ASSERT_EQ(1u, counter.flipped);
    }
};

TEST_F(Primer, test_pairs_one_fail) {
    std::istringstream in1(read_seq_2f);
    std::istringstream in2(read_seq_2r);
    InputReader<PairedEndRead, PairedEndReadFastqImpl> ifp(in1, in2);
    SeqMap primer5p = p.fasta2dict(primer_5p, "seq");
    SeqMap primer3p = p.fasta2dict(primer_3p, "seq");

    PrimerCounters counter("hts_Primers", vm);
    while(ifp.has_next()) {
        auto i = ifp.next();
        PairedEndRead* per = dynamic_cast<PairedEndRead*>(i.get());
        counter.input(*per);
        // check_read_pe(PairedEndRead &pe, PrimerCounters &counter, SeqMap &primer5p, SeqMap &primer3p, const size_t pMismatches, const size_t pEndMismatches, const size_t pfloat, const size_t flip, const size_t keep, const size_t mpmatches
        p.check_read_pe(*per, counter, primer5p, primer3p, 4, 4, 5, true, false, 2);
        counter.output(*per);
        ASSERT_EQ("", (per->non_const_read_one()).get_sub_seq());
        ASSERT_EQ(301u, (per->non_const_read_one()).getRTrim());
        ASSERT_EQ("N", (per->non_const_read_two()).get_sub_seq());
        ASSERT_EQ(301u, (per->non_const_read_two()).getRTrim());
        ASSERT_EQ(0u, counter.flipped);
    }
};

TEST_F(Primer, test_pairs_one_fail_ok) {
    std::istringstream in1(read_seq_2f);
    std::istringstream in2(read_seq_2r);
    InputReader<PairedEndRead, PairedEndReadFastqImpl> ifp(in1, in2);
    SeqMap primer5p = p.fasta2dict(primer_5p, "seq");
    SeqMap primer3p = p.fasta2dict(primer_3p, "seq");

    PrimerCounters counter("hts_Primers", vm);
    while(ifp.has_next()) {
        auto i = ifp.next();
        PairedEndRead* per = dynamic_cast<PairedEndRead*>(i.get());
        counter.input(*per);
        // check_read_pe(PairedEndRead &pe, PrimerCounters &counter, SeqMap &primer5p, SeqMap &primer3p, const size_t pMismatches, const size_t pEndMismatches, const size_t pfloat, const size_t flip, const size_t keep, const size_t mpmatches
        p.check_read_pe(*per, counter, primer5p, primer3p, 4, 4, 5, true, false, 1);
        counter.output(*per);
        ASSERT_EQ("GTCTTGTGTCAGCCGCCGCGTTAATACAGAGGTCCCAAGCGTTGTTCGGATTCATTGGGCGTAAAGGGTGCGTAGGTGGCGCCGTAAGTGGGGTGTGAAATTTCGGAGCTTAACTCCGAAACTGCATTCCATACTGCGGTGCTTGAGGACTGGAGAGGAGACCGGAATTCATGGTGTAGCAGTGAAATGCGTAGAGATCATGAGGAAGACCAGTGGCGAAGGCGGGTCTCTGGACAGTTCCTGACACTGAGGCACGAAGGCTAGGGGAGCGAACGGGATTAGATACCCTAGTAGTCCAGAA", (per->non_const_read_one()).get_sub_seq());
        ASSERT_EQ(0u, (per->non_const_read_one()).getLTrim());
        ASSERT_EQ(0u, (per->non_const_read_one()).getRTrim());
        ASSERT_EQ("CCCGTTCGCTCCCCTAGCCTTCGTGCCTCAGTGTCAGGAACTGTCCAGAGACCCGCCTTCGCCACTGGTCTTCCTCATGATCTCTACGCATTTCACTGCTACACCATGAATTCCGGTCTCCTCTCCAGTCCTCAAGCACCGCAGTATGGAATGCAGTTTCGGAGTTAAGCTCCGAAATTTCACACCCCACTTACGGCGCCACCTACGCACCCTTTACGCCCAATGAATCCGAACAACGCTTGGGCCATCTGTATCACCGCGGCGGCTGACACAGACCACAA", (per->non_const_read_two()).get_sub_seq());
        ASSERT_EQ(20u, (per->non_const_read_two()).getLTrim());
        ASSERT_EQ(0u, (per->non_const_read_two()).getRTrim());
        ASSERT_EQ(0u, counter.flipped);
    }
};

TEST_F(Primer, test_pairs_flipped_5p_match_only) {
    std::istringstream in1(read_seq_3f);
    std::istringstream in2(read_seq_3r);
    InputReader<PairedEndRead, PairedEndReadFastqImpl> ifp(in1, in2);
    SeqMap primer5p = p.fasta2dict(primer_5p, "seq");
    SeqMap primer3p;

    PrimerCounters counter("hts_Primers", vm);
    while(ifp.has_next()) {
        auto i = ifp.next();
        PairedEndRead* per = dynamic_cast<PairedEndRead*>(i.get());
        counter.input(*per);
        // check_read_pe(PairedEndRead &pe, PrimerCounters &counter, SeqMap &primer5p, SeqMap &primer3p, const size_t pMismatches, const size_t pEndMismatches, const size_t pfloat, const size_t flip, const size_t keep, const size_t mpmatches
        p.check_read_pe(*per, counter, primer5p, primer3p, 4, 4, 5, true, false, 1);
        counter.output(*per);
        ASSERT_EQ("TACGGAGGGTGCAAGCGTTGTTCGGATTTATTGGGCGTAAAGCGCGTGTAGGCGGTTTTCTAAGTCTGATGTGAAAGCCCTGGGCTCAACCCAGGAAGTGCATTGGATACTGGGAGACTTGAATACGGGAGAGGGTAGTGGAATTCCTAGTGTAGGAGTGAAATCCGTAGATATTAGGAGGAACACCGGTGGCGAAGGCGGCTACCTGGACCGATATTGACGCTGAGACGCGAAAGCGTGGGGAGCAAACAGGATTAGATACCCCAGTAGTCCAAGATCGGA", (per->non_const_read_one()).get_sub_seq());
        ASSERT_EQ(19u, (per->non_const_read_one()).getLTrim());
        ASSERT_EQ("TGGACTACTGGGGTATCTAATCCTGTTTGCTCCCCACGCTTTCGCGTCTCAGCGTCAATATCGGTCCAGGTAGCCGCCTTCGCCACCGGTGTTCCTCCTAATATCTACGGATTTCACTCCTACACTAGGAATTCCACTACCCTCTCCCGTATTCAAGTCTCCCAGTATCCAATGCACTTCCTGGGTTGAGCCCAGGGCTTTCACATCAGACTTAGAAAACCGCCTACACGCGCTTTCCGCCACATAAATCCGAACAACGCTTGCCCCCTCCGTATTACCTCGGCGGCTGCCACCGACCGGA", (per->non_const_read_two()).get_sub_seq());
        ASSERT_EQ(0u, (per->non_const_read_two()).getLTrim());
        ASSERT_EQ(1u, counter.flipped);
    }
};

TEST_F(Primer, test_pairs_flipped_3p_match_only) {
    std::istringstream in1(read_seq_3f);
    std::istringstream in2(read_seq_3r);
    InputReader<PairedEndRead, PairedEndReadFastqImpl> ifp(in1, in2);
    SeqMap primer5p;
    SeqMap primer3p = p.fasta2dict(primer_3p, "seq");

    PrimerCounters counter("hts_Primers", vm);
    while(ifp.has_next()) {
        auto i = ifp.next();
        PairedEndRead* per = dynamic_cast<PairedEndRead*>(i.get());
        counter.input(*per);
        // check_read_pe(PairedEndRead &pe, PrimerCounters &counter, SeqMap &primer5p, SeqMap &primer3p, const size_t pMismatches, const size_t pEndMismatches, const size_t pfloat, const size_t flip, const size_t keep, const size_t mpmatches
        p.check_read_pe(*per, counter, primer5p, primer3p, 4, 4, 5, true, false, 1);
        counter.output(*per);
        ASSERT_EQ("GTGTCAGCCGCCGCGGTAATACGGAGGGTGCAAGCGTTGTTCGGATTTATTGGGCGTAAAGCGCGTGTAGGCGGTTTTCTAAGTCTGATGTGAAAGCCCTGGGCTCAACCCAGGAAGTGCATTGGATACTGGGAGACTTGAATACGGGAGAGGGTAGTGGAATTCCTAGTGTAGGAGTGAAATCCGTAGATATTAGGAGGAACACCGGTGGCGAAGGCGGCTACCTGGACCGATATTGACGCTGAGACGCGAAAGCGTGGGGAGCAAACAGGATTAGATACCCCAGTAGTCCAAGATCGGA", (per->non_const_read_one()).get_sub_seq());
        ASSERT_EQ(0u, (per->non_const_read_one()).getRTrim());
        ASSERT_EQ("CCTGTTTGCTCCCCACGCTTTCGCGTCTCAGCGTCAATATCGGTCCAGGTAGCCGCCTTCGCCACCGGTGTTCCTCCTAATATCTACGGATTTCACTCCTACACTAGGAATTCCACTACCCTCTCCCGTATTCAAGTCTCCCAGTATCCAATGCACTTCCTGGGTTGAGCCCAGGGCTTTCACATCAGACTTAGAAAACCGCCTACACGCGCTTTCCGCCACATAAATCCGAACAACGCTTGCCCCCTCCGTATTACCTCGGCGGCTGCCACCGACCGGA", (per->non_const_read_two()).get_sub_seq());
        ASSERT_EQ(21u, (per->non_const_read_two()).getLTrim());
        ASSERT_EQ(1u, counter.flipped);
    }
};

TEST_F(Primer, test_single_match_only) {
    std::istringstream in1(read_seq_se);
    InputReader<SingleEndRead, SingleEndReadFastqImpl> ifs(in1);
    SeqMap primer5p = p.fasta2dict(primer_5p, "seq");
    SeqMap primer3p = p.fasta2dict(primer_3p, "seq");

    PrimerCounters counter("hts_Primers", vm);
    while(ifs.has_next()) {
        auto i = ifs.next();
        SingleEndRead* ser = dynamic_cast<SingleEndRead*>(i.get());
        counter.input(*ser);
        // check_read_pe(PairedEndRead &pe, PrimerCounters &counter, SeqMap &primer5p, SeqMap &primer3p, const size_t pMismatches, const size_t pEndMismatches, const size_t pfloat, const size_t flip, const size_t keep, const size_t mpmatches
        p.check_read_se(*ser, counter, primer5p, primer3p, 4, 4, 5, true, false, 2);
        counter.output(*ser);
        ASSERT_EQ("TACGGAGGGTGCAAGCGTTGTTCGGATTTATTGGGCGTAAAGCGCGTGTAGGCGGTTTTCTAAGTCTGATGTGAAAGCCCTGGGCTCAACCCAGGAAGTGCATTGGATACTGGGAGACTTGAATACGGGAGAGGGTAGTGGAATTCCTAGTGTAGGAGTGAAATCCGTAGATATTAGGAGGAACACCGGTGGCGAAGGCGGCTACCTGGACCGATATTGACGCTGAGACGCGAAAGCGTGGGGAGCAAACAGG", (ser->non_const_read_one()).get_sub_seq());
        ASSERT_EQ(19u, (ser->non_const_read_one()).getLTrim());
        ASSERT_EQ(21u, (ser->non_const_read_one()).getRTrim());
        ASSERT_EQ(0u, counter.flipped);
    }
};

TEST_F(Primer, test_single_match_flip) {
    std::istringstream in1(read_seq_se);
    InputReader<SingleEndRead, SingleEndReadFastqImpl> ifs(in1);
    SeqMap primer3p = p.fasta2dict(primer_5p, "seq");
    SeqMap primer5p = p.fasta2dict(primer_3p, "seq");

    PrimerCounters counter("hts_Primers", vm);
    while(ifs.has_next()) {
        auto i = ifs.next();
        SingleEndRead* ser = dynamic_cast<SingleEndRead*>(i.get());
        counter.input(*ser);
        // check_read_pe(PairedEndRead &pe, PrimerCounters &counter, SeqMap &primer5p, SeqMap &primer3p, const size_t pMismatches, const size_t pEndMismatches, const size_t pfloat, const size_t flip, const size_t keep, const size_t mpmatches
        p.check_read_se(*ser, counter, primer5p, primer3p, 4, 4, 5, true, false, 2);
        counter.output(*ser);
        ASSERT_EQ("CCTGTTTGCTCCCCACGCTTTCGCGTCTCAGCGTCAATATCGGTCCAGGTAGCCGCCTTCGCCACCGGTGTTCCTCCTAATATCTACGGATTTCACTCCTACACTAGGAATTCCACTACCCTCTCCCGTATTCAAGTCTCCCAGTATCCAATGCACTTCCTGGGTTGAGCCCAGGGCTTTCACATCAGACTTAGAAAACCGCCTACACGCGCTTTACGCCCAATAAATCCGAACAACGCTTGCACCCTCCGTA", (ser->non_const_read_one()).get_sub_seq());
        ASSERT_EQ(21u, (ser->non_const_read_one()).getLTrim());
        ASSERT_EQ(19u, (ser->non_const_read_one()).getRTrim());
        ASSERT_EQ(1u, counter.flipped);
    }
};

TEST_F(Primer, test_single_match_flip_5p_only) {
    std::istringstream in1(read_seq_se);
    InputReader<SingleEndRead, SingleEndReadFastqImpl> ifs(in1);
    SeqMap primer3p;
    SeqMap primer5p = p.fasta2dict(primer_3p, "seq");

    PrimerCounters counter("hts_Primers", vm);
    while(ifs.has_next()) {
        auto i = ifs.next();
        SingleEndRead* ser = dynamic_cast<SingleEndRead*>(i.get());
        counter.input(*ser);
        // check_read_pe(PairedEndRead &pe, PrimerCounters &counter, SeqMap &primer5p, SeqMap &primer3p, const size_t pMismatches, const size_t pEndMismatches, const size_t pfloat, const size_t flip, const size_t keep, const size_t mpmatches
        p.check_read_se(*ser, counter, primer5p, primer3p, 4, 4, 5, true, false, 1);
        counter.output(*ser);
        ASSERT_EQ("CCTGTTTGCTCCCCACGCTTTCGCGTCTCAGCGTCAATATCGGTCCAGGTAGCCGCCTTCGCCACCGGTGTTCCTCCTAATATCTACGGATTTCACTCCTACACTAGGAATTCCACTACCCTCTCCCGTATTCAAGTCTCCCAGTATCCAATGCACTTCCTGGGTTGAGCCCAGGGCTTTCACATCAGACTTAGAAAACCGCCTACACGCGCTTTACGCCCAATAAATCCGAACAACGCTTGCACCCTCCGTATTACCGCGGCGGCTGACAC", (ser->non_const_read_one()).get_sub_seq());
        ASSERT_EQ(21u, (ser->non_const_read_one()).getLTrim());
        ASSERT_EQ(0u, (ser->non_const_read_one()).getRTrim());
        ASSERT_EQ(1u, counter.flipped);
    }
};

TEST_F(Primer, test_single_match_flip_3p_only) {
    std::istringstream in1(read_seq_se);
    InputReader<SingleEndRead, SingleEndReadFastqImpl> ifs(in1);
    SeqMap primer3p = p.fasta2dict(primer_5p, "seq");
    SeqMap primer5p;

    PrimerCounters counter("hts_Primers", vm);
    while(ifs.has_next()) {
        auto i = ifs.next();
        SingleEndRead* ser = dynamic_cast<SingleEndRead*>(i.get());
        counter.input(*ser);
        // check_read_pe(PairedEndRead &pe, PrimerCounters &counter, SeqMap &primer5p, SeqMap &primer3p, const size_t pMismatches, const size_t pEndMismatches, const size_t pfloat, const size_t flip, const size_t keep, const size_t mpmatches
        p.check_read_se(*ser, counter, primer5p, primer3p, 4, 4, 5, true, false, 1);
        counter.output(*ser);
        ASSERT_EQ("TGGACTACTGGGGTATCTAATCCTGTTTGCTCCCCACGCTTTCGCGTCTCAGCGTCAATATCGGTCCAGGTAGCCGCCTTCGCCACCGGTGTTCCTCCTAATATCTACGGATTTCACTCCTACACTAGGAATTCCACTACCCTCTCCCGTATTCAAGTCTCCCAGTATCCAATGCACTTCCTGGGTTGAGCCCAGGGCTTTCACATCAGACTTAGAAAACCGCCTACACGCGCTTTACGCCCAATAAATCCGAACAACGCTTGCACCCTCCGTA", (ser->non_const_read_one()).get_sub_seq());
        ASSERT_EQ(0u, (ser->non_const_read_one()).getLTrim());
        ASSERT_EQ(19u, (ser->non_const_read_one()).getRTrim());
        ASSERT_EQ(1u, counter.flipped);
    }
};
