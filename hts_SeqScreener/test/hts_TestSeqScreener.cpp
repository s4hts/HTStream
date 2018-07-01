#include "gtest/gtest.h"
#include <sstream>
#include <iostream>
#include "hts_SeqScreener.h"

class SeqScreener : public ::testing::Test {
    public:
        const std::string phixTest = "ACTGACTGACTGACTGACTGACTGACTG";
        const std::string readData_1 = "@R1\nAAAAACTGACTGACTGTTTT\n+\nAAAAACTGACTGACTGTTTT\n";
        const size_t lookup_kmer_test = 2;
};

TEST_F(SeqScreener, all_from_fastq) {
    const std::string faFile = ">1\nACGT\nACGT\n>2\nTTTT\n";
    std::istringstream fa(faFile);
    InputReader<SingleEndRead, FastaReadImpl> f(fa);
    Read r = fasta_set_to_one_read( f );
    std::cout << r.get_seq() << '\n';
    ASSERT_EQ("ACGTACGTTTTT", r.get_seq());
};

TEST_F(SeqScreener, check_check_read) {
    std::string s("AAAAAAAGCT");
    Read readPhix = Read(s, "", "");
    Read testRead = Read(s, "", "");
    kmerSet lookup;
    size_t true_kmer = 5;
    setLookup(lookup, readPhix, true_kmer);

    boost::dynamic_bitset<> fLu(true_kmer * 2);
    boost::dynamic_bitset<> rLu(true_kmer * 2);
    size_t lookup_loc = 0;
    size_t lookup_loc_rc = (true_kmer * 2) -2;
    double val = check_read(lookup, testRead, true_kmer * 2, lookup_loc, lookup_loc_rc, fLu, rLu );
    std::cout << "Hits should equal 6 == " << val << '\n';
    ASSERT_EQ(6, val);
};;


TEST_F(SeqScreener, setLookupTestOrderedVec) {
    std::string s("AAAAAAAGCT");
    std::cout << "Testing Building Lookup Table with sequence " << s << '\n';
    Read readPhix = Read(s, "", "");
    kmerSet lookup;

    setLookup(lookup, readPhix, 5);
    std::cout << lookup.size() << '\n';
    ASSERT_EQ(4, lookup.size());
};

TEST_F(SeqScreener, setLookupTest) {
    std::string s("AAAAAAAGCT");
    std::cout << "Testing Building Lookup Table with sequence " << s << '\n';
    Read readPhix = Read(s, "", "");
    kmerSet lookup;
    setLookup(lookup, readPhix, 5);
    //ASSERT_EQ(true , lookup.end != lookup.find(boost::dynamic_bitset<>(10, "1001111111")))
    //std::string s("1001111111");

};

TEST_F(SeqScreener, setLookupTestAmbiguities) {
    std::string s("GAAAAAAGCVM");
    std::cout << "Testing Building Lookup Table with sequence ambiguities " << s << '\n';
    Read readAmbiguities = Read(s, "", "");
    std::cout << "Seq in read " << readAmbiguities.get_seq() << '\n';
    kmerSet lookup;
    setLookup(lookup, readAmbiguities, 5);
    std::cout << "myset contains:";
    for ( auto it = lookup.begin(); it != lookup.end(); ++it )
        std::cout << " " << *it;
    std::cout << std::endl;
    std::cout << lookup.size() << '\n';
    ASSERT_EQ(4, lookup.size());
};

TEST_F(SeqScreener, setLookupTestAmbiguities2_nokmers) {
    std::string s("GAANAAGCVM");
    std::cout << "Testing Building Lookup Table with sequence ambiguities " << s << '\n';
    Read readAmbiguities = Read(s, "", "");
    std::cout << "Seq in read " << readAmbiguities.get_seq() << '\n';
    kmerSet lookup;
    setLookup(lookup, readAmbiguities, 5);
    std::cout << lookup.size() << '\n';
    ASSERT_EQ(0, lookup.size());
};
