#include "gtest/gtest.h"
#include <sstream>
#include <iostream>
#include "phix_remover.h"

class PhixRemover : public ::testing::Test {
    public:
        const std::string phixTest = "ACTGACTGACTGACTGACTGACTGACTG";
        const std::string readData_1 = "@R1\nAAAAACTGACTGACTGTTTT\n+\nAAAAACTGACTGACTGTTTT\n";
        const size_t lookup_kmer_test = 2;
};

TEST_F(PhixRemover, check_check_read) {
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


TEST_F(PhixRemover, setLookupTestOrderedVec) {
    std::string s("AAAAAAAGCT");
    std::cout << "Testing Building Lookup Table with sequence " << s << '\n';
    Read readPhix = Read(s, "", ""); 
    kmerSet lookup;
    setLookup(lookup, readPhix, 5);
    ASSERT_EQ(4, lookup.size());
};

