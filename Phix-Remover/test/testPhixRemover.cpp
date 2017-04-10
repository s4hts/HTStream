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
    std::cout << "Testing Building Lookup Table with sequence " << s << '\n';
    Read readPhix = Read(s, "", ""); 
    Read testRead = Read(s, "", ""); 
    
    size_t lookup_kmer_test = 5;
    size_t true_kmer = 7; 
    firstLookup fl( new std::shared_ptr<Lookup> [1UL << (lookup_kmer_test * 2)] );
    firstLookupPointer l = fl.get();
    setLookup(l, readPhix, true_kmer, lookup_kmer_test );
    

    size_t bitKmer = true_kmer* 2;
    size_t bitKmerLookupSize = lookup_kmer_test * 2;

    size_t lookup_loc = 0; // change bit 0 and 1 the << 2
    size_t lookup_loc_rc = bitKmerLookupSize - 2; // change bit 31 and 30 then >> 2
    size_t rest_loc = 0;
    size_t rest_loc_rc = 0;

    size_t diff = 0;
    /*Particullary time consuming in each read*/
    boost::dynamic_bitset <> forwardLookup(bitKmerLookupSize);
    boost::dynamic_bitset <> reverseLookup(bitKmerLookupSize);

    boost::dynamic_bitset <> forwardRest;
    boost::dynamic_bitset <> reverseRest;

    /*If lookup is less the Rest then we have a special case*/
    if (bitKmer > bitKmerLookupSize) {
        diff = bitKmer - bitKmerLookupSize;
        forwardRest = boost::dynamic_bitset<>(  diff  );
        reverseRest = boost::dynamic_bitset<>(  diff );
        rest_loc = 0;
        rest_loc_rc = diff - 2;
    }

    

    double val = check_read(l, testRead, true_kmer, lookup_kmer_test, rest_loc, rest_loc_rc, bitKmer, bitKmerLookupSize, lookup_loc, lookup_loc_rc, diff, forwardLookup, reverseLookup, forwardRest, reverseRest );
    std::cout << "Hits should equal 4 == " << val << '\n';
    ASSERT_EQ(4, val);
};;


TEST_F(PhixRemover, setLookupTestOrderedVec) {
    std::string s("AAAAAAAGCT");
    std::cout << "Testing Building Lookup Table with sequence " << s << '\n';
    Read readPhix = Read(s, "", ""); 


    size_t lookup_kmer_test = 5;
    size_t true_kmer = 7; 
    firstLookup fl( new std::shared_ptr<Lookup> [1UL << (lookup_kmer_test * 2)] );
    firstLookupPointer l = fl.get();

    setLookup(l, readPhix, true_kmer, lookup_kmer_test );
    auto lastPosition =  l[((1 << (lookup_kmer_test * 2)) - 1)];
    auto a = lastPosition->vecBits;
    std::cout << "Should be in order\n"; 
    lastPosition->print();
    std::cout << "Test ordered Vector Words" << '\n';
    ASSERT_EQ(true, a[0] < a[1]);
    ASSERT_EQ(true, a[1] < a[2]);
    std::cout << "Test Binary Search" << '\n';
    boost::dynamic_bitset<> bs_1(4, 7);
    boost::dynamic_bitset<> bs_2(4, 9);
    boost::dynamic_bitset<> bs_3(4, 15);
    boost::dynamic_bitset<> bs_bad(4, 1);
    ASSERT_EQ(1,lastPosition->check(bs_1));
    ASSERT_EQ(1,lastPosition->check(bs_2));
    ASSERT_EQ(1,lastPosition->check(bs_3));
    ASSERT_EQ(0,lastPosition->check(bs_bad));
    //totoal locaitons
    int z = 0;
    for (int i = 0 ; i < (1 << (lookup_kmer_test * 2)) ; ++i) {
        if (l[i]) {
           ++z;
        } 
    }
    std::cout << "Number of locations in lookup should be 2 == "  << z << '\n';
    ASSERT_EQ(2, z);
};

