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

TEST_F(PhixRemover, PhixTest) {
    std::istringstream in1(readData_1);
    std::istringstream in2(readData_1);

    InputReader<PairedEndRead, PairedEndReadFastqImpl> ifp(in1, in2);

    Read readPhix = Read(phixTest, "", ""); 
    
    kmerSet lookup;
    kmerSet lookup_rc;
    setLookup(lookup, lookup_rc, readPhix, 8);

    while(ifp.has_next()) {
        auto i = ifp.next();
        PairedEndRead *per = dynamic_cast<PairedEndRead*>(i.get());
        Read &rb1 = per->non_const_read_one();
        Read rb1RC = Read(rb1.get_seq_rc(), "", "");
        std::cout << rb1.get_seq_rc() << '\n';
        size_t val = check_read(rb1, lookup, lookup_rc, 8);
        size_t val2 = check_read(rb1RC, lookup, lookup_rc, 8);
        std::cout << val << "\n";
        std::cout << val2 << "\n";

        ASSERT_EQ(5, val);
        ASSERT_EQ(5, val2);
    }
};

TEST_F(PhixRemover, setLookupTest) {
    Read readPhix = Read("AAAAAAAAAAAAAAAAG", "", ""); 
    Read read = Read("AAAAAAAAAAAAAAAAG", "", ""); 
    size_t lookup_kmer_test = 7;
    size_t true_kmer = 7; 
    firstLookup fl( new std::shared_ptr<Lookup> [1UL << (lookup_kmer_test * 2)] );
    firstLookupPointer l = fl.get();

    setLookup(l, readPhix, true_kmer, lookup_kmer_test );
    boost::dynamic_bitset<> location_equals_to(32); // (16 * 2) - (2 * 2) for the size;

    std::cout << "Done\n"; 
    l[ (1UL << (lookup_kmer_test * 2)) - 1]->print();
    exit(0);
    l[7]->print(); 
    l[(1 << (lookup_kmer_test * 2)) - 1]->print();
    std::cout << "HERE\n";

    std::cout << check_read(l, read, 16, 2) << '\n';
};
