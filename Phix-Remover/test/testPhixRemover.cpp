#include "gtest/gtest.h"
#include <sstream>
#include <iostream>
#include "phix_remover.h"

class PhixRemover : public ::testing::Test {
    public:
        const std::string phixTest = "ACTGACTGACTGACTGACTGACTGACTG";
        const std::string readData_1 = "@R1\nAAAAACTGACTGACTGTTTT\n+\nAAAAACTGACTGACTGTTTT\n";
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

TEST_F(PhixRemover, LookupTable) {
    
    const size_t kmer = 8;
    Read readPhix = Read("AAAAAAAACTG", "", ""); 

    kmerSet lookup;
    kmerSet lookup_rc;
    setLookup(lookup, lookup_rc, readPhix, 8);
    boost::dynamic_bitset<> test(16);
    for (auto it = lookup.begin(); it != lookup.end(); ++it) {
        std::cout << *it << '\n';
    }
    for (auto it = lookup_rc.begin(); it != lookup_rc.end(); ++it) {
        std::cout << *it << '\n';
    }
    int lookupTest =   lookup.find( test ) != lookup.end(); // ALL A's
    ASSERT_EQ(1, lookupTest );

    test.set(0);
    lookupTest =   lookup.find( test ) != lookup.end(); // C added to the end
    ASSERT_EQ(1,  lookupTest  );
    
    test <<= 2;
    test.set(0);
    test.set(1);
    lookupTest = lookup.find( test ) != lookup.end(); // T added to the end
    ASSERT_EQ(1,  lookupTest  );
    
    test <<= 2;
    test.set(1);
    lookupTest = lookup.find( test ) != lookup.end(); // T added to the end
    ASSERT_EQ(1,  lookupTest  );

    test.reset(); // test rc
    test.flip();  //start with ALL T's
    lookupTest = lookup_rc.find( test ) != lookup_rc.end(); // T added to the end

    ASSERT_EQ(1,  lookupTest  );
    
    test.flip(14); // flips 14th pos 1011111... For "G" but in string "C"
    lookupTest = lookup_rc.find( test ) != lookup_rc.end(); 
    ASSERT_EQ(1,  lookupTest  );
    
    test >>= 2; // 0 in 14 and 15 for T
    lookupTest = lookup_rc.find( test ) != lookup_rc.end(); 
    
    ASSERT_EQ(1,  lookupTest );
    
    test >>= 2; // 1 at 14 for "C" but acutally "G"
    test.set(14);
    lookupTest = lookup_rc.find( test ) != lookup_rc.end(); 

    ASSERT_EQ(1,  lookupTest );

}
