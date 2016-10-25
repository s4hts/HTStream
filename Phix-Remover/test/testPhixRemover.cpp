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
    const size_t kmer = 8;
    Read readPhix = Read(phixTest, "", ""); 
    
    std::array<size_t, 1<<kmer*2> lookup;
    std::array<size_t, 1<<kmer*2> lookup_rc;
    setLookup(lookup, lookup_rc, readPhix, kmer);

    while(ifp.has_next()) {
        auto i = ifp.next();
        PairedEndRead *per = dynamic_cast<PairedEndRead*>(i.get());
        Read &rb1 = per->non_const_read_one();
        Read rb1RC = Read(rb1.get_seq_rc(), "", "");
        std::cout << rb1.get_seq_rc() << '\n';
        size_t val = check_read(rb1, lookup, lookup_rc, 8);
        size_t val2 = check_read(rb1RC, lookup, lookup_rc, 8);
        ASSERT_EQ(5, val);
        ASSERT_EQ(5, val2);
    }
};

TEST_F(PhixRemover, LookupTable) {
    
    const size_t kmer = 8;
    Read readPhix = Read("AAAAAAAACTG", "", ""); 
                        //0123456789
    std::array<size_t, 1<<kmer*2> lookup;
    std::array<size_t, 1<<kmer*2> lookup_rc;
    setLookup(lookup, lookup_rc, readPhix, kmer);

    std::cout << lookup[1] << '\n';
    std::cout << lookup[7] << '\n';
    std::cout << lookup[30] << '\n';
    
    std::cout << lookup_rc[12287] << '\n';
    std::cout << lookup_rc[49151] << '\n';
    std::cout << lookup_rc[65535] << '\n';
    
    
    
    ASSERT_EQ(1,  lookup[0]  );
    ASSERT_EQ(1,  lookup[1]  );
    ASSERT_EQ(1,  lookup[7]  );
    ASSERT_EQ(1,  lookup[30]  );
    ASSERT_EQ(1,  lookup_rc[19455]  );
    ASSERT_EQ(1,  lookup_rc[12287]  );
    ASSERT_EQ(1,  lookup_rc[49151]  );
    ASSERT_EQ(1,  lookup_rc[65535]  );

}
