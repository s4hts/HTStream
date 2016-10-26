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
    
    kmerArray lookup;
    kmerArray lookup_rc;
    setLookup(lookup, lookup_rc, readPhix);

    while(ifp.has_next()) {
        auto i = ifp.next();
        PairedEndRead *per = dynamic_cast<PairedEndRead*>(i.get());
        Read &rb1 = per->non_const_read_one();
        Read rb1RC = Read(rb1.get_seq_rc(), "", "");
        std::cout << rb1.get_seq_rc() << '\n';
        size_t val = check_read(rb1, lookup, lookup_rc);
        size_t val2 = check_read(rb1RC, lookup, lookup_rc);
        ASSERT_EQ(5, val);
        ASSERT_EQ(5, val2);
    }
};

TEST_F(PhixRemover, LookupTable) {
    
    const size_t kmer = 8;
    Read readPhix = Read("AAAAAAAACTG", "", ""); 

    kmerArray lookup;
    kmerArray lookup_rc;
    setLookup(lookup, lookup_rc, readPhix);

    ASSERT_EQ(1,  lookup[0]  );
    ASSERT_EQ(1,  lookup[1]  );
    ASSERT_EQ(1,  lookup[7]  );
    ASSERT_EQ(1,  lookup[30]  );
    ASSERT_EQ(1,  lookup_rc[19455]  );
    ASSERT_EQ(1,  lookup_rc[12287]  );
    ASSERT_EQ(1,  lookup_rc[49151]  );
    ASSERT_EQ(1,  lookup_rc[65535]  );

}
