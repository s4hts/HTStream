#include "minimizer.h"
#include <gtest/gtest.h>

#include <iostream>

class MinimizerTest: public ::testing::Test {
public:
    const std::string testseq = "GTCATGCACGTTCAC";
};

TEST_F(MinimizerTest, basic_minimizers) {
    Minimizer<ReferencePtr> min(4, 3);

    std::istringstream fa_to_read(string2fasta(testseq, "seq"));
    InputReader<Reference, FastaReadImpl> faReader(fa_to_read);
    auto refu = faReader.next();

    //ReadPtr = std::make_shared<Read>(testseq, "##", "");
    ReferencePtr ref = ReferencePtr(std::move(refu));
    min.find_mins(ref);
    
    ASSERT_EQ(min.get_kmers().size(), 3u);
    auto i = min.get_kmers().find(*min.two_bit("AAC"));
    ASSERT_EQ(i->second.ref->get_seq().substr(i->second.offset, i->second.length), "CACGTTCAC");
    i = min.get_kmers().find(*min.two_bit("ACG"));
    ASSERT_EQ(i->second.ref->get_seq().substr(i->second.offset, i->second.length), "TGCACGT");
    i = min.get_kmers().find(*min.two_bit("ATG"));
    ASSERT_EQ(i->second.ref->get_seq().substr(i->second.offset, i->second.length), "GTCATGCAC");

    // for(auto i = min.kmers.begin(); i != min.kmers.end(); ++i) {
    //     std::cout << i->first << " : " << i->second.read->get_seq().substr(i->second.offset, i->second.length) << std::endl;
    // }

    Minimizer<ReferencePtr> min2(3, 3);
    min2.find_mins(ref);
    for(auto i = min2.get_kmers().begin(); i != min2.get_kmers().end(); ++i) {
        std::cout << i->first << " : " << i->second.ref->get_seq().substr(i->second.offset, i->second.length) << std::endl;
    }

}
