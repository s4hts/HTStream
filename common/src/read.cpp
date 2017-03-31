#include "read.h"
#include <boost/dynamic_bitset.hpp>
#include <numeric>

std::string ReadBase::bit_to_str(const BitSet &bits) {
    size_t str_len = bits.size()/2;
    std::string out;
    out.resize(str_len);
    size_t i = bits.size() -1;
    for(size_t stridx = 0; stridx < str_len; ++stridx) {
        if (bits[i] == 0 && bits[i-1] == 0) {
            out[stridx] = 'A';
        } else if (bits[i] == 0 && bits[i-1] == 1) {
            out[stridx] = 'C';
        } else if (bits[i] == 1 && bits[i-1] == 0) {
            out[stridx] = 'G';
        } else {
            out[stridx] = 'T';
        }
        i -= 2;
    }
    return out;
}

// Read
Read Read::subread(size_t _start, size_t _length){
    return Read(seq.substr(_start, _length), qual.substr(_start,_length), id);
}

std::string Read::subseq(size_t _start, size_t _length){
    return seq.substr(_start, _length);
}

//PairedEndRead
boost::optional<BitSet> PairedEndRead::get_key(size_t start, size_t length){
    return std::max(str_to_bit(one.subseq(start, length) + two.subseq(start, length)),
                    str_to_bit(two.subseq(start, length) + one.subseq(start, length)));
}

boost::optional<BitSet> ReadBase::reverse_complement(const std::string& str, int start, int length) {
    auto rstart = str.rbegin() + start;
    auto rend = str.rbegin() + start + length;
    auto rv = str_to_bit(std::string(rstart, rend));
    if (rv) {
        return ~(*rv);
    } else {
        return rv;
    }
}

//SingleEndRead
boost::optional<BitSet> SingleEndRead::get_key(size_t _start, size_t _length){
    //The C ensures no PE and SE are mapped to the same locaitn
    return str_to_bit("C" + one.subseq(_start, _length*2));
}

inline size_t qual_sum(size_t s, const char c) {
    return size_t(c) + s;
}

double SingleEndRead::avg_q_score()
{
    size_t sum = std::accumulate(one.get_qual().begin(), one.get_qual().end(), size_t(0), qual_sum);
    return sum/double(one.get_qual().length());

}

double PairedEndRead::avg_q_score()
{
    size_t sum = std::accumulate(one.get_qual().begin(), one.get_qual().end(), size_t(0), qual_sum);
    sum += std::accumulate(two.get_qual().begin(), two.get_qual().end(), size_t(0), qual_sum);
    return sum/double(one.get_qual().length() + two.get_qual().length());
}

char Read::complement(char bp) {
    
    switch(bp) {
        case 'A':
            return 'T';
        case 'T':
            return 'A';
        case 'G':
            return 'C';
        case 'C':
            return 'G';
    }
    return 'N';
}

void PairedEndRead::setStats(Counter &c) { //Could take stats on one that is discard PE_OUt is not incremented here
    if (!one.getDiscard()) {
        c["R1_Length"] += one.getLengthTrue();
        c["R1_Left_Trim"] += one.getLTrim();
        c["R1_Right_Trim"] += one.getRTrim();
    }
    if (!two.getDiscard()) {
        c["R2_Length"] += two.getLengthTrue();
        c["R2_Left_Trim"] += two.getLTrim();
        c["R2_Right_Trim"] += two.getRTrim(); 
    }

}

void SingleEndRead::setStats(Counter &c) {
    if (!one.getDiscard()) {
        c["SE_Length"] += one.getLengthTrue();
        c["SE_Left_Trim"] += one.getLTrim();
        c["SE_Right_Trim"] += one.getRTrim();
    }
}

