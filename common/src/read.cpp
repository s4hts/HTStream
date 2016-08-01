#include "read.h"
#include <boost/dynamic_bitset.hpp>
#include <numeric>

boost::dynamic_bitset<> ReadBase::str_to_bit(const std::string& StrKey){
  // converts a string to a 2bit representation: A:00, T:11, C:01, G:10
  // ~ will then convert to the complimentary bp
  boost::dynamic_bitset<> bit(2 * StrKey.length());
  size_t i = (2 * StrKey.length()) - 1;

  for(const char &c : StrKey){
      //  non branching conversion
      bit[i-1] = !!((c + 10) & ( 1 << 2));
      bit[i] = !!((c + 10) & ( 1 << 4));
      i -= 2;
  }
  return bit;
}

std::string ReadBase::bit_to_str(const boost::dynamic_bitset<> &bits) {
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
Read Read::subread(size_t start, size_t length){
    return Read(seq.substr(start, length), qual.substr(start,length), id);
}

std::string Read::subseq(size_t start, size_t length){
    return seq.substr(start, length);
}

//PairedEndRead
boost::dynamic_bitset<> PairedEndRead::get_key(size_t start, size_t length){
    return str_to_bit(one.subseq(start, length) + two.subseq(start, length));
}

BitSet ReadBase::reverse_complement(const std::string& str, int start, int length) {
    auto rstart = str.rbegin() + start;
    auto rend = str.rbegin() + start + length;
    return ~str_to_bit(std::string(rstart, rend));
}

//SingleEndRead
boost::dynamic_bitset<> SingleEndRead::get_key(size_t start, size_t length){
    
    return std::max (str_to_bit(one.subseq(start, length)), reverse_complement(one.get_seq(), int(start), int(length)));
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
