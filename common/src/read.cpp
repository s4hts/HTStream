#include "read.h"
#include <boost/bind.hpp>
#include <numeric>


std::string strjoin(const std::vector <std::string>& v, const std::string& delim) {
    std::ostringstream s;
    for (const auto& i : v) {
        if (&i != &v[0]) {
            s << delim;
        }
        s << i;
    }
    return s.str();
}

std::string ReadBase::bit_to_str(const BitSet &bits) {
    return bits.to_string();
}

boost::optional<BitSet> ReadBase::bitjoin(const boost::optional<BitSet> &bit1, const boost::optional<BitSet> &bit2) {
    if ((bit1 == boost::none) || (bit2 == boost::none)) {
        return boost::none;
    }
    BitSet bittag(*bit1);
    bittag.append(*bit2);
    return bittag;
}

// Read
Read Read::subread(size_t _start, size_t _length){
    return Read(seq.substr(_start, _length), qual.substr(_start,_length), id);
}

std::string Read::subseq(size_t _start, size_t _length){
    return seq.substr(_start, _length);
}

//PairedEndRead
boost::optional<BitSet> PairedEndRead::get_key(size_t _start, size_t _length){
    if (std::min(one->getLength(), two->getLength()) <= _start+_length){
      return boost::none;
    } else {
      auto key1 = BitSet::from_subsequence(one->get_seq(), _start, _length);
      auto key2 = BitSet::from_subsequence(two->get_seq(), _start, _length);
      if (!key1 || !key2) {
          return boost::none;
      }
      BitSet one_two(*key1);
      one_two.append(*key2);
      BitSet two_one(*key2);
      two_one.append(*key1);
      return std::max(one_two, two_one);
    }
}

boost::optional<BitSet> ReadBase::reverse_complement(const std::string& str, int start, int length) {
    BitSet out;
    auto rstart = str.rbegin() + start;
    auto rend = str.rbegin() + start + length;
    for (auto bp = rstart; bp != rend; ++bp) {
        switch (*bp) {
        case 'A':
        case 'a':
            out.append_code(3);
            break;
        case 'C':
        case 'c':
            out.append_code(2);
            break;
        case 'G':
        case 'g':
            out.append_code(1);
            break;
        case 'T':
        case 't':
            out.append_code(0);
            break;
        default:
            return boost::none;
        }
    }
    return out;
}

//SingleEndRead
boost::optional<BitSet> SingleEndRead::get_key(size_t _start, size_t _length){
    //The C ensures no PE and SE are mapped to the same location
    if (one->getLength() <= (_start+_length*2)){
      return boost::none;
    } else {
      auto key = BitSet::from_subsequence(one->get_seq(), _start, _length * 2);
      if (!key) {
          return boost::none;
      }
      BitSet se_key;
      se_key.append_code(1);
      se_key.append(*key);
      return se_key;
    }
}

inline double qual_sum(double s, const char c, size_t offset) {
    return (double(c) - offset) + s; //need ascii offset
}

double SingleEndRead::avg_q_score(const size_t qual_offset)
{
    double sum = std::accumulate(one->get_qual().begin(), one->get_qual().end(), double(0), boost::bind(&qual_sum, _1, _2, qual_offset));
    return sum/double(one->get_qual().length());

}

double PairedEndRead::avg_q_score(const size_t qual_offset)
{
    double sum = std::accumulate(one->get_qual().begin(), one->get_qual().end(), double(0), boost::bind(&qual_sum, _1, _2, qual_offset));
    sum += std::accumulate(two->get_qual().begin(), two->get_qual().end(), double(0), boost::bind(&qual_sum, _1, _2, qual_offset));
    return sum/double(one->get_qual().length() + two->get_qual().length());
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

void SingleEndRead::accept(ReadVisitor &rv) {
    rv.visit(this);
}

void PairedEndRead::accept(ReadVisitor &rv) {
    rv.visit(this);
}
