#include "read.h"
#include <boost/dynamic_bitset.hpp>
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
    size_t str_len = bits.size();
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

boost::optional<BitSet> ReadBase::bitjoin(const boost::optional<BitSet> &bit1, const boost::optional<BitSet> &bit2 ) {
    size_t bits = bit1 -> size() + bit2 -> size();
    BitSet bittag(bits);
    for(size_t i = 0; i < bit1 -> size(); i++) { bittag[i] = (int)(*bit1)[i]; }
    for(size_t i = 0; i < bit2 -> size(); i++) { bittag[i + (bit1 -> size())] = (int)(*bit2)[i]; }
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
      return std::max(str_to_bit(one->subseq(_start, _length) + two->subseq(_start, _length)),
                      str_to_bit(two->subseq(_start, _length) + one->subseq(_start, _length)));
    }
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
    //The C ensures no PE and SE are mapped to the same location
    if (one->getLength() <= (_start+_length*2)){
      return boost::none;
    } else {
      return str_to_bit("C" + one->subseq(_start, _length*2));
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
