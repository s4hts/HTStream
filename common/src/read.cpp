#include "read.h"
#include <boost/dynamic_bitset.hpp>


boost::dynamic_bitset<> ReadBase::strToBit(const std::string& StrKey){
  // converts a string to a 2bit representation: A:00, C:01, T:10, G:11
  boost::dynamic_bitset<> bit(2 * StrKey.length());
  size_t i = (2 * StrKey.length()) - 1;
  for(const char &c : StrKey){
    //bit <<= 2;
    switch(c) {
      case 'A': break;
      case 'C': bit[i-1] = 1;
      break;
      case 'T': bit[i] = 1;
      break;
      case 'G': bit[i] = 1; bit[i-1] = 1;
      break;
    }
    i -= 2;
  }
  return bit;
}

// Read
Read Read::subread(size_t start, size_t length){
    return Read(seq.substr(start, length), qual.substr(start,length), id);
}

std::string Read::subseq(size_t start, size_t length){
  return seq.substr(start, length);
}

//PairedEndRead
boost::dynamic_bitset<> PairedEndRead::getKey(size_t start, size_t length){
  return strToBit(one.subseq(start, length) + two.subseq(start, length));
}

//SingleEndRead
boost::dynamic_bitset<> SingleEndRead::getKey(size_t start, size_t length){
  return strToBit(one.subseq(start, 2*length));
}

size_t qual_sum(const char c, size_t s) {
    return size_t(c) + s;
}

double SingleEndRead::avg_q_score()
{
    size_t sum = std::accumulate(one.get_qual().begin(), one.get_qual().end(), 0, qual_sum);
    return sum/double(one.get_qual().length());
    
}

double PairedEndRead::avg_q_score()
{
    size_t sum = std::accumulate(one.get_qual().begin(), one.get_qual().end(), 0, qual_sum);
    sum += std::accumulate(two.get_qual().begin(), two.get_qual().end(), 0, qual_sum);
    return sum/double(one.get_qual().length() + two.get_qual().length());
}
