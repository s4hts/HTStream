#include "read.h"
#include <boost/dynamic_bitset.hpp>


Read Read::subread(int start, int length){
    //if(seq.length >= start + length){
    // Test for a throw an exception or just let string through an out_of_range exception?
    return Read(seq.substr(start, length), qual.substr(start,length));
    //}else{
    //}
}

std::string subseq(int start, int length){
  return seq.substr(start, length);
}

std::string PairedEndRead::getStrKey(int start, int length){
  return one.subseq(start, length) + two.subseq(start, length);
}

std::string SingleEndRead::getStrKey(int start, int length){
  return one.subseq(start, 2*length);
}
