#include "read.h"
#include <boost/dynamic_bitset.hpp>


// ReadBase
//ReadBase::~ReadBase(){}

boost::dynamic_bitset<> ReadBase::strToBit(std::string& StrKey){
  boost::dynamic_bitset<> bit(2 * StrKey.length());
  //std::string seq("ACTG");

  for(char &c : StrKey){
    bit <<= 2;
    switch(c) {
      case 'A': break;
      case 'C': bit[0] = 1;
      break;
      case 'T': bit[1] = 1;
      break;
      case 'G': bit[0] = 1; bit[1] = 1;
      break;
    }
  }
  return bit;

}


// Read
Read Read::subread(size_t start, size_t length){
    //if(seq.length >= start + length){
    // Test for a throw an exception or just let string through an out_of_range exception?
    return Read(seq.substr(start, length), qual.substr(start,length));
    //}else{
    //}
}

std::string Read::subseq(size_t start, size_t length){
  return seq.substr(start, length);
}

//PairedEndRead
std::string PairedEndRead::getStrKey(size_t start, size_t length){
  return one.subseq(start, length) + two.subseq(start, length);
}


//SingleEndRead
std::string SingleEndRead::getStrKey(size_t start, size_t length){
  return one.subseq(start, 2*length);
}
