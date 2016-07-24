#ifndef READ_H
#define READ_H

#include <boost/dynamic_bitset.hpp>

class Read {
private:
    std::string seq;
    std::string qual;

public:
    Read(const std::string& seq_, const std::string& qual_) :
        seq(seq_), qual(qual_) { }
    Read subread(int start, int length);
    std::string subseq(int start, int length);
};


class ReadBase {
public:
    //virtual ~ReadBase();
    virtual std::string getStrKey(int start, int length) = 0;
    boost::dynamic_bitset<> strToBit(std::string& StrKey);
};


class PairedEndRead: public ReadBase {
private:

    Read one;
    Read two;
    std::string id;
public:
  //~PairedEndRead();
    PairedEndRead(const Read& one, const Read& two, const std::string& id) :
        one(one), two(two), id(id) { }
    virtual std::string getStrKey(int start, int length);

};

class SingleEndRead: public ReadBase {
private:
    Read one;
    std::string id;

public:
    SingleEndRead(const Read& one, const std::string& id) :
        one(one), id(id) { }
    virtual std::string getStrKey(int start, int length);

};

#endif
