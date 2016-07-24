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
    virtual std::string getStrKey(int start, int length) = 0;
    virtual ~ReadBase();

    boost::dynamic_bitset strToBit(const std::string& key);
};

class PairedEndRead: public ReadBase {
private:

    Read one;
    Read two;
    std::string id;
public:
    PairedEndRead(const read& one, const read& two, const std::string& id) :
        one(one), two(two), id(id) { }

};

class SingleEndRead: public ReadBase {
private:
    Read one;
    std::string id;

public:
    SingleEndRead(const read& one, const std::string& id) :
        one(one), id(id) { }

};

#endif
