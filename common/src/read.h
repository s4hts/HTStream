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
    Read subread(size_t start, size_t length);
    std::string subseq(size_t start, size_t length);
};


class ReadBase {
private:
    std::string id;
    
public:
    ReadBase(const std::string& id) : id(id) {}
    virtual ~ReadBase() {};
    virtual std::string getStrKey(size_t start, size_t length) = 0;
    boost::dynamic_bitset<> strToBit(const std::string& StrKey);
    std::string getId() { return id; }

};


class PairedEndRead: public ReadBase {
private:

    Read one;
    Read two;
public:
    PairedEndRead(const Read& one, const Read& two, const std::string& id) :
        ReadBase(id), one(one), two(two) { }
    std::string getStrKey(size_t start, size_t length);

};

class SingleEndRead: public ReadBase {
private:
    Read one;
public:
    SingleEndRead(const Read& one, const std::string& id) :
        ReadBase(id), one(one) { }
    std::string getStrKey(size_t start, size_t length);

};

#endif
