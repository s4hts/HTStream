#ifndef READ_H
#define READ_H

#include <boost/dynamic_bitset.hpp>

class Read {
private:
    std::string seq;
    std::string qual;
    std::string id;

public:
    Read(const std::string& seq_, const std::string& qual_, const std::string& id_) :
        seq(seq_), qual(qual_), id(id_) { }
    Read subread(size_t start, size_t length);
    std::string subseq(size_t start, size_t length);
    const std::string& get_seq() const { return seq; }
    const std::string& get_qual() const { return qual; }
    const std::string& get_id() const { return id; }
};


class ReadBase {
public:
    virtual ~ReadBase() {}
    virtual boost::dynamic_bitset<> getKey(size_t start, size_t length) = 0;
    static boost::dynamic_bitset<> strToBit(const std::string& StrKey);
    virtual double avg_q_score() = 0;
};


class PairedEndRead: public ReadBase {
private:
    Read one;
    Read two;
public:
    PairedEndRead(const Read& one_, const Read& two_) :
        one(one_), two(two_) { }
    boost::dynamic_bitset<> getKey(size_t start, size_t length);
    const Read& get_read_one() const { return one; }
    const Read& get_read_two() const { return two; }
    double avg_q_score();
};

class SingleEndRead: public ReadBase {
private:
    Read one;
public:
    SingleEndRead(const Read& one_) :
        one(one_) { }
    boost::dynamic_bitset<> getKey(size_t start, size_t length);
    const Read& get_read() const { return one; }
    double avg_q_score();
};

#endif
