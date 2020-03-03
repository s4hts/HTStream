#ifndef READ_H
#define READ_H

#include <boost/dynamic_bitset.hpp>
#include <boost/optional.hpp>
#include <memory>
#include <iostream>
#include <regex>
#include <unordered_map>
#include "typedefs.h"

typedef boost::dynamic_bitset<> BitSet;

class ReadBase {
public:
    virtual ~ReadBase() {}
    ReadBase(ReadBase const&) = default;
    ReadBase() = default;
    virtual boost::optional<boost::dynamic_bitset<>> get_key(size_t start, size_t length) = 0;
    static boost::optional<BitSet> str_to_bit(const std::string& StrKey) {
          // converts a string to a 2bit representation: A:00, T:11, C:01, G:10
        // ~ will then convert to the complimentary bp
        BitSet bit(2 * StrKey.length());
        size_t i = (2 * StrKey.length()) -1;
        for (const char &c : StrKey) {
            switch(c) {
            case 'A':
                break;
            case 'C':
                bit[i-1] = 1;
                break;
            case 'G':
                bit[i] = 1;
                break;
            case 'T':
                bit[i] = 1;
                bit[i-1] = 1;
                break;
            case 'N':
                return boost::none;
                break;
            }
            i -= 2;
        }
        return bit;
    }
    static std::string bit_to_str(const BitSet &bits);
    static boost::optional<BitSet> reverse_complement(const std::string& str, int start, int length);
    virtual double avg_q_score() = 0;
    bool rc;
};
typedef std::shared_ptr<ReadBase> ReadBasePtr;

class Read {
private:
    std::string seq;
    std::string qual;
    std::string id_orig;
    std::string id;
    std::string id2;
    std::vector<std::string> comments;
    size_t length;
    size_t cut_R;
    size_t cut_L;
    bool discard;
    size_t minLength;
public:
    Read(const std::string& seq_, const std::string& qual_, const std::string& id_) :
        seq(seq_), qual(qual_), id_orig(id_), length(seq_.length()), cut_R(seq_.length()), cut_L(0), discard(false), minLength(1) {
           std::regex re("\\s{1,}");
           std::string fmt = " ";
           std::string s = regex_replace(id_orig, re, fmt);
           std::string fastq_delimiter = " ";
           // split the id into 2 delimited on the first space
           size_t pos = 0;
           pos = s.find(fastq_delimiter);
           if (pos != std::string::npos){
             id = s.substr(0, pos);
             id2 = s.erase(0, pos + fastq_delimiter.length());
           } else {
             id = s;
             id2 = "";
           }
         }
    Read() :
        seq(""), qual(""), id_orig(""), id(""), id2(""), length(0), cut_R(seq.length()), cut_L(0), discard(false), minLength(1) {
        }
    Read subread(size_t _start, size_t _length);
    std::string subseq(size_t _start, size_t _length);
    const std::string& get_seq() const  { return seq; }
    const std::string& get_qual() const { return qual; }
    const std::string get_id() const { std::string tmp = (id2 == "") ? id : id + get_comment() + " " + id2; return tmp; }
    const std::string get_id_first() const { return id; }
    const std::string get_comment(bool sam=false) const {
      std::string result = "";
      const char delim = (sam) ? '\t' : '_';
      for (auto const& s : comments) { result = result + delim + s; }
      return result;
    }
    static char complement(char bp);

    const std::string get_sub_seq() const { return cut_R < cut_L ? "" : seq.substr(cut_L, cut_R - cut_L); }
    const std::string get_sub_qual() const { return cut_R < cut_L ? "" : qual.substr(cut_L, cut_R - cut_L); }


    const std::string get_seq_rc() const { if (cut_R < cut_L) { return ""; }
                                           std::string s = seq.substr(cut_L, cut_R - cut_L) ;
                                           std::transform(begin(s), end(s), begin(s), complement);
                                           std::reverse(begin(s), end(s));
                                           return s; }


    const std::string get_qual_rc() const { if (cut_R < cut_L) { return ""; }
                                            std::string q = qual.substr(cut_L, cut_R - cut_L);
                                            std::reverse(begin(q), end(q));
                                            return q;  }


    void add_comment( std::string tag ) { comments.push_back(tag);}
    void set_read_rc() {
        if (cut_R < cut_L) {
            return;
        }
        std::string s = seq.substr(cut_L, cut_R - cut_L) ;
        std::transform(begin(s), end(s), begin(s), complement);
        std::reverse(begin(s), end(s));
        std::string q = qual.substr(cut_L, cut_R - cut_L);
        std::reverse(begin(q), end(q));
        seq = s;
        qual = q;
    }

    void changeSeq( size_t loc, char bp ) { seq[loc] = bp; }
    void changeQual( size_t loc, char score ) {qual[loc] = score; }

    void setRCut( size_t cut_R_ ) { cut_R = cut_R_;}
    void setLCut( size_t cut_L_ ) { cut_L = cut_L_; }
    bool getDiscard() { discard = int(minLength) > int(cut_R) - int(cut_L); return discard; }
    void setDiscard(size_t minLength_) { minLength = minLength_; }
    size_t getLength() const { return length; }
    size_t getLengthTrue() { return cut_R < cut_L ? 0 : cut_R - cut_L; }
    size_t getLTrim() { return cut_L; }
    size_t getRTrim() { return length - cut_R; }
};

class PairedEndRead: public ReadBase {
private:
    Read one;
    Read two;
public:
    PairedEndRead(const Read& one_, const Read& two_) :
        one(one_), two(two_) { }
    boost::optional<BitSet> get_key(size_t start, size_t length);
    Read& non_const_read_one() { return one; }
    Read& non_const_read_two() { return two; }
    void checkDiscarded(size_t minLength) {one.setDiscard(minLength); two.setDiscard(minLength);}
    const Read& get_read_one() const { return one; }
    const Read& get_read_two() const { return two; }
    double avg_q_score();

    std::shared_ptr<ReadBase> convert(bool stranded);
};
typedef std::shared_ptr<PairedEndRead> PairedEndReadPtr;

/* start, finish, discarded */
class SingleEndRead: public ReadBase {
private:
    Read one;
public:
    SingleEndRead(const Read& one_) :
        one(one_) { }
    boost::optional<BitSet> get_key(size_t start, size_t length);
    Read& non_const_read_one() { return one; }
    const Read& get_read() const { return one; }
    void checkDiscarded(size_t minLength) {one.setDiscard(minLength);}
    double avg_q_score();
    std::shared_ptr<ReadBase> convert(bool stranded);
    void set_read_rc() { one.set_read_rc();}
};
typedef std::shared_ptr<SingleEndRead> SingleEndReadPtr;

#endif
