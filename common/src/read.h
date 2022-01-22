#ifndef READ_H
#define READ_H
#define NDEBUG

#include <boost/dynamic_bitset.hpp>
#include <boost/optional.hpp>
#include <boost/algorithm/string.hpp>
#include <memory>
#include <iostream>
#include <sstream>
#include <unordered_map>
#include "typedefs.h"

typedef boost::dynamic_bitset<> BitSet;
std::string strjoin(const std::vector <std::string>& v, const std::string& delim);

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

public:
    Read(const std::string& seq_, const std::string& qual_, const std::string& id_) :
        seq(seq_), qual(qual_), id_orig(id_), length(seq_.length()), cut_R(seq_.length()), cut_L(0), discard(false) {
           std::string fastq_delimiter =  " \t\f\n\r\v";
           // split the id into 2 delimited on the first space
           size_t pos = 0;
           pos = id_orig.find_first_of(fastq_delimiter);
           if (pos != std::string::npos && pos+2 < id_orig.size()){
             id = id_orig.substr(0, pos);
             // id2 remove the read designation, add back in on write out
             id2 = id_orig.substr(pos+2);
           } else {
             id = id_orig;
             id2 = "";
           }
           // Identifies any comments in the read id (id) and puts them into the comment vector
           std::vector <std::string> comment_tmp;
           boost::split(comment_tmp, id, boost::is_any_of("|"));
           id = comment_tmp[0];
           join_comment(comment_tmp, 1);
         }
    Read() :
        seq(""), qual(""), id_orig(""), id(""), id2(""), length(0), cut_R(seq.length()), cut_L(0), discard(false) {
        }
    Read subread(size_t _start, size_t _length);
    std::string subseq(size_t _start, size_t _length);
    const std::string& get_seq() const  { return seq; }
    const std::string& get_qual() const { return qual; }
    const std::string get_id_orig() const { return id_orig; }
    const std::string get_id_fastq(const std::string& read="", const std::string& sam_comment="") const {
        std::string tmp = id + sam_comment;
        if (!(id2 == "")) tmp = tmp + ' ' + read + id2;
        return tmp;
    }
    const std::string get_id_tab(const std::string& read="") const {
        std::string tmp = id;
        if (!(id2 == "")) tmp = tmp + ' ' + read + id2;
        return tmp;
    }
    const std::string get_id_first() const { return id; }
    const std::string get_id_second() const { return id2; }
    std::vector<std::string> get_comment() const { return comments;}
    static char complement(char bp);

    const std::string get_sub_seq() const { return cut_R <= cut_L ? "N" : seq.substr(cut_L, cut_R - cut_L); }
    const std::string get_sub_qual() const { return cut_R <= cut_L ? "#" : qual.substr(cut_L, cut_R - cut_L); }


    const std::string get_seq_rc() const { if (cut_R <= cut_L) { return "N"; }
                                           std::string s = seq.substr(cut_L, cut_R - cut_L) ;
                                           std::transform(begin(s), end(s), begin(s), complement);
                                           std::reverse(begin(s), end(s));
                                           return s; }


    const std::string get_qual_rc() const { if (cut_R <= cut_L) { return "#"; }
                                            std::string q = qual.substr(cut_L, cut_R - cut_L);
                                            std::reverse(begin(q), end(q));
                                            return q;  }


    void add_comment( const std::string& tag ) { if (tag != "") comments.push_back(tag);}
    void join_comment( const std::vector<std::string>& new_comments, size_t offset = 0) { comments.insert(comments.end(), new_comments.begin() + offset, new_comments.end()); }
    void set_read_rc() {
        if (cut_R <= cut_L) {
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

    void setRCut( size_t cut_R_ ) { assert(cut_R_ <= length); cut_R = cut_R_; }
    void setLCut( size_t cut_L_ ) { assert(cut_L_ <= length); cut_L = cut_L_; }
    bool getDiscard() const { return discard; }
    void setDiscard() { discard = true; }
    size_t getLength() const { return length; }
    size_t getLengthTrue() const { return cut_R <= cut_L ? 1 : cut_R - cut_L; }
    // number of bp that are trimmed off left side
    size_t getLTrim() const { return cut_L; }
    // number of bp that are trimmed off right side
    size_t getRTrim() const { return length - cut_R; }
};

typedef std::shared_ptr<Read> ReadPtr;
typedef std::vector<ReadPtr> Reads;

class ReadVisitor;

class ReadBase {
public:
    virtual ~ReadBase() {}
    ReadBase(ReadBase const&) = default;
    ReadBase() = default;

    Reads& get_reads_non_const() { return reads; }
    const Reads& get_reads() const { return reads; }

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

    virtual void accept(ReadVisitor &rv) = 0;

    bool rc;
    Reads reads;
};
typedef std::shared_ptr<ReadBase> ReadBasePtr;

class PairedEndRead: public ReadBase {
public:
    PairedEndRead() { reads.resize(2); }
    PairedEndRead(const Read& one_, const Read& two_) : PairedEndRead() {
        one = std::make_shared<Read>(one_);
        two = std::make_shared<Read>(two_);
        reads[0] = one;
        reads[1] = two;
    }

    PairedEndRead(const ReadPtr& one_, const ReadPtr& two_) : PairedEndRead() {
        one = one_;
        two = two_;
        reads[0] = one;
        reads[1] = two;
    }
    PairedEndRead(const std::vector<ReadPtr>& reads_) {
        reads = reads_;
        one = reads[0];
        two = reads[1];
    }

    virtual boost::optional<BitSet> get_key(size_t start, size_t length);
    Read& non_const_read_one() { return *one; }
    Read& non_const_read_two() { return *two; }
    const Read& get_read_one() const { return *one; }
    const Read& get_read_two() const { return *two; }
    virtual double avg_q_score();
    std::shared_ptr<ReadBase> convert(bool stranded);
    virtual void accept(ReadVisitor &rv);

private:
    ReadPtr one;
    ReadPtr two;
};
typedef std::shared_ptr<PairedEndRead> PairedEndReadPtr;

/* start, finish, discarded */
class SingleEndRead: public ReadBase {
public:
    SingleEndRead() { reads.resize(1); }
    SingleEndRead(const Read& one_) : SingleEndRead() {
        one = std::make_shared<Read>(one_);
        reads[0] = one;
    }
    SingleEndRead(const ReadPtr& one_) : SingleEndRead() {
        one = one_;
        reads[0] = one;
    }
    SingleEndRead(const std::vector<ReadPtr> reads_) {
        reads = reads_;
        one = reads[0];
    }

    virtual boost::optional<BitSet> get_key(size_t start, size_t length);
    Read& non_const_read_one() { return *one; }
    const Read& get_read() const { return *one; }
    virtual double avg_q_score();
    std::shared_ptr<ReadBase> convert(bool stranded);
    void set_read_rc() { one->set_read_rc();}
    virtual void accept(ReadVisitor &rv);

private:
    ReadPtr one;
};
typedef std::shared_ptr<SingleEndRead> SingleEndReadPtr;

class ReadVisitor {
public:
    virtual ~ReadVisitor() {};
    virtual void visit(SingleEndRead* ser) = 0;
    virtual void visit(PairedEndRead* per) = 0;
};

template <typename ST, typename PT>
class ReadVisitorFunc : public ReadVisitor {
public:
    virtual ~ReadVisitorFunc() {}

    ReadVisitorFunc(ST&& single_visit_, PT&& paired_visit_) :
        single_visit(std::move(single_visit_)),
        paired_visit(std::move(paired_visit_)) {}


    virtual void visit(SingleEndRead* ser) {
        single_visit(ser);
    }

    virtual void visit(PairedEndRead* per) {
        paired_visit(per);
    }

private:
    ST single_visit;
    PT paired_visit;
};

// We use this make_read_visitor_func because class template arguments cannot be inferred.
template <typename ST, typename PT>
ReadVisitorFunc<ST, PT> make_read_visitor_func(ST&& single_visit, PT&& paired_visit) {
    return ReadVisitorFunc<ST, PT>(std::forward<ST>(single_visit), std::forward<PT>(paired_visit));
}

#endif
