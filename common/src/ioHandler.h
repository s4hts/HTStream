#ifndef IOHANDLER_H
#define IOHANDLER_H

#include <istream>
#include <memory>
#include "read.h"
#include <boost/iostreams/concepts.hpp>
#include <boost/algorithm/string.hpp>
#include <iostream>
// ### input ###
template <class T, class Impl>
class InputReader : Impl {
public:
    typedef std::unique_ptr<T> value_type;
    using Impl::Impl;
    
    bool has_next();
    value_type next();
};

class InputFastq {
protected:
    Read load_read(std::istream *input);
    
    std::string id, seq, id2, qual;
};

class SingleEndReadFastqImpl : public InputFastq{
public:
    SingleEndReadFastqImpl(std::istream& in) : input(&in) {}
    
protected:
    std::istream* input = 0;
};

class PairedEndReadFastqImpl : public InputFastq {
public:
    PairedEndReadFastqImpl(std::istream& in1_, std::istream& in2_) : in1(&in1_), in2(&in2_) {}
  
protected:
    std::istream* in1, * in2 = 0;
};

class TabReadImpl : public InputFastq {
public:
    TabReadImpl(std::istream& in1_) : in1(&in1_) {}
    std::vector<Read> load_read(std::istream *input);
protected:
    std::istream* in1;
    //to read the line
    std::string tabLine;
};

class InterReadImpl : public InputFastq {
public:
    InterReadImpl(std::istream& in1_) : in1(&in1_) {}
protected:
    std::istream *in1;
};

// ### output ###

/*template <class T, class Impl>
class OutputWriter : Impl {
public:
    using Impl::Impl;
    void write(const T& data);
};*/

class OutputWriter {
public:
    virtual ~OutputWriter() { }
    virtual void write(const PairedEndRead& read) { throw std::runtime_error("No PE implementation of write (Probably a SE read)"); }
    virtual void write(const SingleEndRead& read) { throw std::runtime_error("No SE implementaiton of write (Probably a PE read)"); }
    virtual void write(const ReadBase &read) { throw std::runtime_error("No ReadBase class, only accessable with tab"); } //maybe typecase eventually 
};

class SingleEndReadOutFastq : public OutputWriter {
public:
    SingleEndReadOutFastq(std::ostream& out) : output(out) {}
    ~SingleEndReadOutFastq() { output.flush(); }
    void write(const SingleEndRead &read) { format_writer(read.get_read()); }
protected:
    void format_writer(const Read &read) { output << "@" << read.get_id() << '\n' << read.get_seq() << "\n+\n" << read.get_qual() << '\n'; }
    std::ostream& output;
};

class PairedEndReadOutFastq : public OutputWriter {
public:
    PairedEndReadOutFastq(std::ostream &out1_, std::ostream &out2_) : out1(out1_), out2(out2_) {}
    ~PairedEndReadOutFastq() { out1.flush(); out2.flush(); }
    void write(const PairedEndRead &read) { format_writer(read.get_read_one(), read.get_read_two()); }
protected:
    std::ostream &out1, &out2;
    void format_writer(const Read &read1, const Read &read2) { out1 << "@" << read1.get_id() << '\n' << read1.get_seq() << "\n+\n" << read1.get_qual() << '\n' ; out2 << "@" << read2.get_id() << '\n' << read2.get_seq() << "\n+\n" << read2.get_qual() << '\n'; }
};

class PairedEndReadOutInter : public OutputWriter {
public:
    PairedEndReadOutInter(std::ostream& out1_) : out1(out1_) {}
    ~PairedEndReadOutInter() { out1.flush(); }
    void write(const PairedEndRead &read) { format_writer(read.get_read_one(), read.get_read_two()); }
protected:
    std::ostream &out1;
    void format_writer(const Read &read1, const Read &read2) { out1 << "@" << read1.get_id() << '\n' << read1.get_seq() << "\n+\n" << read1.get_qual() << '\n' ; out1 << "@" << read2.get_id() << '\n' << read2.get_seq() << "\n+\n" << read2.get_qual() << '\n'; }
};

class ReadBaseOutTab : public OutputWriter {
public:
    ReadBaseOutTab(std::ostream& _output) : output(_output) {}
    ~ReadBaseOutTab() { output.flush(); }
    void write(const PairedEndRead &read) { format_writer(read.get_read_one(), read.get_read_two()); }
    void write(const SingleEndRead &read) { format_writer(read.get_read()); }
    void write(const ReadBase &read) {  
        const PairedEndRead *per = dynamic_cast<const PairedEndRead*>(&read);
        if (per) {
            format_writer(per->get_read_one(), per->get_read_two());
        } else {
            const SingleEndRead *ser = dynamic_cast<const SingleEndRead*>(&read);
            if (ser == NULL) {
                throw std::runtime_error("output tab could not cast to SE or PE read");
            }
            format_writer(ser->get_read());
        }    
    }
protected:
    void format_writer(const Read &read) { output << read.get_id() << '\t' << read.get_seq() << '\t' << read.get_qual() << '\n' ; }
    void format_writer(const Read &read1, const Read &read2) { output << read1.get_id() << '\t' << read1.get_seq() << '\t' << read1.get_qual() << '\t' << read2.get_seq() << '\t' << read2.get_qual() << '\n'; }
    std::ostream &output;
};


#endif
