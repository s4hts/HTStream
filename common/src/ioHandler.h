#ifndef IOHANDLER_H
#define IOHANDLER_H

#include <istream>
#include <memory>
#include "read.h"
#include <boost/iostreams/concepts.hpp>

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
    PairedEndReadFastqImpl(std::istream& in1_, std::istream& in2_) : in1(&in1_), in2(&in2_) {
    }
  
protected:
    std::istream* in1, * in2 = 0;
};

// ### output ###

template <class T, class Impl>
class OutputWriter : Impl {
public:
    using Impl::Impl;
    void write(const T& data);
};

class OutputFastq {
protected:
    void write_read(const Read& read, std::ostream &output);
};

class SingleEndReadOutFastq : public OutputFastq {
public:
    SingleEndReadOutFastq(std::ostream& out) : output(out) {}
    ~SingleEndReadOutFastq() { output.flush(); }
protected:
    std::ostream& output;
};

class PairedEndReadOutFastq : public OutputFastq {
public:
    PairedEndReadOutFastq(std::ostream& out1_, std::ostream& out2_) : out1(out1_), out2(out2_) {}
    ~PairedEndReadOutFastq() { out1.flush(); out2.flush(); }
protected:
    std::ostream &out1, &out2;
};

#endif
