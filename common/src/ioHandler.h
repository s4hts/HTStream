#ifndef IOHANDLER_H
#define IOHANDLER_H

#include <istream>
#include <memory>
#include "read.h"


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

class SingleEndReadImpl : public InputFastq{
public:
    SingleEndReadImpl(std::istream& in) : input(&in) {}
    
protected:
    std::istream* input = 0;
};


class PairedEndReadImpl : public InputFastq {
public:
    PairedEndReadImpl(std::istream& in1_, std::istream& in2_) : in1(&in1_), in2(&in2_) {}
  
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

class SingleEndReadOut : public OutputFastq {
public:
    SingleEndReadOut(std::ostream& out) : output(out) {}
    ~SingleEndReadOut() { output.flush(); }
protected:
    std::ostream& output;
};

class PairedEndReadOut : public OutputFastq {
public:
    PairedEndReadOut(std::ostream& out1_, std::ostream& out2_) : out1(out1_), out2(out2_) {}
    ~PairedEndReadOut() { out1.flush(); out2.flush(); }
protected:
    std::ostream &out1, &out2;
};

#endif
