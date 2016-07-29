#ifndef IOHANDLER_H
#define IOHANDLER_H

#include <istream>
#include <memory>
#include "read.h"
/*
class inputBase {
public:
    inputBase(std::istream& in) : input(in) {};

    
private:
    std::istream& input;
};
*/

class inputTab {
public:
    inputTab(std::istream& in);
    
};


template <class T, class Impl>
class InputReader : Impl {
public:
    typedef std::unique_ptr<T> value_type;
    using Impl::Impl;
    
    bool has_next();
    value_type next();
        
};


class SingleEndReadImpl {
public:
    SingleEndReadImpl(std::istream& in) : input(&in) {}
    
protected:
    std::istream* input = 0;
};


class PairedEndReadImpl {
public:
    PairedEndReadImpl(std::istream& in1, std::istream& in2) : in1(&in1), in2(&in2) {}
  
protected:
    std::istream* in1, * in2 = 0;
};

/*
class outputTab {
};

class outputFastq
{
}
*/

#endif
