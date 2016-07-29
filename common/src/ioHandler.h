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

class inputFastqSingle {
public:
    inputFastqSingle(std::istream& in) : input(&in) {}

    class iterator {
    public:
        typedef SingleEndRead value_type;
        typedef SingleEndRead& reference;
        typedef SingleEndRead* pointer;
        typedef std::input_iterator_tag iterator_category; //or another tag

        iterator(std::istream *in);
        bool operator==(const iterator&) const;
        bool operator!=(const iterator&) const;
   /*
        iterator (const iterator&);
        ~iterator();
        iterator& operator=(const iterator&);
        */

        iterator& operator++();
        reference operator*() const;
        pointer operator->() const;

    private:
        void get_read();

        std::istream* input;
        std::unique_ptr<SingleEndRead> read = 0;

    };

    iterator begin();
    iterator end();

private:
    std::istream* input = 0;
};


class inputFastqPaired {
public:
    inputFastqPaired(std::istream& in1, std::istream& in2) : in1(&in1), in2(&in2) {}

    class iterator {
    public:
        typedef PairedEndRead value_type;
        typedef PairedEndRead& reference;
        typedef PairedEndRead* pointer;
        typedef std::input_iterator_tag iterator_category; //or another tag

        iterator(std::istream *in1, std::istream *in2);
        bool operator==(const iterator&) const;
        bool operator!=(const iterator&) const;
   /*
        iterator (const iterator&);
        ~iterator();
        iterator& operator=(const iterator&);
        */

        iterator& operator++();
        reference operator*() const;
        pointer operator->() const;

    private:
        void get_read();

        std::istream* in1, *in2 = 0;
        std::unique_ptr<PairedEndRead> read = 0;

    };

    iterator begin();
    iterator end();

private:
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
