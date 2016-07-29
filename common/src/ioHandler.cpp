#include "ioHandler.h"
#include <exception>

/*
inputBase::iterator inputBase::begin() {
    return inputBase::iterator(&input);
}

inputBase::iterator inputBase::end() {
    return inputBase::iterator(0);
}
*/

inputFastqSingle::iterator inputFastqSingle::begin() {
    return inputFastqSingle::iterator(input);
}

inputFastqSingle::iterator inputFastqSingle::end() {
    return inputFastqSingle::iterator(0);
}

inputFastqSingle::iterator::iterator(std::istream *in) : input(in) {
    if (input and input->good()) {
        get_read();
    } else {
        input = nullptr;
    }
};

Read load_read(std::istream *input) {
    std::string id, seq, id2, qual;
    while(std::getline(*input, id) && id.size() < 1) {
    }
    if (id.size() < 1) {
        throw std::runtime_error("invalid id line empty");
    }
    if (id[0] != '@') {
        throw std::runtime_error("id line did not begin with @");
    }
    std::getline(*input, seq);
    if (seq.size() < 1) {
        throw std::runtime_error("invalid seq line empty");
    }
    std::getline(*input, id2);
    if (id2.size() < 1) {
        throw std::runtime_error("invalid id2 line empty");
    }
    if (id2[0] != '+') {
        throw std::runtime_error("invalid id2 line did not begin with +");
    }
    std::getline(*input, qual);
    if (qual.size() != seq.size()) {
        throw std::runtime_error("qual string not the same length as sequence");
    }

    // ignore extra lines at end of file
    while(input->good() and input->peek() == '\n') {
        input->get();
    }
    return Read(seq, qual, id.substr(1));
}

void inputFastqSingle::iterator::get_read() {
    this->read.reset(new SingleEndRead(load_read(this->input)));
}
        
inputFastqSingle::iterator& inputFastqSingle::iterator::operator++(
    ) {
    if(input and input->good()) {
        get_read();
        return *this;
    } else {
        read.reset();
        input = nullptr;
        return *this;
    }
}

SingleEndRead& inputFastqSingle::iterator::operator*() const {
    return *read;
}

SingleEndRead* inputFastqSingle::iterator::operator->() const {
    return read.get();
}

bool inputFastqSingle::iterator::operator==(const iterator& other) const {
    return (other.read == this->read && other.input == this->input);
}

bool inputFastqSingle::iterator::operator!=(const iterator& other) const {
    return !(other == *this);
}

// inputFastqPaired 

inputFastqPaired::iterator inputFastqPaired::begin() {
    return inputFastqPaired::iterator(in1, in2);
}

inputFastqPaired::iterator inputFastqPaired::end() {
    return inputFastqPaired::iterator(0, 0);
}


void inputFastqPaired::iterator::get_read() {
    Read r1 = load_read(in1);
    Read r2 = load_read(in2);
    this->read.reset(new PairedEndRead(r1, r2));
}


inputFastqPaired::iterator::iterator(std::istream *in1, std::istream *in2) : in1(in1), in2(in2) {
    if (in1 and in1->good() and in2 and in2->good()) {
        get_read();
    } else {
        in1 = nullptr;
        in2 = nullptr;
    }
};

        
inputFastqPaired::iterator& inputFastqPaired::iterator::operator++(
    ) {
    if(in1 and in1->good() and in2 and in2->good()) {
        get_read();
        return *this;
    } else {
        read.reset();
        in1 = nullptr;
        in2 = nullptr;
        return *this;
    }
}

PairedEndRead& inputFastqPaired::iterator::operator*() const {
    return *read;
}

PairedEndRead* inputFastqPaired::iterator::operator->() const {
    return read.get();
}

bool inputFastqPaired::iterator::operator==(const iterator& other) const {
    return (other.read == this->read && other.in1 == this->in1 && other.in2 == this->in2);
}

bool inputFastqPaired::iterator::operator!=(const iterator& other) const {
    return !(other == *this);
}
