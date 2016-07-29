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
        read = get_read();
    }
};

std::unique_ptr<SingleEndRead> inputFastqSingle::iterator::get_read() {
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
    return std::unique_ptr<SingleEndRead>(new SingleEndRead(Read(seq, qual), id.substr(1)));
}
        
inputFastqSingle::iterator& inputFastqSingle::iterator::operator++(
    ) {
    if(input and input->good()) {
        read = get_read();
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
