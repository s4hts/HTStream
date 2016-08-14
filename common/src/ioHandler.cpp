#include "ioHandler.h"
#include <exception>
    
void skip_lr(std::istream *input) {
    while(input and input->good() and (input->peek() == '\n' || input->peek() == '\r')) {
        input->get();
    }
}

Read InputFastq::load_read(std::istream *input) {
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

//Overrides load_read for tab delimited reads
std::vector<Read> TabReadImpl::load_read(std::istream *input) {
    //to read the line
    std::string tabLine;

    std::vector <Read> reads(2);
    while(std::getline(*input, tabLine) && tabLine.size() < 1) {
    }

    std::vector <std::string> parsedRead;
    boost::split(parsedRead, tabLine, boost::is_any_of("\t"));

    if (parsedRead.size() != 3 || parsedRead.size() != 5) {
        throw std::runtime_error("There are not either 3 or 5 elements within a tab delimited file line");
    }
    if (parsedRead[1].size() != parsedRead[2].size()) {
        throw std::runtime_error("sequence and qualities are not the same length");
    }
    
    reads[0] = Read(parsedRead[0], parsedRead[1], parsedRead[2]);
    
    if (parsedRead.size() == 5) {
        reads[1] = Read(parsedRead[0], parsedRead[3], parsedRead[4]);
    }
    
    if (parsedRead[3].size() != parsedRead[4].size()) {
        throw std::runtime_error("sequence and qualities are not the same length");
    }
    // ignore extra lines at end of file
    while(input->good() and input->peek() == '\n') {
        input->get();
    }

    return reads;
}

template <>
bool InputReader<SingleEndRead, SingleEndReadFastqImpl>::has_next() {
    // ignore extra lines at end of file
    skip_lr(input);
    return (input and input->good());
};


template <>
InputReader<SingleEndRead, SingleEndReadFastqImpl>::value_type InputReader<SingleEndRead, SingleEndReadFastqImpl>::next() {
    return InputReader<SingleEndRead, SingleEndReadFastqImpl>::value_type(new SingleEndRead(load_read(input)));
}

template <>
InputReader<PairedEndRead, PairedEndReadFastqImpl>::value_type InputReader<PairedEndRead, PairedEndReadFastqImpl>::next() {
    Read r1 = load_read(in1);
    Read r2 = load_read(in2);
    return InputReader<PairedEndRead, PairedEndReadFastqImpl>::value_type(new PairedEndRead(r1, r2));

}

template <>
bool InputReader<PairedEndRead, PairedEndReadFastqImpl>::has_next() {
    // ignore extra lines at end of file
    skip_lr(in1);
    skip_lr(in2);
    return (in1 and in1->good() and in2 and in2->good());
};

template <>
InputReader<TabRead, TabReadImpl>::value_type InputReader<TabRead, TabReadImpl>::next() {
    /*Read 2 could be null*/
    std::vector<Read> rs = load_read(in1);

    return InputReader<TabRead, TabReadImpl>::value_type(new TabRead(rs[0], rs[1]));
}

template <>
bool InputReader<TabRead, TabReadImpl>::has_next() {
    // ignore extra lines at end of file
    skip_lr(in1);
    return (in1 and in1->good());
}
// ### output ###
void OutputFastq::write_read(const Read& read, std::ostream &output) {
    output << '@' << read.get_id() << "\n";
    output << read.get_seq() << "\n";
    output << '+' << "\n";
    output << read.get_qual() << "\n";
}

template <>
void OutputWriter<SingleEndRead, SingleEndReadOutFastq>::write(const SingleEndRead& data) {
    write_read(data.get_read(), output);
}

template <>
void OutputWriter<PairedEndRead, PairedEndReadOutFastq>::write(const PairedEndRead& data) {
    write_read(data.get_read_one(), out1);
    write_read(data.get_read_two(), out2);
}
