#include "ioHandler.h"
#include <exception>
#include <cerrno>

void writer_helper(ReadBase *r, std::shared_ptr<OutputWriter> pe, std::shared_ptr<OutputWriter> se, bool stranded, bool no_orphans ) { //class writer
    PairedEndRead *per = dynamic_cast<PairedEndRead*>(r);
    if (per) {
        Read &one = per->non_const_read_one();
        Read &two = per->non_const_read_two();

        if (!one.getDiscard() && !two.getDiscard()) {
            pe->write(*per);
        } else if (!one.getDiscard() && !no_orphans) { //if stranded RC
            se->write_read(one, false);
        } else if (!two.getDiscard() && !no_orphans) { // Will never be RC
            se->write_read((per->get_read_two()), stranded);
        } else {

        }
    } else {
        SingleEndRead *ser = dynamic_cast<SingleEndRead*>(r);
        if (!ser) {
            throw std::runtime_error("Unknow read found");
        }
        if (! (ser->non_const_read_one()).getDiscard() ) {
            se->write(*ser);
        } else {
            
        }
            
    }
}

void writer_helper(ReadBase *r, std::shared_ptr<OutputWriter> pe, std::shared_ptr<OutputWriter> se, bool stranded, Counter &c, bool no_orphans ) {
    PairedEndRead *per = dynamic_cast<PairedEndRead*>(r);
    if (per) {
        if (!(per->non_const_read_one()).getDiscard() && !(per->non_const_read_two()).getDiscard()) {
            ++c["PE_Out"];
            per->setStats(c);
            pe->write(*per);
        } else if (!(per->non_const_read_one()).getDiscard() && !no_orphans) { //if stranded RC
            ++c["SE_Out"];
            ++c["R2_Discarded"];
            per->setStats(c);
            se->write_read((per->get_read_one()), false);
        } else if (!(per->non_const_read_two()).getDiscard() && !no_orphans) { // Will never be RC
            ++c["SE_Out"];
            ++c["R1_Discarded"];
            per->setStats(c);
            se->write_read((per->get_read_two()), stranded);
        } else {
            ++c["R1_Discarded"];
            ++c["R2_Discarded"];
            per->setStats(c);
        }
    } else {
        SingleEndRead *ser = dynamic_cast<SingleEndRead*>(r);
        if (!ser) {
            throw std::runtime_error("Unknow read found");
        }
        if (! (ser->non_const_read_one()).getDiscard() ) {
            ++c["SE_Out"];
            ser->setStats(c);
            se->write(*ser);
        } else {
            ++c["SE_Discarded"];
            ser->setStats(c);
        }
            
    }
}

void skip_lr(std::istream *input) {
    while(input and input->good() and (input->peek() == '\n' || input->peek() == '\r')) {
        input->get();
    }
}

void throw_error(const std::string& filename) {
    throw std::runtime_error(filename + ": " +  std::strerror( errno ));
}


int check_open_r(const std::string& filename) {
    FILE* f = NULL;

    bf::path p(filename);
    if (!bf::exists(p)) {
        throw std::runtime_error("File " + filename + " was not found.");
    }

    if (p.extension() == ".gz") {
        f = popen(("gunzip -c " + filename).c_str(), "r");
    } else {
        f = fopen(filename.c_str(), "r");
    }

    if (!f) {
        throw_error(filename);
    }
    return fileno(f);
}


int check_exists(const std::string& filename, bool force, bool gzip, bool std_out) {
    FILE* f = NULL;
    
    if (std_out) {
        return fileno(stdout);
    }
    bf::path p(filename);

    if (force || !bf::exists(p)) {
        if (gzip) {
            f = popen(("gzip > " + filename + ".gz").c_str(), "w");
        } else {
            f = fopen(filename.c_str(), "w");
        }
        if (!f) {
            throw_error(filename);
        }
        return fileno(f);
    } else {
        throw std::runtime_error("File " + filename + " all ready exists. Please use -F or delete it\n");
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
    while(input->good() and (input->peek() == '\n' || input->peek() == '\r')) {
        input->get();
    }
    return Read(seq, qual, id.substr(1));
}

Read InputFasta::load_read(std::istream *input) {
    while(std::getline(*input, id) && id.size() < 1) {
    }
    if (id.size() < 1) {
        throw std::runtime_error("invalid id - line empty");
    }
    if (id[0] != '>') {
        throw std::runtime_error("id line did not begin with >");
    }
    seq = "";
    while (std::getline(*input, tmpSeq)) {
        seq += tmpSeq;
        if (input->peek() == '>') {
            break;
        }
    }
    if (seq.size() < 1) {
        throw std::runtime_error("no sequence");  
    }
    while(input->good() and (input->peek() == '\n' || input->peek() == '\r')) {
        input->get();
    }
    return Read(seq, "", id);

}

//Overrides load_read for tab delimited reads
std::vector<Read> TabReadImpl::load_read(std::istream *input) {

    std::vector <Read> reads(1);
    while(std::getline(*input, tabLine) && tabLine.size() < 1) {
    }

    std::vector <std::string> parsedRead;
    boost::split(parsedRead, tabLine, boost::is_any_of("\t"));

    
    if (parsedRead.size() != 3 && parsedRead.size() != 5) {
        throw std::runtime_error("There are not either 3 or 5 elements within a tab delimited file line");
    }

    if (parsedRead[1].size() != parsedRead[2].size()) {
        throw std::runtime_error("sequence and qualities are not the same length 1");
    }
    
    reads[0] = Read(parsedRead[1], parsedRead[2], parsedRead[0]);
   
    if (parsedRead.size() != 3) {
        
        if (parsedRead[3].size() != parsedRead[4].size()) {
            throw std::runtime_error("sequence and qualities are not the same length 2");
        }
   
        reads.push_back(Read(parsedRead[3], parsedRead[4], parsedRead[0]));
    }

    // ignore extra lines at end of file
    while(input->good() and (input->peek() == '\n' || input->peek() == '\r')) {
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
bool InputReader<SingleEndRead, FastaReadImpl>::has_next() {
    // ignore extra lines at end of file
    skip_lr(input);
    return (input and input->good());
};

template<>
InputReader<SingleEndRead, FastaReadImpl>::value_type InputReader<SingleEndRead, FastaReadImpl>::next() {
    return InputReader<SingleEndRead, FastaReadImpl>::value_type(new SingleEndRead(load_read(input)));
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

template<>
InputReader<PairedEndRead, InterReadImpl>::value_type InputReader<PairedEndRead, InterReadImpl>::next() {
    Read r1 = load_read(in1);
    Read r2;
    try {
        r2 = load_read(in1);
    } catch (const std::exception&) {
        throw std::runtime_error("odd number of sequences in interleaved file");
    }

    return InputReader<PairedEndRead, InterReadImpl>::value_type(new PairedEndRead(r1, r2));
}

template<>
bool InputReader<PairedEndRead, InterReadImpl>::has_next() {
    skip_lr(in1);
    return(in1 && in1->good());
}

template <>
InputReader<ReadBase, TabReadImpl>::value_type InputReader<ReadBase, TabReadImpl>::next() {
    std::vector<Read> rs = load_read(in1);
    if (rs.size() == 1) {
        return InputReader<SingleEndRead, TabReadImpl>::value_type(new SingleEndRead(rs[0]));
    }

    return InputReader<PairedEndRead, TabReadImpl>::value_type(new PairedEndRead(rs[0], rs[1]));
}

template <>
bool InputReader<ReadBase, TabReadImpl>::has_next() {
    // ignore extra lines at end of file
    skip_lr(in1);
    return (in1 and in1->good());
}
