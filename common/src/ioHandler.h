#ifndef IOHANDLER_H
#define IOHANDLER_H

#include <istream>
#include <fstream>
#include <memory>
#include <utility>

#include "read.h"
#include <boost/iostreams/concepts.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/device/file_descriptor.hpp>
#include <boost/iostreams/stream.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/device/mapped_file.hpp>

#include <boost/filesystem/path.hpp>
#include <boost/filesystem.hpp>


#include <iostream>
#include <string>


namespace bf = boost::filesystem;
namespace bi = boost::iostreams;

int check_open_r(const std::string& filename) ;
int check_exists(const std::string& filename, bool force, bool gzip, bool std_out) ;

class HtsOfstream {
private:
    bool gzip;
    bool force;
    bool std_out;
    std::string filename;
    std::shared_ptr<std::ostream> out = nullptr;

    void create_out() {
        out.reset(new bi::stream<bi::file_descriptor_sink> {check_exists(filename, force, gzip, std_out), bi::close_handle});
    }

public:
    ~HtsOfstream() {
        if (out) {
            std::flush(*out);
        }
    }
    
    HtsOfstream(std::string filename_, bool force_, bool gzip_, bool stdout_) : filename(filename_), force(force_), gzip(gzip_),
                                                                                std_out(stdout_)  { }
    
    HtsOfstream(std::shared_ptr<std::ostream> out_) : out(out_) { }

    template<class T>
    HtsOfstream& operator<< (T s) {
        if (!out) {
            create_out();
        }
        *out << s;
        return *this;
    } 

    void flush() {
        if (out) {
            std::flush(*out);
        }
    }
};


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

class OutputWriter {
public:
    virtual ~OutputWriter() {  }
    virtual void write(const PairedEndRead& read) { throw std::runtime_error("No PE implementation of write (Probably a SE read)"); }
    virtual void write(const SingleEndRead& read) { throw std::runtime_error("No SE implementaiton of write (Probably a PE read)"); }
    virtual void write(const ReadBase &read) { throw std::runtime_error("No ReadBase class, only accessable with tab"); } //maybe typecase eventually 
};

class SingleEndReadOutFastq : public OutputWriter {
public:
    SingleEndReadOutFastq(std::shared_ptr<HtsOfstream> &out_) : output(out_) { }
    ~SingleEndReadOutFastq() { output->flush(); }
    void write(const SingleEndRead &read) { format_writer(read.get_read()); }
protected:
    std::shared_ptr<HtsOfstream> output = nullptr;
    
    void format_writer(const Read &read) { 
       *output << "@" << read.get_id() << '\n' << read.get_seq() << "\n+\n" << read.get_qual() << '\n'; 
    }
    
};

class PairedEndReadOutFastq : public OutputWriter {
public:
    PairedEndReadOutFastq(std::shared_ptr<HtsOfstream> &out1_, std::shared_ptr<HtsOfstream> &out2_) : out1(out1_), out2(out2_) { }
    ~PairedEndReadOutFastq() { out1->flush(); out2->flush(); }
    void write(const PairedEndRead &read) { format_writer(read.get_read_one(), read.get_read_two());  }
protected:
    std::shared_ptr<HtsOfstream> out1 = nullptr;
    std::shared_ptr<HtsOfstream> out2 = nullptr;
    
    void format_writer(const Read &read1, const Read &read2) { 
        *out1 << "@" << read1.get_id() << '\n' << read1.get_seq() << "\n+\n" << read1.get_qual() << '\n'; 
        *out2 << "@" << read2.get_id() << '\n' << read2.get_seq() << "\n+\n" << read2.get_qual() << '\n'; 
    }
};

class PairedEndReadOutInter : public OutputWriter {
public:
    PairedEndReadOutInter(std::shared_ptr<HtsOfstream> &out_) : out1(out_) { }
    ~PairedEndReadOutInter() { out1->flush(); }
    void write(const PairedEndRead &read) { format_writer(read.get_read_one(), read.get_read_two()); }
protected:
    std::shared_ptr<HtsOfstream> out1 = nullptr;
    void format_writer(const Read &read1, const Read &read2) { 
        *out1 << "@" << read1.get_id() << '\n' << read1.get_seq() << "\n+\n" << read1.get_qual() << '\n';
        *out1 << "@" << read2.get_id() << '\n' << read2.get_seq() << "\n+\n" << read2.get_qual() << '\n'; 
    }
};

class ReadBaseOutTab : public OutputWriter {
public:
    ReadBaseOutTab(std::shared_ptr<HtsOfstream> &out_) : output(out_) { }
    ~ReadBaseOutTab() { output->flush(); }
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
    std::shared_ptr<HtsOfstream> output = nullptr;
    
    void format_writer(const Read &read) { 
        *output << read.get_id() << '\t' << read.get_seq() << '\t' << read.get_qual() << '\n'; 
    }

    void format_writer(const Read &read1, const Read &read2) {
        *output << read1.get_id() << '\t' << read1.get_seq() << '\t' << read1.get_qual() << '\t' << read2.get_seq() << '\t' << read2.get_qual() << '\n';
    }
    
};


#endif
