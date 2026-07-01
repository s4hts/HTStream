#ifndef IOHANDLER_H
#define IOHANDLER_H

#define BOOST_DYNAMIC_BITSET_DONT_USE_FRIENDS

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
#include <cstdio>

#include "counters.h"
#include "hts_exception.h"

namespace bf = boost::filesystem;
namespace bi = boost::iostreams;

int check_open_r(const std::string& filename);
void set_gzip_compression_level(size_t level);
std::string string2fasta(std::string seqstring, std::string prefix, const char delim=',');
Read fasta_to_read(std::string fasta_file);

class HtsOfstream {
private:
    std::string filename;
    bool force;
    bool gzip;
    bool std_out;
    std::shared_ptr<std::ostream> out = nullptr;
    FILE* gzfile = 0;

    int check_exists(const std::string& filename, bool force, bool gzip, bool std_out) ;

    void create_out() {
        out.reset(new bi::stream<bi::file_descriptor_sink> (
                check_exists(filename, force, gzip, std_out), bi::close_handle));
    }

public:
    ~HtsOfstream() {
        if (out) {
            out->flush();
            if (gzfile) {
                pclose(gzfile);
            }
            out.reset();
        }
    }

    HtsOfstream(std::string filename_, bool force_, bool gzip_, bool stdout_) :
        filename(filename_), force(force_), gzip(gzip_), std_out(stdout_)  { }

    HtsOfstream(std::shared_ptr<std::ostream> out_) : out(out_) { }

    template<class T>
    HtsOfstream& operator<< (T s) {
        if (!out) {
            create_out();
        }
        *out << s;
        return *this;
    }

    void write(const std::string& s) {
        if (!out) {
            create_out();
        }
        out->write(s.data(), static_cast<std::streamsize>(s.size()));
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
    ReadPtr load_read(std::istream *input);

    std::string id, seq, id2, qual;
};

class InputFasta {
protected:
    ReadPtr load_read(std::istream *input);
    std::string id, seq;
    std::string tmpSeq;
};

class FastaReadImpl : public InputFasta {
public:
    FastaReadImpl(std::istream& input_) : input(&input_) {}
protected:
    std::istream* input = 0;
};

class PairedEndReadFastqImpl : public InputFastq {
public:
    PairedEndReadFastqImpl(std::istream& in1_, std::istream& in2_) : in1(&in1_), in2(&in2_) {}
    PairedEndReadFastqImpl(std::vector<std::string> in1_, std::vector<std::string> in2_) : fin1(in1_), fin2(in2_) {}
protected:
    std::istream* in1, * in2 = 0;
    std::vector<std::string> fin1, fin2;
    bi::stream<bi::file_descriptor_source> fs1;
    bi::stream<bi::file_descriptor_source> fs2;
};

class InterReadImpl : public InputFastq {
public:
    InterReadImpl(std::istream& in1_) : in1(&in1_) {}
    InterReadImpl(std::vector<std::string> in_) : fin(in_) {}
protected:
    std::istream *in1;
    std::vector<std::string> fin;
    bi::stream<bi::file_descriptor_source> inter;
};

class SingleEndReadFastqImpl : public InputFastq {
public:
    SingleEndReadFastqImpl(std::istream& in_) : input(&in_) {}
    SingleEndReadFastqImpl(std::vector<std::string> in_) : finput(in_) {}
protected:
    std::istream* input = 0;
    std::vector<std::string> finput;
    bi::stream<bi::file_descriptor_source> fs;
};

class TabReadImpl : public InputFastq {
public:
    TabReadImpl(std::istream& in1_) : in1(&in1_) {}
    TabReadImpl(std::vector<std::string> in_) : fin(in_) {}
    std::vector<ReadPtr> load_read(std::istream *input);
protected:
    std::istream* in1;
    std::vector<std::string> fin;
    bi::stream<bi::file_descriptor_source> tabin;
    //to read the line
    std::string tabLine;
};


class inputReaders {
private:
  std::vector<std::string> r1_input;
  std::vector<std::string> r2_input;
  std::vector<std::string> se_input;
  std::vector<std::string> interleaved_input;
  std::vector<std::string> tab_input;

public:
    inputReaders(std::vector<std::string> r1_input_, std::vector<std::string> r2_input_,
                 std::vector<std::string> se_input_,
                 std::vector<std::string> interleaved_input_,
                 std::vector<std::string> tab_input_) :
        r1_input(r1_input_), r2_input(r2_input_), se_input(se_input_), interleaved_input(interleaved_input_), tab_input(tab_input_) { }

};

class OutputWriter {
public:
    virtual ~OutputWriter() {  }
    virtual void write(const PairedEndRead& ) { throw HtsRuntimeException("No PE implementation of write (Probably a SE read)"); }
    virtual void write(const SingleEndRead& ) { throw HtsRuntimeException("No SE implementaiton of write (Probably a PE read)"); }
    virtual void write_read(const Read &, bool ) { throw HtsRuntimeException("No write_read class, only accessable with SE"); } //only SE
};

class SingleEndReadOutFastq : public OutputWriter {
public:
    SingleEndReadOutFastq(std::shared_ptr<HtsOfstream> &out_) : output(out_) { }
    virtual ~SingleEndReadOutFastq() {};
    virtual void write(const SingleEndRead &read) { format_writer(read.get_read()); }
    virtual void write_read(const Read &read, bool rc) { if (rc) { format_writer_rc(read); } else { format_writer(read); } }
protected:
    std::shared_ptr<HtsOfstream> output = nullptr;

    void format_writer_rc(const Read &read) {
        std::string record;
        record.reserve(read.getLengthTrue() * 2 + read.get_id_first().size() + 8);
        record += '@';
        record += read.get_id_fastq('1');
        record += '\n';
        record += read.get_seq_rc();
        record += "\n+\n";
        record += read.get_qual_rc();
        record += '\n';
        output->write(record);
    }
    void format_writer(const Read &read) {
        std::string record;
        record.reserve(read.getLengthTrue() * 2 + read.get_id_first().size() + 8);
        record += '@';
        record += read.get_id_fastq('1');
        record += '\n';
        read.append_sub_seq(record);
        record += "\n+\n";
        read.append_sub_qual(record);
        record += '\n';
        output->write(record);
    }

};

class PairedEndReadOutFastq : public OutputWriter {
public:
    PairedEndReadOutFastq(std::shared_ptr<HtsOfstream> &out1_, std::shared_ptr<HtsOfstream> &out2_) : out1(out1_), out2(out2_) { }
    virtual ~PairedEndReadOutFastq() {};
    virtual void write(const PairedEndRead &read) { format_writer(read.get_read_one(), read.get_read_two()); }
protected:
    std::shared_ptr<HtsOfstream> out1 = nullptr;
    std::shared_ptr<HtsOfstream> out2 = nullptr;

    void format_writer(const Read &read1, const Read &read2) {
        std::string record1;
        record1.reserve(read1.getLengthTrue() * 2 + read1.get_id_first().size() + 8);
        record1 += '@';
        record1 += read1.get_id_fastq('1');
        record1 += '\n';
        read1.append_sub_seq(record1);
        record1 += "\n+\n";
        read1.append_sub_qual(record1);
        record1 += '\n';
        out1->write(record1);

        std::string record2;
        record2.reserve(read2.getLengthTrue() * 2 + read2.get_id_first().size() + 8);
        record2 += '@';
        record2 += read2.get_id_fastq('2');
        record2 += '\n';
        read2.append_sub_seq(record2);
        record2 += "\n+\n";
        read2.append_sub_qual(record2);
        record2 += '\n';
        out2->write(record2);
    }
};

class PairedEndReadOutInter : public OutputWriter {
public:
    PairedEndReadOutInter(std::shared_ptr<HtsOfstream> &out_) : out1(out_) { }
    virtual ~PairedEndReadOutInter() {};
    virtual void write(const PairedEndRead &read) { format_writer(read.get_read_one(), read.get_read_two()); }
protected:
    std::shared_ptr<HtsOfstream> out1 = nullptr;

    void format_writer(const Read &read1, const Read &read2) {
        std::string record;
        record.reserve((read1.getLengthTrue() + read2.getLengthTrue()) * 2 + read1.get_id_first().size() + read2.get_id_first().size() + 16);
        record += '@';
        record += read1.get_id_fastq('1');
        record += '\n';
        read1.append_sub_seq(record);
        record += "\n+\n";
        read1.append_sub_qual(record);
        record += "\n@";
        record += read2.get_id_fastq('2');
        record += '\n';
        read2.append_sub_seq(record);
        record += "\n+\n";
        read2.append_sub_qual(record);
        record += '\n';
        out1->write(record);
    }
};

/*Unmapped reads*/
class ReadBaseOutUnmapped : public OutputWriter {
public:
    ReadBaseOutUnmapped(std::shared_ptr<HtsOfstream> &out_) : output(out_) { }
    virtual ~ReadBaseOutUnmapped() {};
    virtual void write(const PairedEndRead &read) { format_writer(read.get_read_one(), read.get_read_two()); }
    virtual void write(const SingleEndRead &read) { format_writer(read.get_read()); }
    virtual void write_read(const Read &read, bool rc) { if (rc) { format_writer_rc(read); } else { format_writer(read); } }

protected:
    std::shared_ptr<HtsOfstream> output = nullptr;

    /*sam format specs spaces are for readability
     * id \t bitwise flag \t rname \t pos \t mapQ \t CIGAR \t RNEXT \t PNEXT \t TLEN \t SEQ \t QUAL\n
     *
     * id = id
     * bitwas flag
     * SE - 4
     * PE R1 - 77
     * PE R2 - 141
     * RNAME - *
     * POS - 0
     * MAPQ - 0
     * CIGAR - *
     * RNEXT - *
     * PNEXT - 0
     * TLEN - 0
     * SEQ - seq
     * QUAL - qual */
    const size_t se_bitwise = 4;
    const size_t pe1_bitwise = 77;
    const size_t pe2_bitwise = 141;

    void samout(const Read &read, size_t bitwiseflag) {
        std::string record;
        record.reserve(read.getLengthTrue() * 2 + read.get_id_first().size() + 64);
        record += read.get_id_first();
        record += '\t';
        record += std::to_string(bitwiseflag);
        record += "\t*\t0\t0\t*\t*\t0\t0\t";
        read.append_sub_seq(record);
        record += '\t';
        read.append_sub_qual(record);
        for (auto const& s : read.get_comment()) {
            record += '\t';
            record += s;
        }
        record += '\n';
        output->write(record);
    }

    void samout_rc(const Read &read, size_t bitwiseflag) {
        std::string record;
        record.reserve(read.getLengthTrue() * 2 + read.get_id_first().size() + 64);
        record += read.get_id_first();
        record += '\t';
        record += std::to_string(bitwiseflag);
        record += "\t*\t0\t0\t*\t*\t0\t0\t";
        record += read.get_seq_rc();
        record += '\t';
        record += read.get_qual_rc();
        for (auto const& s : read.get_comment()) {
            record += '\t';
            record += s;
        }
        record += '\n';
        output->write(record);
    }

    /*Unmapped specs for SE reads*/
    void format_writer(const Read &read) {
        samout(read, se_bitwise);
    }

    void format_writer(const Read &read1, const Read &read2) {
        samout(read1, pe1_bitwise);
        samout(read2, pe2_bitwise);
    }

    void format_writer_rc(const Read &read) {
       samout_rc(read, se_bitwise);
    }
};


class ReadBaseOutTab : public OutputWriter {
public:
    ReadBaseOutTab(std::shared_ptr<HtsOfstream> &out_) : output(out_) { }
    virtual ~ReadBaseOutTab() {};
    virtual void write(const PairedEndRead &read) { format_writer(read.get_read_one(), read.get_read_two()); }
    virtual void write(const SingleEndRead &read) { format_writer(read.get_read()); }
    virtual void write_read(const Read &read, bool rc) { if (rc) { format_writer_rc(read); } else { format_writer(read); } }

protected:
    std::shared_ptr<HtsOfstream> output = nullptr;

    void format_writer(const Read &read) {
        std::string record;
        record.reserve(read.getLengthTrue() * 2 + read.get_id_first().size() + 8);
        record += read.get_id_tab('1');
        record += '\t';
        read.append_sub_seq(record);
        record += '\t';
        read.append_sub_qual(record);
        if (read.get_comment().size() > 0){
            record += '\t';
            record += strjoin(read.get_comment(), "|");
        }
        record += '\n';
        output->write(record);
    }

    void format_writer(const Read &read1, const Read &read2) {
        std::string record;
        record.reserve((read1.getLengthTrue() + read2.getLengthTrue()) * 2 + read1.get_id_first().size() + read2.get_id_first().size() + 16);
        record += read1.get_id_tab('1');
        record += '\t';
        read1.append_sub_seq(record);
        record += '\t';
        read1.append_sub_qual(record);
        record += '\t';
        record += read2.get_id_tab('2');
        record += '\t';
        read2.append_sub_seq(record);
        record += '\t';
        read2.append_sub_qual(record);

        if (read1.get_comment().size() > 0 || read2.get_comment().size() > 0){
            const std::vector <std::string>& comment1 = read1.get_comment();
            const std::vector <std::string>& comment2 = read2.get_comment();
            std::string strComment;
            if (comment1.size() > 0){
                strComment += strjoin(comment1, "|");
            }
            strComment += '\t';
            if (comment2.size() > 0){
                strComment += strjoin(comment2, "|");
            }
            record += '\t';
            record += strComment;
        }
        record += '\n';
        output->write(record);
    }

    void format_writer_rc(const Read &read) {
        std::string record;
        record.reserve(read.getLengthTrue() * 2 + read.get_id_first().size() + 8);
        record += read.get_id_tab('1');
        record += '\t';
        record += read.get_seq_rc();
        record += '\t';
        record += read.get_qual_rc();
        if (read.get_comment().size() > 0){
            record += '\t';
            record += strjoin(read.get_comment(), "|");
        }
        record += '\n';
        output->write(record);
    }
};

class WriterHelper : public ReadVisitor {
public:
    virtual ~WriterHelper() {}
    WriterHelper(std::shared_ptr<OutputWriter> pe_, std::shared_ptr<OutputWriter> se_,
                 bool stranded_ = false, bool no_orphans_ = false) :
        stranded(stranded_), no_orphans(no_orphans_), pe(pe_), se(se_) {}

    void operator() (ReadBase &read) {
        read.accept(*this);
    }

    void operator() (SingleEndRead &read) {
        visit(&read);
    }

    void operator() (PairedEndRead &read) {
        visit(&read);
    }

    virtual void visit(PairedEndRead *per) {
        Read &one = per->non_const_read_one();
        Read &two = per->non_const_read_two();

        if (!one.getDiscard() && !two.getDiscard()) {
            pe->write(*per);
        } else if (!one.getDiscard() && !no_orphans) { // Will never be RC
            one.join_comment(two.get_comment());
            se->write_read(one, false);
        } else if (!two.getDiscard() && !no_orphans) { // if stranded RC
            two.join_comment(one.get_comment());
            se->write_read(two, stranded);
        }

    }

    virtual void visit(SingleEndRead *ser) {
        if (! (ser->non_const_read_one()).getDiscard() ) {
            se->write(*ser);
        }
    }

private:
    bool stranded;
    bool no_orphans;
    std::shared_ptr<OutputWriter> pe;
    std::shared_ptr<OutputWriter> se;
};

#endif
