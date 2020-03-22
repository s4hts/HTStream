#ifndef STATS_H
#define STATS_H
//  this is so we can implment hash function for dynamic_bitset
#define BOOST_DYNAMIC_BITSET_DONT_USE_FRIENDS

#include "ioHandler.h"
#include "utils.h"
#include "main_template.h"

#include <map>
#include <unordered_map>
#include <boost/dynamic_bitset.hpp>
#include <boost/functional/hash.hpp>
#include <algorithm>

extern template class InputReader<SingleEndRead, SingleEndReadFastqImpl>;
extern template class InputReader<PairedEndRead, PairedEndReadFastqImpl>;
extern template class InputReader<PairedEndRead, InterReadImpl>;
extern template class InputReader<ReadBase, TabReadImpl>;

class StatsCounters : public Counters {

public:
    Vec R1_Length;
    Vec R2_Length;
    Vec SE_Length;

    Mat R1_bases;
    Mat R2_bases;
    Mat SE_bases;

    Mat R1_qualities;
    Mat R2_qualities;
    Mat SE_qualities;

    std::vector<Label> bases;

    uint64_t A = 0;
    uint64_t C = 0;
    uint64_t G = 0;
    uint64_t T = 0;
    uint64_t N = 0;

    uint64_t SE_BpLen = 0;
    uint64_t SE_bQ30 = 0;

    uint64_t R1_BpLen = 0;
    uint64_t R1_bQ30 = 0;
    uint64_t R2_BpLen = 0;
    uint64_t R2_bQ30 = 0;

    StatsCounters(const std::string &statsFile, bool force, bool appendStats, const std::string &program_name, const std::string &notes) : Counters::Counters(statsFile, force, appendStats, program_name, notes) {
        se.push_back(std::forward_as_tuple("SE_bpLen", SE_BpLen));
        se.push_back(std::forward_as_tuple("SE_bQ30", SE_bQ30));

        pe.push_back(std::forward_as_tuple("R1_bpLen", R1_BpLen));
        pe.push_back(std::forward_as_tuple("R1_bQ30", R1_bQ30));
        pe.push_back(std::forward_as_tuple("R2_bpLen", R2_BpLen));
        pe.push_back(std::forward_as_tuple("R2_bQ30", R2_bQ30));

        bases.push_back(std::forward_as_tuple("A", A));
        bases.push_back(std::forward_as_tuple("C", C));
        bases.push_back(std::forward_as_tuple("G", G));
        bases.push_back(std::forward_as_tuple("T", T));
        bases.push_back(std::forward_as_tuple("N", N));
    }
    virtual ~StatsCounters() {}

    void read_stats(Read &r) {
        std::string seq = r.get_seq();
        size_t i = 0;

        for (std::string::iterator bp = seq.begin(); bp != seq.end(); ++bp) {
            i = static_cast<size_t>(bp - seq.begin() );

            switch (*bp) {
            case 'A':
                ++A;
                break;
            case 'C':
                ++C;
                break;
            case 'G':
                ++G;
                break;
            case 'T':
                ++T;
                break;
            case 'N':
                ++N;
                break;
            default:
                throw std::runtime_error("Unknown bp in stats counter");
            }
        }
    }

    using Counters::output;
    virtual void output(PairedEndRead &per, bool no_orphans = false) {
        size_t index = 0;

        Counters::output(per, no_orphans);
        Read &one = per.non_const_read_one();
        Read &two = per.non_const_read_two();

        // Total length of bases per read
        R1_BpLen += one.getLength();
        R2_BpLen += two.getLength();
        // Size histogram per read
        if ( one.getLength() + 1 > R1_Length.size() ) {
            R1_Length.resize(one.getLength() + 1);
        }
        ++R1_Length[one.getLength()];

        if ( two.getLength() + 1 > R2_Length.size() ) {
            R2_Length.resize(two.getLength() + 1);
        }
        ++R2_Length[two.getLength()];

        // READ 1 Base and Quality stats
        // update size of base and Q score matrix if needed
        size_t current_size = R1_bases.size();
        for( size_t gap = 0 ; gap < (one.getLength() - current_size); gap++ ) {
            Vec bases(5,0); // A,C,T,G,N
            Vec qualities(43,0); // quality score 0 to 42
            R1_bases.push_back(bases);
            R1_qualities.push_back(qualities);
        }
        std::string seq = one.get_seq();
        std::string qual = one.get_qual();
        uint64_t r1_q30bases=0;
        for (size_t index = 0; index < one.getLength(); ++index) {
            // bases
            char bp = seq[index];
            switch (bp) {
              case 'A':
                  ++A;
                  ++R1_bases[index][0];
                  break;
              case 'C':
                  ++C;
                  ++R1_bases[index][1];
                  break;
              case 'G':
                  ++G;
                  ++R1_bases[index][2];
                  break;
              case 'T':
                  ++T;
                  ++R1_bases[index][3];
                  break;
              case 'N':
                  ++N;
                  ++R1_bases[index][4];
                  break;
              default:
                  throw std::runtime_error("Unknown bp in stats counter");
            }
            // qualities
            size_t qscore = qual[index];
            r1_q30bases += (qscore - 33) >= 30;
            ++R1_qualities[index][(qscore - 33)];
        }
        R1_bQ30 += r1_q30bases;

        // READ 2 Base and Quality stats
        // update size of base and Q score matrix if needed
        current_size = R2_bases.size();
        for( size_t gap = 0 ; gap < (two.getLength() - current_size); gap++ ) {
            Vec bases(5,0); // A,C,T,G,N
            Vec qualities(43,0); // quality score 0 to 42
            R2_bases.push_back(bases);
            R2_qualities.push_back(qualities);
        }
        seq = two.get_seq();
        qual = two.get_qual();
        uint64_t r2_q30bases=0;
        for (size_t index = 0; index < two.getLength(); ++index) {
            // bases
            char bp = seq[index];
            switch (bp) {
              case 'A':
                  ++A;
                  ++R2_bases[index][0];
                  break;
              case 'C':
                  ++C;
                  ++R2_bases[index][1];
                  break;
              case 'G':
                  ++G;
                  ++R2_bases[index][2];
                  break;
              case 'T':
                  ++T;
                  ++R2_bases[index][3];
                  break;
              case 'N':
                  ++N;
                  ++R2_bases[index][4];
                  break;
              default:
                  throw std::runtime_error("Unknown bp in stats counter");
            }
            // qualities
            size_t qscore = qual[index];
            r2_q30bases += (qscore - 33) >= 30;
            ++R2_qualities[index][(qscore - 33)];
        }
        R2_bQ30 += r2_q30bases;
    }

    void output(SingleEndRead &ser) {
        Counters::output(ser);
        Read &one = ser.non_const_read_one();

        SE_BpLen += one.getLength();
        if ( one.getLength() + 1 > SE_Length.size() ) {
            SE_Length.resize(one.getLength() + 1);
        }
        ++SE_Length[one.getLength()];
        // Single end Base and Quality stats
        // update size of base and Q score matrix if needed
        size_t current_size = SE_bases.size();
        for( size_t gap = 0 ; gap < (one.getLength() - current_size); gap++ ) {
            Vec bases(5,0); // A,C,T,G,N
            Vec qualities(43,0); // quality score 0 to 42
            SE_bases.push_back(bases);
            SE_qualities.push_back(qualities);
        }
        std::string seq = one.get_seq();
        std::string qual = one.get_qual();
        uint64_t q30bases=0;
        for (size_t index = 0; index < one.getLength(); ++index) {
            // bases
            char bp = seq[index];
            switch (bp) {
              case 'A':
                  ++A;
                  ++SE_bases[index][0];
                  break;
              case 'C':
                  ++C;
                  ++SE_bases[index][1];
                  break;
              case 'G':
                  ++G;
                  ++SE_bases[index][2];
                  break;
              case 'T':
                  ++T;
                  ++SE_bases[index][3];
                  break;
              case 'N':
                  ++N;
                  ++SE_bases[index][4];
                  break;
              default:
                  throw std::runtime_error("Unknown bp in stats counter");
            }
            // qualities
            size_t qscore = qual[index];
            q30bases += (qscore - 33) >= 30;
            ++SE_qualities[index][(qscore - 33)];
        }
        SE_bQ30 += q30bases;
    }

    virtual void write_out() {
        std::vector<Vector> iSE_Length;
        for (size_t i = 1; i < SE_Length.size(); ++i) {
            if (SE_Length[i] > 0) {
                iSE_Length.push_back(std::forward_as_tuple(i, SE_Length[i]));
            }
        }

        std::vector<Vector> iR1_Length;
        for (size_t i = 1; i < R1_Length.size(); ++i) {
            if (R1_Length[i] > 0) {
                iR1_Length.push_back(std::forward_as_tuple(i, R1_Length[i]));
            }
        }

        std::vector<Vector> iR2_Length;
        for (size_t i = 1; i < R2_Length.size(); ++i) {
            if (R2_Length[i] > 0) {
                iR2_Length.push_back(std::forward_as_tuple(i, R2_Length[i]));
            }
        }

        std::vector<std::string> ind_se;
        for (size_t j = 1; j <= SE_bases.size(); j++){
          ind_se.push_back(std::to_string((int)j));
        }
        std::vector<std::string> ind_pe1;
        for (size_t j = 1; j <= R1_bases.size(); j++){
          ind_pe1.push_back(std::to_string((int)j));
        }
        std::vector<std::string> ind_pe2;
        for (size_t j = 1; j <= R2_bases.size(); j++){
          ind_pe2.push_back(std::to_string((int)j));
        }
        std::vector<std::string> b{ "A", "C", "G", "T", "N"};
        std::vector<std::string> q;
        for (size_t j = 0; j < 43; j++){
          q.push_back(std::to_string((int)j));
        }

        initialize_json();

        write_labels(generic);
        write_sublabels("Base_composition", bases);
        start_sublabel("Single_end");
        write_labels(se, 2);
        write_vector("SE_readlength_histogram",iSE_Length, 2);
        write_matrix("SE_Base_by_Cycle",SE_bases, b, ind_se, 0, 2);
        write_matrix("SE_Qualities_by_Cycle",SE_qualities, q, ind_se, 0, 2);
        end_sublabel();
        start_sublabel("Paired_end");
        write_labels(pe, 2);
        write_vector("R1_readlength_histogram",iR1_Length, 2);
        write_matrix("R1_Base_by_Cycle",R1_bases, b, ind_pe1, 0, 2);
        write_matrix("R1_Qualities_by_Cycle",R1_qualities, q, ind_pe1, 0, 2);
        write_vector("R2_readlength_histogram",iR2_Length, 2);
        write_matrix("R2_Base_by_Cycle",R2_bases, b, ind_pe2, 0, 2);
        write_matrix("R2_Qualities_by_Cycle",R2_qualities, q, ind_pe2, 0, 2);
        end_sublabel();

        finalize_json();
    }

};

class Stats: public MainTemplate<StatsCounters, Stats> {
public:

    Stats() {
        program_name = "hts_Stats";
        app_description =
            "The hts_Stats app produce basic statistics about the reads in a dataset.\n";
        app_description += "  Including the basepair composition and number of bases Q30.";
    }

    void add_extra_options(po::options_description &) {
    }

    template <class T, class Impl>
    void do_app(InputReader<T, Impl> &reader, std::shared_ptr<OutputWriter> pe, std::shared_ptr<OutputWriter> se, StatsCounters& counters, const po::variables_map &) {
        while(reader.has_next()) {
            auto i = reader.next();
            PairedEndRead* per = dynamic_cast<PairedEndRead*>(i.get());
            if (per) {
                counters.input(*per);
                counters.output(*per);
                writer_helper(per, pe, se, false);
            } else {
                SingleEndRead* ser = dynamic_cast<SingleEndRead*>(i.get());
                if (ser) {
                    counters.input(*ser);
                    counters.output(*ser);
                    writer_helper(ser, pe, se, false);
                } else {
                    throw std::runtime_error("Unknown read type");
                }
            }
        }
    }
};
#endif
