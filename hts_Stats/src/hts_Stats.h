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

const uint64_t QUAL_MAX = 94;

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

    uint64_t SE_bQ30 = 0;

    uint64_t R1_bQ30 = 0;
    uint64_t R2_bQ30 = 0;

    StatsCounters(const std::string &program_name, const po::variables_map &vm) : Counters::Counters(program_name, vm) {
        R1_Length.resize(1);
        R2_Length.resize(1);
        SE_Length.resize(1);
        se.push_back(std::forward_as_tuple("total_Q30_basepairs", SE_bQ30));
        r1.push_back(std::forward_as_tuple("total_Q30_basepairs", R1_bQ30));
        r2.push_back(std::forward_as_tuple("total_Q30_basepairs", R2_bQ30));

        bases.push_back(std::forward_as_tuple("A", A));
        bases.push_back(std::forward_as_tuple("C", C));
        bases.push_back(std::forward_as_tuple("G", G));
        bases.push_back(std::forward_as_tuple("T", T));
        bases.push_back(std::forward_as_tuple("N", N));
    }
    virtual ~StatsCounters() {}

    void read_stats(Read &r, Vec &Length, Mat &read_bases, Mat &read_qualities, uint64_t &read_bQ30) {
        // Size histogram per read
        if ( r.getLength() + 1 > Length.size() ) {
            Length.resize(r.getLength() + 1);
        }
        ++Length[r.getLength()];
        // READ Base and Quality stats
        // update size of base and Q score matrix if needed
        for( size_t gap = 0 ; read_bases.size() < r.getLength(); gap++ ) {
            Vec bases(5,0); // A,C,T,G,N
            Vec qualities(QUAL_MAX,0); // quality score 0 to MAX
            read_bases.push_back(bases);
            read_qualities.push_back(qualities);
        }
        std::string seq = r.get_seq();
        std::string qual = r.get_qual();
        uint64_t q30bases=0;
        for (size_t index = 0; index < r.getLength(); ++index) {
            // bases
            char bp = seq[index];
            switch (bp) {
              case 'A':
                  ++A;
                  ++read_bases[index][0];
                  break;
              case 'C':
                  ++C;
                  ++read_bases[index][1];
                  break;
              case 'G':
                  ++G;
                  ++read_bases[index][2];
                  break;
              case 'T':
                  ++T;
                  ++read_bases[index][3];
                  break;
              case 'N':
                  ++N;
                  ++read_bases[index][4];
                  break;
              default:
                  throw HtsRuntimeException(std::string("Unknown bp in stats counter: ") + bp);
            }
            // qualities
            size_t qscore = qual[index];
            uint_fast64_t qscore_int = qscore - 33;
            if (qscore_int < QUAL_MAX) {
                q30bases += (qscore_int) >= 30;
                ++read_qualities[index][qscore_int];
            }
        }
        read_bQ30 += q30bases;
    }

    using Counters::output;
    void output(PairedEndRead &per) {
        Counters::output(per);
        read_stats(per.non_const_read_one(), R1_Length, R1_bases, R1_qualities, R1_bQ30);
        read_stats(per.non_const_read_two(), R2_Length, R2_bases, R2_qualities, R2_bQ30);
    }

    void output(SingleEndRead &ser) {
        Counters::output(ser);
        read_stats(ser.non_const_read_one(), SE_Length, SE_bases, SE_qualities, SE_bQ30);
    }

    void write_out() {
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
        for (size_t j = 0; j < QUAL_MAX; j++){
          q.push_back(std::to_string((int)j));
        }

        initialize_json();

        start_sublabel("Program_details");
        write_values(pd, 2);
        start_sublabel("options",2);
        write_options(3);
        end_sublabel(2);
        end_sublabel();

        start_sublabel("Fragment");
        write_values(fragment, 2);
        start_sublabel("base_composition", 2);
        write_values(bases, 3);
        end_sublabel(2);
        end_sublabel();

        start_sublabel("Single_end");
        write_values(se, 2);
        write_vector("readlength_histogram",iSE_Length, 2);
        write_matrix("base_by_cycle",SE_bases, b, ind_se, 0, 2);
        write_matrix("qualities_by_cycle",SE_qualities, q, ind_se, 0, 2);
        end_sublabel();

        start_sublabel("Paired_end");
        write_values(pe, 2);
        start_sublabel("Read1",2);
        write_values(r1, 3);
        write_vector("readlength_histogram",iR1_Length, 3);
        write_matrix("base_by_cycle",R1_bases, b, ind_pe1, 0, 3);
        write_matrix("qualities_by_cycle",R1_qualities, q, ind_pe1, 0, 3);
        end_sublabel(2);
        start_sublabel("Read2",2);
        write_values(r2, 3);
        write_vector("readlength_histogram",iR2_Length, 3);
        write_matrix("base_by_cycle",R2_bases, b, ind_pe2, 0, 3);
        write_matrix("qualities_by_cycle",R2_qualities, q, ind_pe2, 0, 3);
        end_sublabel(2);
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

        WriterHelper writer(pe, se, false);

        while(reader.has_next()) {
            auto i = reader.next();
            counters.input(*i);
            counters.output(*i);
            writer(*i);
        }
    }
};
#endif
