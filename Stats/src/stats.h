#ifndef STATS_H
#define STATS_H
//  this is so we can implment hash function for dynamic_bitset
#define BOOST_DYNAMIC_BITSET_DONT_USE_FRIENDS

#include "ioHandler.h"
#include <map>
#include <unordered_map>
#include <boost/dynamic_bitset.hpp>
#include <boost/functional/hash.hpp>
#include <algorithm>
#include "utils.h"

extern template class InputReader<SingleEndRead, SingleEndReadFastqImpl>;
extern template class InputReader<PairedEndRead, PairedEndReadFastqImpl>;
extern template class InputReader<PairedEndRead, InterReadImpl>;
extern template class InputReader<ReadBase, TabReadImpl>;

class StatsCounters : public Counters {

public:
    std::vector<uint_fast64_t> R1_Length;
    std::vector<uint_fast64_t> R2_Length;
    std::vector<uint_fast64_t> SE_Length;

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

    StatsCounters(const std::string &statsFile, bool appendStats, const std::string &program_name, const std::string &notes) : Counters::Counters(statsFile, appendStats, program_name, notes) {
        R1_Length.resize(1);
        R2_Length.resize(1);
        SE_Length.resize(1);

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


    void read_stats(Read &r) {
        for (auto bp : r.get_seq()) {
            switch (bp) {
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
    void output(PairedEndRead &per) {
        Counters::output(per);
        Read &one = per.non_const_read_one();
        Read &two = per.non_const_read_two();
        read_stats(one);
        read_stats(two);
        uint64_t r1_q30bases=0;
        for (auto q : one.get_qual()) {
            r1_q30bases += (q - 33) >= 30;
        }
        uint64_t r2_q30bases=0;
        for (auto q : two.get_qual()) {
            r2_q30bases += (q - 33) >= 30;
        }
        R1_bQ30 += r1_q30bases;
        R2_bQ30 += r2_q30bases;
        R1_BpLen += one.getLength();
        R2_BpLen += two.getLength();
        if ( one.getLength() + 1 > R1_Length.size() ) {
            R1_Length.resize(one.getLength() + 1);
        }
        ++R1_Length[one.getLength()];
        if ( two.getLength() + 1 > R2_Length.size() ) {
            R2_Length.resize(two.getLength() + 1);
        }
        ++R2_Length[two.getLength()];
    }
    
    void output(SingleEndRead &ser) {
        Counters::output(ser);
        Read &one = ser.non_const_read_one();
        read_stats(one);
        uint_fast64_t q30bases=0;
        for (auto q : one.get_qual()) {
            q30bases += (q - 33) >= 30;
        }
        SE_bQ30 += q30bases;
        SE_BpLen += one.getLength();
        if ( one.getLength() + 1 > SE_Length.size() ) {
            SE_Length.resize(one.getLength() + 1);
        }
        ++SE_Length[one.getLength()];
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

        initialize_json();

        write_labels(generic);
        write_vector("SE_readlength_histogram",iSE_Length);
        write_vector("R1_readlength_histogram",iR1_Length);
        write_vector("R2_readlength_histogram",iR2_Length);
        write_sublabels("Base_composition", bases);
        write_sublabels("Single_end", se);
        write_sublabels("Paired_end", pe);

        finalize_json();        
    }
 
};

template <class T, class Impl>
void helper_stats(InputReader<T, Impl> &reader, std::shared_ptr<OutputWriter> pe, std::shared_ptr<OutputWriter> se, StatsCounters& counters) {
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

#endif
