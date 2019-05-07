#ifndef PRIMERS_H
#define PRIMERS_H
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

#ifndef MIN3
# define MIN3(a, b, c) ((a) < (b) ? ((a) < (c) ? (a) : (c)) : ((b) < (c) ? (b) : (c)))
#endif

class PrimerCounters : public Counters {

public:

    PrimerCounters(const std::string &statsFile, bool force, bool appendStats, const std::string &program_name, const std::string &notes) : Counters::Counters(statsFile, force, appendStats, program_name, notes) {
    }

    using Counters::input;
    virtual void input(const ReadBase &read) {
        Counters::input(read);
    }

    virtual void write_out() {

        initialize_json();

        write_labels(generic);
        write_sublabels("Single_end", se);
        write_sublabels("Paired_end", pe);

        finalize_json();
    }
};

typedef struct AlignPos {
    int dist; // The edit distance
    size_t spos; // Matching Start Position
    size_t epos; // Matching End Position
}ALIGNPOS;

/*
compute the Levenstein distance between a (primer) and b (sequence read)
pegged to the 5' end and bounded by edit distance k

float_bp is the pre-primer allowed float_bp basepair
max error is the max number of mismaches between primer and seq, not including primer_end_mismatches
end matches is the number of required perfect match basepairs at the end of primer
*/
ALIGNPOS
bounded_edit_distance(const std::string &primer, const std::string &seq, size_t float_bp, int max_error, size_t end_matches)
{
    size_t lastdiag, olddiag, cmin, endmatchcount;
    size_t primerlen = primer.length();
    size_t seqlen = seq.length();
    size_t column[primerlen - end_matches + 1];

    ALIGNPOS val = { max_error+1, 0, primerlen}; // dist and positions

    /* (primer) should always be greater than end_matches */
    if (primerlen < end_matches){
        val.dist = -2;
        return (val);
    }
    /* (primer) should always be < (read) */

    if (primerlen > seqlen) { // primer should never be greater than the seq
        val.dist = -2;
        return (val);
    }

    for (size_t x = 0; x <= float_bp; x++) {
        for (size_t i = 1; i <= primerlen - end_matches; i++)
            column[i] = i;
        for (size_t i = 1; i <= primerlen - end_matches + max_error ; i++) { // outer loop is the read
            column[0] = i;
            cmin = i;
            for (size_t j = 1, lastdiag = i-1; j <= primerlen - end_matches ; j++) { // inner loop is the primer
                olddiag = column[j];
                column[j] = MIN3(column[j] + 1, column[j-1] + 1, lastdiag + (primer[j-1] == seq[x+i-1] ? 0 : 1));
                lastdiag = olddiag;
                if (column[j] < cmin) cmin = column[j];
            }
            if (cmin > max_error) break; // if the smallest value in the column is > max error break
            if (column[primerlen - end_matches] <= val.dist ){
                endmatchcount=0;
                for (size_t l = 1; l <= end_matches; l++){
                    if (primer[primerlen - l] != seq[x + i + end_matches - l]){
                        break;
                    }
                    endmatchcount++;
                }
                if (endmatchcount == end_matches){

                    val.dist = column[primerlen - end_matches ]; // bottom right node, global alignment
                    val.spos = x;
                    val.epos = x + i + end_matches;
                }
            }
        }
    }
    return (val);
}

void check_read_pe(PairedEndRead &pe, const int pMismatches, const size_t pEndMismatches, const size_t pfloat, const size_t flip, const size_t keep) {

    ALIGNPOS best_val;
    Read &r1 = pe.non_const_read_one();
    Read &r2 = pe.non_const_read_two();

    const std::string &seq1 = r1.get_seq();
    const std::string &seq2 = r2.get_seq();

    const std::string p5Primer;
    best_val = bounded_edit_distance(p5Primer,  seq1,  pfloat,  pMismatches, pEndMismatches);
}

void check_read_se(SingleEndRead &se, const int pMismatches, const size_t pEndMismatches, const size_t pfloat, const size_t flip, const size_t keep) {

    Read &r1 = se.non_const_read_one();

    const std::string &seq1 = r1.get_seq();

    const std::string &seq2 = r1.get_seq_rc();

}

/* This is the helper class for Primer
 * The idea is ...
 * */
template <class T, class Impl>
void helper_Primers(InputReader<T, Impl> &reader, std::shared_ptr<OutputWriter> pe, std::shared_ptr<OutputWriter> se, PrimerCounters &counter, po::variables_map vm) {

    const int pMismatches = vm["primer_mismatches"].as<int>();
    const size_t pEndMismatches = vm["primer_end_mismatches"].as<size_t>();
    const size_t pfloat = vm["float"].as<size_t>();
    const size_t flip = vm["flip"].as<size_t>();
    const size_t keep = vm["keep"].as<size_t>();

    while(reader.has_next()) {
        auto i = reader.next();
        PairedEndRead* per = dynamic_cast<PairedEndRead*>(i.get());
        if (per) {
            counter.input(*per);
            counter.output(*per);
            writer_helper(per, pe, se);
        } else {
            SingleEndRead* ser = dynamic_cast<SingleEndRead*>(i.get());
            if (ser) {
                counter.input(*ser);
                counter.output(*ser);
                writer_helper(ser, pe, se);
            } else {
                throw std::runtime_error("Unknown read type");
            }
        }
    }
}

#endif
