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

typedef std::unordered_map <std::string, std::string> SeqMap;

class PrimerCounters : public Counters {

public:

    uint64_t keep_primer = 0;
    uint64_t flipped = 0;

    uint64_t SE_Discarded = 0;

    uint64_t R1_Discarded = 0;
    uint64_t R2_Discarded = 0;
    uint64_t PE_Discarded = 0;

    PrimerCounters(const std::string &statsFile, bool force, bool appendStats, const std::string &program_name, const std::string &notes) : Counters::Counters(statsFile, force, appendStats, program_name, notes) {
        generic.push_back(std::forward_as_tuple("Keep_primer", keep_primer));
        generic.push_back(std::forward_as_tuple("Flipped", flipped));

        se.push_back(std::forward_as_tuple("SE_discarded", SE_Discarded));

        pe.push_back(std::forward_as_tuple("R1_discarded", R1_Discarded));
        pe.push_back(std::forward_as_tuple("R2_discarded", R2_Discarded));
        pe.push_back(std::forward_as_tuple("PE_discarded", PE_Discarded));
    }

    void set_keep_primer(bool keep) {
        keep_primer = (uint64_t)keep;
    }
    void increment_flipped() {
        flipped++;
    }

    using Counters::input;
    virtual void input(const ReadBase &read) {
        Counters::input(read);
    }
    virtual void output(SingleEndRead &ser)  {
        if (ser.non_const_read_one().getDiscard()) {
            ++SE_Discarded;
        } else {
            ++SE_Out;
            ++TotalFragmentsOutput;
        }
    }

    virtual void output(PairedEndRead &per, bool no_orphans = false)  {
        Read &one = per.non_const_read_one();
        Read &two = per.non_const_read_two();
        if (!one.getDiscard() && !two.getDiscard()) {
            ++PE_Out;
            ++TotalFragmentsOutput;
        } else if (!one.getDiscard() && !no_orphans) { //if stranded RC
            ++SE_Out;
            ++R2_Discarded;
            ++TotalFragmentsOutput;
        } else if (!two.getDiscard() && !no_orphans) { // Will never be RC
            ++SE_Out;
            ++R1_Discarded;
            ++TotalFragmentsOutput;
        } else {
            ++PE_Discarded;
        }
    }

    virtual void write_out() {

        initialize_json();

        write_labels(generic);
        write_sublabels("Single_end", se);
        write_sublabels("Paired_end", pe);

        finalize_json();
    }
};

SeqMap fasta2dict(std::string primers){
    SeqMap primerMap;
    bf::path p(primers);
    if (bf::exists(p)) {
      // fastq file
      bi::stream <bi::file_descriptor_source> fa_to_read{check_open_r(primers), bi::close_handle};
      InputReader<SingleEndRead, FastaReadImpl> fp(fa_to_read);
      while(fp.has_next()) {
          auto i = fp.next();
          SingleEndRead* ser = dynamic_cast<SingleEndRead*>(i.get());
          Read r1 = ser->non_const_read_one();
          primerMap[r1.get_id()] = r1.get_seq();
      }
    } else {
      // comma seperated
      std::istringstream fa_to_read(string2fasta(primers));
      InputReader<SingleEndRead, FastaReadImpl> fp(fa_to_read);
      while(fp.has_next()) {
          auto i = fp.next();
          SingleEndRead* ser = dynamic_cast<SingleEndRead*>(i.get());
          Read r1 = ser->non_const_read_one();
          primerMap[r1.get_id()] = r1.get_seq();
      }
    }
    return (primerMap);
}

/*
Ambiguity matching, return 0 on match and 1 on no-match
A |	A	T
C |	C	G
G |	G	C
M |	A or C
R |	A or G
W |	A or T
S |	C or G
Y |	C or T
K |	G or T
V |	A or C or G
H |	A or C or T
D |	A or G or T
B |	C or G or T
N |	G or A or T or C
*/
int charMatch(const char p, const char s){
    if (p == s){
      return 0;
    } else {
      switch(p) {
        case 'M' : if (s == 'A' || s == 'C' ) return  0;
                 break;
        case 'R' : if (s == 'A' || s == 'G' ) return  0;
                 break;
        case 'W' : if (s == 'A' || s == 'T' ) return  0;
                 break;
        case 'S' : if (s == 'C' || s == 'G' ) return  0;
                 break;
        case 'Y' : if (s == 'C' || s == 'T' ) return  0;
                 break;
        case 'K' : if (s == 'G' || s == 'T' ) return  0;
                 break;
        case 'V' : if (s == 'A' || s == 'C' || s == 'G' ) return  0;
                 break;
        case 'H' : if (s == 'A' || s == 'C' || s == 'T' ) return  0;
                 break;
        case 'D' : if (s == 'A' || s == 'G' || s == 'T' ) return  0;
                 break;
        case 'B' : if (s == 'C' || s == 'G' || s == 'T' ) return  0;
                 break;
        case 'N' : if (s == 'A' || s == 'C' || s == 'G' || s == 'T' ) return  0;
                 break;
      }
    }
    return 1;
}

typedef struct AlignPos {
    int dist; // The edit distance
    size_t spos; // Matching Start Position
    size_t epos; // Matching End Position
    std::string name;
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
                column[j] = MIN3(column[j] + 1, column[j-1] + 1, lastdiag + charMatch(primer[j-1],seq[x+i-1]));
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

void check_read_pe(PairedEndRead &pe, PrimerCounters &counter, SeqMap &primer5p, SeqMap &primer3p, const int pMismatches, const size_t pEndMismatches, const size_t pfloat, const size_t flip, const size_t keep) {

    ALIGNPOS test_val, best_val;
    Read &r1 = pe.non_const_read_one();
    Read &r2 = pe.non_const_read_two();

    best_val.dist = pMismatches + 1;
    const std::string &seq1 = r1.get_seq();
    for ( auto it = primer5p.begin(); it != primer5p.end(); ++it ){
      const std::string p5Primer = it->second;
      test_val = bounded_edit_distance(p5Primer,  seq1,  pfloat,  pMismatches, pEndMismatches);
      if (test_val.dist < best_val.dist){
        best_val = test_val;
        best_val.name = it->first;
      }
      if (best_val.dist == 0) break;
    }
    if (best_val.dist <= pMismatches){
      r1.add_comment("P5:Z:" + best_val.name);
      if (!keep) r1.setLCut(best_val.epos);
    } else if (flip) {
      best_val.dist = pMismatches + 1;
      const std::string &seq1 = r2.get_seq();
      for ( auto it = primer5p.begin(); it != primer5p.end(); ++it ){
        const std::string p5Primer = it->second;
        test_val = bounded_edit_distance(p5Primer,  seq1,  pfloat,  pMismatches, pEndMismatches);
        if (test_val.dist < best_val.dist){
          best_val = test_val;
          best_val.name = it->first;
        }
        if (best_val.dist == 0) break;
      }
      if (best_val.dist <= pMismatches){
        std::swap(r1, r2);
        counter.increment_flipped();
        r1.add_comment("P5:Z:" + best_val.name + "-FLIP");
        if (!keep) r1.setLCut(best_val.epos);
      }
    }
    best_val.dist = pMismatches + 1;
    const std::string &seq2 = r2.get_seq();
    for ( auto it = primer3p.begin(); it != primer3p.end(); ++it ){
      const std::string p3Primer = it->second;
      test_val = bounded_edit_distance(p3Primer,  seq2,  pfloat,  pMismatches, pEndMismatches);
      if (test_val.dist < best_val.dist){
        best_val = test_val;
        best_val.name = it->first;
      }
      if (best_val.dist == 0) break;
    }
    if (best_val.dist <= pMismatches){
      r2.add_comment("P3:Z:" + best_val.name);
      if (!keep) r2.setLCut(best_val.epos);
    }

}

void check_read_se(SingleEndRead &se, PrimerCounters &counter, SeqMap &primer5p, SeqMap &primer3p, const int pMismatches, const size_t pEndMismatches, const size_t pfloat, const size_t flip, const size_t keep) {

    ALIGNPOS test_val, best_val;
    Read &r1 = se.non_const_read_one();

    best_val.dist = pMismatches + 1;
    const std::string &seq1 = r1.get_seq();
    for ( auto it = primer5p.begin(); it != primer5p.end(); ++it ){
      const std::string p5Primer = it->second;
      test_val = bounded_edit_distance(p5Primer,  seq1,  pfloat,  pMismatches, pEndMismatches);
      if (test_val.dist < best_val.dist){
        best_val = test_val;
        best_val.name = it->first;
      }
      if (best_val.dist == 0) break;
    }
    if (best_val.dist <= pMismatches){
      r1.add_comment("P5:Z:" + best_val.name);
      if (!keep) r1.setLCut(best_val.epos);
    } else if (flip) {
      Read &temp = r1;
      temp.set_read_rc();
      best_val.dist = pMismatches + 1;
      const std::string &seq1 = temp.get_seq();
      for ( auto it = primer5p.begin(); it != primer5p.end(); ++it ){
        const std::string p5Primer = it->second;
        test_val = bounded_edit_distance(p5Primer,  seq1,  pfloat,  pMismatches, pEndMismatches);
        if (test_val.dist < best_val.dist){
          best_val = test_val;
          best_val.name = it->first;
        }
        if (best_val.dist == 0) break;
      }
      if (best_val.dist <= pMismatches){
        r1 = temp;
        counter.increment_flipped();
        r1.add_comment("P5:Z:" + best_val.name + "-FLIP");
        if (!keep) r1.setLCut(best_val.epos);
      }
    }
    best_val.dist = pMismatches + 1;
    const std::string &seq2 = r1.get_seq_rc();
    for ( auto it = primer3p.begin(); it != primer3p.end(); ++it ){
      const std::string p3Primer = it->second;
      test_val = bounded_edit_distance(p3Primer,  seq2,  pfloat,  pMismatches, pEndMismatches);
      if (test_val.dist < best_val.dist){
        best_val = test_val;
        best_val.name = it->first;
      }
      if (best_val.dist == 0) break;
    }
    if (best_val.dist <= pMismatches){
      r1.add_comment("P3:Z:" + best_val.name);
      if (!keep) r1.setRCut(r1.getLength() -  best_val.epos);
    }
}

/* This is the helper class for Primer
 * The idea is ...
 * */
template <class T, class Impl>
void helper_Primers(InputReader<T, Impl> &reader, std::shared_ptr<OutputWriter> pe, std::shared_ptr<OutputWriter> se, PrimerCounters &counter, po::variables_map vm, SeqMap &primer5p, SeqMap &primer3p) {

    const int pMismatches = vm["primer_mismatches"].as<int>();
    const size_t pEndMismatches = vm["primer_end_mismatches"].as<size_t>();
    const size_t pfloat = vm["float"].as<size_t>();
    const bool flip = vm["flip"].as<bool>();
    const bool keep = vm["keep"].as<bool>();

    while(reader.has_next()) {
        auto i = reader.next();
        PairedEndRead* per = dynamic_cast<PairedEndRead*>(i.get());
        if (per) {
            counter.input(*per);
            check_read_pe(*per, counter, primer5p, primer3p, pMismatches, pEndMismatches, pfloat, flip, keep);
            per->checkDiscarded(vm["min-length"].as<size_t>());
            counter.output(*per, vm["no-orphans"].as<bool>());
            writer_helper(per, pe, se, vm["stranded"].as<bool>(), vm["no-orphans"].as<bool>());
        } else {
            SingleEndRead* ser = dynamic_cast<SingleEndRead*>(i.get());
            if (ser) {
                counter.input(*ser);
                check_read_se(*ser, counter, primer5p, primer3p, pMismatches, pEndMismatches, pfloat, flip, keep);
                ser->checkDiscarded(vm["min-length"].as<size_t>());
                counter.output(*ser);
                writer_helper(ser, pe, se);
            } else {
                throw std::runtime_error("Unknown read type");
            }
        }
    }
}

#endif
