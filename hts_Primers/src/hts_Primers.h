#ifndef PRIMERS_H
#define PRIMERS_H
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

#ifndef MIN3
# define MIN3(a, b, c) ((a) < (b) ? ((a) < (c) ? (a) : (c)) : ((b) < (c) ? (b) : (c)))
#endif

typedef std::unordered_map <std::string, std::string> SeqMap;

class PrimerCounters : public Counters {

public:

    uint64_t keep_primer = 0;
    uint64_t flipped = 0;

    std::vector<sLabel> primers;
    std::unordered_map < std::string, uint_fast64_t> primers_seen_counter;

    uint64_t SE_Discarded = 0;
    uint64_t SE_Primer_Trim = 0;
    uint64_t SE_Primer_BpTrim = 0;

    uint64_t R1_Discarded = 0;
    uint64_t R1_Primer_Trim = 0;
    uint64_t R1_Primer_BpTrim = 0;
    uint64_t R2_Discarded = 0;
    uint64_t R2_Primer_Trim = 0;
    uint64_t R2_Primer_BpTrim = 0;
    uint64_t PE_Discarded = 0;

    PrimerCounters(const std::string &program_name, const po::variables_map &vm) : Counters::Counters(program_name, vm) {

        fragment.push_back(std::forward_as_tuple("flipped", flipped));

        se.push_back(std::forward_as_tuple("primerTrim", SE_Primer_Trim));
        se.push_back(std::forward_as_tuple("primerBpTrim", SE_Primer_BpTrim));
        se.push_back(std::forward_as_tuple("discarded", SE_Discarded));

        r1.push_back(std::forward_as_tuple("primerTrim", R1_Primer_Trim));
        r1.push_back(std::forward_as_tuple("primerBpTrim", R1_Primer_BpTrim));
        r1.push_back(std::forward_as_tuple("discarded", R1_Discarded));
        r2.push_back(std::forward_as_tuple("primerTrim", R2_Primer_Trim));
        r2.push_back(std::forward_as_tuple("primerBpTrim", R2_Primer_BpTrim));
        r2.push_back(std::forward_as_tuple("discarded", R2_Discarded));
        pe.push_back(std::forward_as_tuple("discarded", PE_Discarded));
    }

    void set_seqmap(SeqMap primerMap_5p, SeqMap primerMap_3p) {
        SeqMap primers_5p = primerMap_5p;
        for (auto &pt: primers_5p) {
            primers.push_back(std::forward_as_tuple(pt.first, pt.second));
        }
        SeqMap primers_3p = primerMap_3p;
        for (auto &pt: primers_3p) {
            primers.push_back(std::forward_as_tuple(pt.first, pt.second));
        }
    }

    void increment_flipped() {
        flipped++;
    }

    using Counters::input;
    virtual void input(const ReadBase &read) {
        Counters::input(read);
    }

    void primer_match_counter(std::string &p5primer, std::string &p3primer){
        std::string primerPair = "\"" + p5primer + "\",\"" + p3primer + "\"";
        primers_seen_counter[primerPair]++;
    }

    virtual void output(SingleEndRead &ser)  {
        if (ser.non_const_read_one().getDiscard()) {
            ++SE_Discarded;
        } else {
            Read &one = ser.non_const_read_one();
            if (one.getLengthTrue() < one.getLength()) {
                ++SE_Primer_Trim;
                SE_Primer_BpTrim += (one.getLength() - one.getLengthTrue());
            }
            ++SE_Out;
            ++TotalFragmentsOutput;
        }
    }

    void output(PairedEndRead &per, bool no_orphans = false)  {
        Read &one = per.non_const_read_one();
        Read &two = per.non_const_read_two();
        if (!one.getDiscard() && !two.getDiscard()) {
            if (one.getLengthTrue() < one.getLength()){
              ++R1_Primer_Trim;
              R1_Primer_BpTrim += (one.getLength() - one.getLengthTrue());
            }
            if (two.getLengthTrue() < two.getLength()) {
              ++R2_Primer_Trim;
              R2_Primer_BpTrim += (two.getLength() - two.getLengthTrue());
            }
            ++PE_Out;
            ++TotalFragmentsOutput;
        } else if (!one.getDiscard() && !no_orphans) { //if stranded RC
            if (one.getLengthTrue() < one.getLength()) {
                ++SE_Primer_Trim;
                SE_Primer_BpTrim += (one.getLength() - one.getLengthTrue());
            }
            ++SE_Out;
            ++R2_Discarded;
            ++TotalFragmentsOutput;
        } else if (!two.getDiscard() && !no_orphans) { // Will never be RC
            if (two.getLengthTrue() < two.getLength()) {
                ++SE_Primer_Trim;
                SE_Primer_BpTrim += (two.getLength() - two.getLengthTrue());
            }
            ++SE_Out;
            ++R1_Discarded;
            ++TotalFragmentsOutput;
        } else {
            ++PE_Discarded;
        }
    }

    virtual void write_out() {

        std::vector<Label> iPrimers;
        for (auto &it: primers_seen_counter) {
            if (it.second > 0) {
                iPrimers.push_back(std::forward_as_tuple(it.first, it.second));
            }
        }

        initialize_json();

        start_sublabel("Program_details");
        write_values(pd, 2);
        start_sublabel("options",2);
        write_options(3);
        end_sublabel(2);
        write_vector("primers",primers, 2);
        end_sublabel();

        start_sublabel("Fragment");
        write_values(fragment, 2);
        write_vector("primers_counts",iPrimers, 2);
        end_sublabel();

        start_sublabel("Single_end");
        write_values(se, 2);
        end_sublabel();

        start_sublabel("Paired_end");
        write_values(pe, 2);
        start_sublabel("Read1",2);
        write_values(r1, 3);
        end_sublabel(2);
        start_sublabel("Read2",2);
        write_values(r2, 3);
        end_sublabel(2);
        end_sublabel();

        finalize_json();
    }
};

typedef struct AlignPos {
    long dist; // The edit distance
    size_t spos; // Matching Start Position
    size_t epos; // Matching End Position
    std::string name;
}ALIGNPOS;

class Primers: public MainTemplate<PrimerCounters, Primers> {
public:

    Primers() {
        program_name = "hts_Primers";
        app_description =
            "The hts_Primers application identifies primer sequences located on the 5' ends of R1 and R2,\n";
        app_description += "    or 5' and 3' end of SE reads, optionally cut/flip and return the the read adding the \n";
        app_description += "    primer to the read id.\n";
    }

    void add_extra_options(po::options_description &desc) {
        setDefaultParamsCutting(desc);
        // no-orphans|n ; stranded|s ; min-length|m

        desc.add_options()
            ("primers_5p,P", po::value<std::string>(), "5' primers, comma separated list of sequences, or fasta file");
        desc.add_options()
            ("primers_3p,Q", po::value<std::string>(), "3' primers, comma separated list of sequences, or fasta file");
        desc.add_options()
            ("primer_mismatches,d", po::value<size_t>()->default_value(4)->notifier(boost::bind(&check_range<size_t>, "primer_mismatches", _1, 0, 10000)), "Max hamming dist from primer (min 0, max 10000)");
        desc.add_options()
                              ("primer_end_mismatches,e", po::value<size_t>()->default_value(4)->notifier(boost::bind(&check_range<size_t>, "primer_end_mismatches", _1, 0, 10000)), "Required number of matching bases at end of primer (min 0, max 10000)");
        desc.add_options()
            ("float,l", po::value<size_t>()->default_value(0)->notifier(boost::bind(&check_range<size_t>, "float", _1, 0, 10000)), "Variable number of bases preceeding primer allowed to float");
        desc.add_options()
            ("flip,x", po::bool_switch()->default_value(false), "Primers can be seen in both orientiations, tests flip and reorients all reads to the same orientation.");
        desc.add_options()
            ("keep,k", po::bool_switch()->default_value(false), "Don't cut off the primer sequence, leave it as a part of the read");
        desc.add_options()
            ("min_primer_matches,r", po::value<size_t>()->default_value(0)->notifier(boost::bind(&check_range<size_t>, "min_primer_matches", _1, 0, 2)), "Minimum number of primers to match to keep the fragment (0, keep all fragments, 1 must match either 5' or 3' primer, 2 must match both 5' and 3' primers)");
    }

    SeqMap fasta2dict(std::string primers, std::string prefix){
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
                primerMap[r1.get_id_orig()] = r1.get_seq();
            }
        } else {
            // comma seperated
            std::istringstream fa_to_read(string2fasta(primers, prefix));
            InputReader<SingleEndRead, FastaReadImpl> fp(fa_to_read);
            while(fp.has_next()) {
                auto i = fp.next();
                SingleEndRead* ser = dynamic_cast<SingleEndRead*>(i.get());
                Read r1 = ser->non_const_read_one();
                primerMap[r1.get_id_orig()] = r1.get_seq();
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
    int charMatch(const char pr, const char se){
        const char p = std::toupper(pr);
        const char s = std::toupper(se);

        if (p == s){
            return 0;
        } else { // ambiguity base
            switch(p) {
            case 'M' : if (s == 'M' || s == 'A' || s == 'C' ) return  0;
                break;
            case 'R' : if (s == 'R' || s == 'A' || s == 'G' ) return  0;
                break;
            case 'W' : if (s == 'W' || s == 'A' || s == 'T' ) return  0;
                break;
            case 'S' : if (s == 'S' || s == 'C' || s == 'G' ) return  0;
                break;
            case 'Y' : if (s == 'Y' || s == 'C' || s == 'T' ) return  0;
                break;
            case 'K' : if (s == 'K' || s == 'G' || s == 'T' ) return  0;
                break;
            case 'V' : if (s == 'V' || s == 'M' || s == 'R' || s == 'S' || s == 'A' || s == 'C' || s == 'G' ) return  0;
                break;
            case 'H' : if (s == 'H' || s == 'M' || s == 'W' || s == 'Y' || s == 'A' || s == 'C' || s == 'T' ) return  0;
                break;
            case 'D' : if (s == 'D' || s == 'R' || s == 'W' || s == 'K' || s == 'A' || s == 'G' || s == 'T' ) return  0;
                break;
            case 'B' : if (s == 'B' || s == 'S' || s == 'Y' || s == 'K' || s == 'C' || s == 'G' || s == 'T' ) return  0;
                break;
            case 'N' : return  0;
                break;
            }
        }
        return 1;
    }


/*
  compute the Levenstein distance between a (primer) and b (sequence read)
  pegged to the 5' end and bounded by edit distance k

  float_bp is the pre-primer allowed float_bp basepair
  max error is the max number of mismaches between primer and seq, not including primer_end_mismatches
  end matches is the number of required perfect match basepairs at the end of primer
*/
    ALIGNPOS
    bounded_edit_distance(const std::string &primer, const std::string &seq, size_t float_bp, size_t max_error, size_t end_matches)
    {
        long lastdiag, olddiag, cmin;
        size_t endmatchcount;
        (void)lastdiag;
        size_t primerlen = primer.length();
        size_t seqlen = seq.length();
        long column[primerlen - end_matches + 1];

        ALIGNPOS val = { long(max_error+1), 0u, primerlen, ""}; // dist and positions

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
                    column[j] = MIN3(column[j] + 1, column[j-1] + 1, long(lastdiag + charMatch(primer[j-1],seq[x+i-1])));
                    lastdiag = olddiag;
                    if (column[j] < cmin) cmin = column[j];
                }
                if (cmin > long(max_error)) break; // if the smallest value in the column is > max error break
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

    void check_read_pe(PairedEndRead &pe, PrimerCounters &counter, SeqMap &primer5p, SeqMap &primer3p, const size_t pMismatches, const size_t pEndMismatches, const size_t pfloat, const size_t flip, const size_t keep, const size_t mpmatches) {

        ALIGNPOS test_val, best_val;
        std::string p5primer = "None", p3primer = "None";
        size_t pmatches = 0;
        bool flipped = false;

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
        if (best_val.dist <= long(pMismatches)){
            p5primer = best_val.name;
            if (!keep) r1.setLCut(best_val.epos);
            pmatches++;
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
            if (best_val.dist <= long(pMismatches)){
                std::swap(r1, r2);
                p5primer = best_val.name;
                flipped = true;
                if (!keep) r1.setLCut(best_val.epos);
                pmatches++;
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
        if (best_val.dist <= long(pMismatches)){
            p3primer = best_val.name;
            if (!keep) r2.setLCut(best_val.epos);
            pmatches++;
        } else if (flip && p5primer=="None") {
            best_val.dist = pMismatches + 1;
            const std::string &seq2 = r1.get_seq();
            for ( auto it = primer3p.begin(); it != primer3p.end(); ++it ){
                const std::string p3Primer = it->second;
                test_val = bounded_edit_distance(p3Primer,  seq2,  pfloat,  pMismatches, pEndMismatches);
                if (test_val.dist < best_val.dist){
                    best_val = test_val;
                    best_val.name = it->first;
                }
                if (best_val.dist == 0) break;
            }
            if (best_val.dist <= long(pMismatches)){
                std::swap(r1, r2);
                p3primer = best_val.name;
                flipped = true;
                if (!keep) r2.setLCut(best_val.epos);
                pmatches++;
            }
        }
        counter.primer_match_counter(p5primer,p3primer);
        if (pmatches < mpmatches) {
            r1.setRCut(0);
            r2.setRCut(0);
        } else {
            if (flipped){
              counter.increment_flipped();
              r1.add_comment("Pf:Z:FLIP");
            }
            r1.add_comment("P5:Z:" + p5primer);
            r2.add_comment("P3:Z:" + p3primer);
        }
    }

    void check_read_se(SingleEndRead &se, PrimerCounters &counter, SeqMap &primer5p, SeqMap &primer3p, const size_t pMismatches, const size_t pEndMismatches, const size_t pfloat, const size_t flip, const size_t keep, const size_t mpmatches) {

        ALIGNPOS test_val, best_val;
        std::string p5primer = "None", p3primer = "None";
        size_t pmatches = 0;
        bool flipped = false;

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
        if (best_val.dist <= long(pMismatches)){
            p5primer = best_val.name;
            if (!keep) r1.setLCut(best_val.epos);
            pmatches++;
        } else if (flip) {
            best_val.dist = pMismatches + 1;
            const std::string &seq1 = r1.get_seq_rc();
            for ( auto it = primer5p.begin(); it != primer5p.end(); ++it ){
                const std::string p5Primer = it->second;
                test_val = bounded_edit_distance(p5Primer,  seq1,  pfloat,  pMismatches, pEndMismatches);
                if (test_val.dist < best_val.dist){
                    best_val = test_val;
                    best_val.name = it->first;
                }
                if (best_val.dist == 0) break;
            }
            if (best_val.dist <= long(pMismatches)){
                r1.set_read_rc();
                p5primer = best_val.name;
                flipped = true;
                if (!keep) r1.setLCut(best_val.epos);
                pmatches++;
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
        if (best_val.dist <= long(pMismatches)){
            p3primer = best_val.name;
            if (!keep) r1.setRCut(r1.getLength() -  best_val.epos);
            pmatches++;
        } else if (flip && p5primer=="None") {
            best_val.dist = pMismatches + 1;
            const std::string &temp = r1.get_seq();
            for ( auto it = primer3p.begin(); it != primer3p.end(); ++it ){
                const std::string p3Primer = it->second;
                test_val = bounded_edit_distance(p3Primer,  temp,  pfloat,  pMismatches, pEndMismatches);
                if (test_val.dist < best_val.dist){
                    best_val = test_val;
                    best_val.name = it->first;
                }
                if (best_val.dist == 0) break;
            }
            if (best_val.dist <= long(pMismatches)){
                r1.set_read_rc();
                p3primer = best_val.name;
                flipped = true;
                if (!keep) r1.setRCut(r1.getLength() -  best_val.epos);
                pmatches++;
            }
        }
        counter.primer_match_counter(p5primer,p3primer);
        if (pmatches < mpmatches) {
            r1.setRCut(0);
        } else {
            if (flipped){
              counter.increment_flipped();
              r1.add_comment("Pf:Z:FLIP");
            }
            r1.add_comment("P5:Z:" + p5primer);
            r1.add_comment("P3:Z:" + p3primer);
        }
    }

/* This is the helper class for Primer
 * The idea is ...
 * */
    template <class T, class Impl>
    void do_app(InputReader<T, Impl> &reader, std::shared_ptr<OutputWriter> pe, std::shared_ptr<OutputWriter> se, PrimerCounters &counter, const po::variables_map &vm) {
        SeqMap primer5p;
        if (vm.count("primers_5p")) {
            primer5p = fasta2dict(vm["primers_5p"].as<std::string>(), "p5primer");
        }

        SeqMap primer3p;
        if (vm.count("primers_3p")) {
            primer3p = fasta2dict(vm["primers_3p"].as<std::string>(), "p3primer");
        }

        counter.set_seqmap(primer5p, primer3p);

        const size_t pMismatches = vm["primer_mismatches"].as<size_t>();
        const size_t min_mprimers = vm["min_primer_matches"].as<size_t>();
        const size_t pEndMismatches = vm["primer_end_mismatches"].as<size_t>();
        const size_t pfloat = vm["float"].as<size_t>();
        const bool flip = vm["flip"].as<bool>();
        const bool keep = vm["keep"].as<bool>();
        const bool stranded = vm["stranded"].as<bool>();
        const size_t min_length = vm["min-length"].as<size_t>();
        bool no_orphan = vm["no-orphans"].as<bool>();

        while(reader.has_next()) {
            auto i = reader.next();
            PairedEndRead* per = dynamic_cast<PairedEndRead*>(i.get());
            if (per) {
                counter.input(*per);
                check_read_pe(*per, counter, primer5p, primer3p, pMismatches, pEndMismatches, pfloat, flip, keep, min_mprimers);
                per->checkDiscarded(min_length);
                counter.output(*per, no_orphan);
                writer_helper(per, pe, se, stranded, no_orphan);
            } else {
                SingleEndRead* ser = dynamic_cast<SingleEndRead*>(i.get());
                if (ser) {
                    counter.input(*ser);
                    check_read_se(*ser, counter, primer5p, primer3p, pMismatches, pEndMismatches, pfloat, flip, keep, min_mprimers);
                    ser->checkDiscarded(min_length);
                    counter.output(*ser);
                    writer_helper(ser, pe, se);
                } else {
                    throw std::runtime_error("Unknown read type");
                }
            }
        }
    }
};
#endif
