#ifndef OVERLAPPER_H
#define OVERLAPPER_H
//  this is so we can implment hash function for dynamic_bitset
#define BOOST_DYNAMIC_BITSET_DONT_USE_FRIENDS
#define STARTS 4098

#include "ioHandler.h"
#include <map>
#include <unordered_map>
#include <boost/dynamic_bitset.hpp>
#include <boost/functional/hash.hpp>
#include <algorithm>
#include <bitset>
#include <utility>
#include "utils.h"
#include "main_template.h"

class OverlappingCounters : public Counters {
public:
    std::vector<uint_fast64_t> insertLength;

    uint64_t sins = 0;
    uint64_t mins = 0;
    uint64_t lins = 0;
    uint64_t Adapter_Trim = 0;
    uint64_t Adapter_BpTrim = 0;

    uint64_t SE_Discarded = 0;

    uint64_t R1_Discarded = 0;
    uint64_t R2_Discarded = 0;
    uint64_t PE_Discarded = 0;

    OverlappingCounters(const std::string &program_name, const po::variables_map &vm) : Counters::Counters(program_name, vm) {

        insertLength.resize(1);

        fragment.push_back(std::forward_as_tuple("short_inserts", sins));
        fragment.push_back(std::forward_as_tuple("medium_inserts", mins));
        fragment.push_back(std::forward_as_tuple("long_inserts", lins));
        fragment.push_back(std::forward_as_tuple("adapterTrim", Adapter_Trim));
        fragment.push_back(std::forward_as_tuple("adapterBpTrim", Adapter_BpTrim));

        se.push_back(std::forward_as_tuple("SE_discarded", SE_Discarded));

        r1.push_back(std::forward_as_tuple("R1_discarded", R1_Discarded));
        r2.push_back(std::forward_as_tuple("R2_discarded", R2_Discarded));
        pe.push_back(std::forward_as_tuple("PE_discarded", PE_Discarded));
    }

    using Counters::output;

    void output(SingleEndRead &ser, uint_fast64_t origLength) {
        Read &one = ser.non_const_read_one();
        if (!one.getDiscard() && !origLength) { // original SE read, passed through
            ++SE_Out;
            ++TotalFragmentsOutput;
        } else if (!one.getDiscard() && origLength) { // overlapped read
            if (one.getLength() < origLength) {
                ++sins; //adapters must be had (short insert)
                ++Adapter_Trim;
                Adapter_BpTrim += (origLength - one.getLength());
            } else {
                ++mins; //must be a long insert
            }
            if ( one.getLength() + 1 > insertLength.size() ) {
                insertLength.resize(one.getLength() + 1);
            }
            ++insertLength[one.getLength()];
            ++SE_Out;
            ++TotalFragmentsOutput;
        } else {
            ++SE_Discarded;
        }
    }

    virtual void output(PairedEndRead &per, bool no_orphans = false)  {
        Read &one = per.non_const_read_one();
        Read &two = per.non_const_read_two();
        if (!one.getDiscard() && !two.getDiscard()) {
            if ((one.getLengthTrue() < one.getLength()) || (two.getLengthTrue() < two.getLength())) {
                // this should never be the case
                std::cerr << "Should never see this 1\n";
            }
            ++lins;
            ++PE_Out;
            ++TotalFragmentsOutput;
        } else if (!one.getDiscard() && !no_orphans) { //if stranded RC
            if (one.getLengthTrue() < one.getLength()) {
                // this should never be the case
                std::cerr << "Should never see this 2\n";
            }
            ++SE_Out;
            ++R2_Discarded;
            ++TotalFragmentsOutput;
        } else if (!two.getDiscard() && !no_orphans) { // Will never be RC
            if (two.getLengthTrue() < two.getLength()) {
                // this should never be the case
                std::cerr << "Should never see this 3\n";
            }
            ++SE_Out;
            ++R1_Discarded;
            ++TotalFragmentsOutput;
        } else {
            ++PE_Discarded;
        }
    }

    virtual void write_out() {

        std::vector<Vector> iLength;
        for (size_t i = 1; i < insertLength.size(); ++i) {
            if (insertLength[i] > 0) {
                iLength.push_back(std::forward_as_tuple(i, insertLength[i]));
            }
        }

        initialize_json();

        start_sublabel("Program_details");
        write_values(pd, 2);
        end_sublabel();

        start_sublabel("Fragment");
        write_values(fragment, 2);
        write_vector("readlength_histogram",iLength, 2);
        end_sublabel();

        start_sublabel("Single_end");
        write_values(se, 2);
        end_sublabel();

        start_sublabel("Paired_end");
        write_values(pe, 2);
        write_values(r1, 2);
        write_values(r2, 2);
        end_sublabel();

        finalize_json();
    }
};

class Overlapper: public MainTemplate<OverlappingCounters, Overlapper> {
public:

    void add_extra_options(po::options_description &desc) {
        setDefaultParamsCutting(desc);
        // no-orphans|n ; stranded|s ; min-length|m
        setDefaultParamsOverlapping(desc);
    }

    Overlapper() {
        program_name = "hts_Overlapper";
        app_description =
            "The hts_Overlapper application attempts to overlap paired end reads\n";
        app_description += "  to produce the original transcript, trim adapters, and in some\n";
        app_description += "  cases, correct sequencing errors.\n";
        app_description += "Reads come in three flavors:\n";
        app_description += "  sins: Reads produced from an insert shorter than the read length\n";
        app_description += "        will result in a single read in the orientation of R1, and have the\n";
        app_description += "        adapter bases trimmed to produce a SE read.\n";
        app_description += "  mins: Reads produced from a medium-insert greater than read length, but\n";
        app_description += "        somewhat shorter than 2x read length will produce a SE read in the\n";
        app_description += "        orientation of R1.\n";
        app_description += "  lins: Reads produced from long-inserts which do not overlap\n";
        app_description += "        significantly, resulting in a PE read.\n";
    }

/* Within the overlap if they are the same bp, then add q scores
 * If they are different bp, subtract q scores and take the larger quality bp*/
    SingleEndReadPtr checkIfOverlap(Read &r1, Read &r2, size_t loc1, size_t loc2, const double misDensity, const size_t &mismatch, const size_t minOverlap) {
        size_t minLoc = std::min(loc1, loc2);
        size_t loc1_t = loc1 - minLoc;
        size_t loc2_t = loc2 - minLoc;
        size_t r1_len = r1.getLength();
        size_t r2_len = r2.getLength();

        size_t maxLoop = std::min(r1_len - loc1_t, r2_len - loc2_t);
        size_t maxMis = std::min(mismatch, static_cast<size_t>(maxLoop * misDensity));

        const std::string &seq1 = r1.get_seq();
        const std::string &seq2 = r2.get_seq_rc();

        auto i1 = seq1.begin();
        std::advance(i1, loc1_t);
        auto i2 = seq2.begin();
        std::advance(i2, loc2_t);
        if (maxLoop < minOverlap || !threshold_mismatches( i1, i2 , maxLoop, maxMis ) ) {
            return nullptr;
        }

        size_t read1_bp;
        size_t read2_bp;

        const std::string &qual1 = r1.get_qual();
        const std::string &qual2 = r2.get_qual_rc();

        std::string finalSeq;
        std::string finalQual;

        char bp;
        char qual;

        for (size_t i = 0; i < maxLoop; ++i) {
            read1_bp = loc1_t + i;
            read2_bp = loc2_t + i;

            if (seq1[read1_bp] == seq2[read2_bp]) {
                bp = seq1[read1_bp];
                qual = static_cast<char>(std::min(qual1[read1_bp] + qual2[read2_bp] - 33, 40 + 33)); //addition of qual (minus just one of the ascii values
            } else {
                bp = qual1[read1_bp] >= qual2[read2_bp] ? seq1[read1_bp] : seq2[read2_bp];
                qual = static_cast<char>(std::max(qual1[read1_bp] - qual2[read2_bp] + 33, 1 + 33));
            }
            finalSeq += bp;
            finalQual += qual;
        }

        /*R2 is always shorter or equal to R1*/

        if ( loc1_t > loc2_t ) { //We are going to snag the start if it is a lin (notice, if they are both zero then we dont' want to snag anything
            finalSeq = seq1.substr(0, loc1_t) + finalSeq;
            finalQual = qual1.substr(0, loc1_t) + finalQual;
        }

        if (r1_len - loc1_t < r2_len - loc2_t) {
            finalSeq += seq2.substr(maxLoop, r2_len - maxLoop);
            finalQual += qual2.substr(maxLoop, r2_len - maxLoop);
        }

        SingleEndReadPtr overlap(new SingleEndRead(Read(finalSeq, finalQual, r1.get_id())));
        return overlap;
    }

/*Because of the way overlapping works, you only need to check the ends of the shorter read*/
    SingleEndReadPtr getOverlappedReads(Read &r1, Read &r2, const seqLookup &seq1Map,  const double misDensity, const size_t mismatch, const size_t &minOver, const size_t &checkLengths, const size_t kmer) {
        std::string seq2 = r2.get_seq_rc();
        for (size_t bp = 0; bp < (checkLengths - kmer); ++bp) {
            /*Do a quick check if the shorter read kmer shows up in longer read (read 2)
             * If it does, then try the brute force approach*/
            auto test = seq1Map.equal_range(seq2.substr(bp, kmer));
            for (auto it = test.first; it != test.second; ++it) {
                SingleEndReadPtr overlapped = checkIfOverlap(r1, r2, it->second, bp, misDensity, mismatch, minOver);
                if (overlapped != nullptr) {
                    return overlapped;
                }
            }
        }
        for (size_t bp = seq2.length() - (checkLengths + kmer); bp <= (seq2.length() - kmer) ; ++bp) {
            /*Do a quick check if the shorter read kmer shows up in longer read (read 2)
             * If it does, then try the brute force approach*/
            auto test = seq1Map.equal_range(seq2.substr(bp, kmer));
            for (auto it = test.first; it != test.second; ++it) {
                SingleEndReadPtr overlapped = checkIfOverlap(r1, r2, it->second, bp, misDensity, mismatch, minOver);
                if (overlapped != nullptr) {
                    return overlapped;
                }
            }
        }
        return nullptr;
    }

    SingleEndReadPtr check_read(PairedEndRead &pe, const double misDensity, const size_t &mismatch, const size_t &minOver, const size_t &checkLengths, const size_t kmer, const size_t kmerOffset) {

        Read &r1 = pe.non_const_read_one();
        Read &r2 = pe.non_const_read_two();

        bool swapped = false;
        /* Read1 needs to always be longer than Read 2 */
        if (r1.getLength() < r2.getLength()) {
            std::swap(r1, r2);
            swapped = true;
        }
        /* checkL needs to be as long as or longer than the shortest read */
        size_t checkL = std::min(r2.getLength(), checkLengths);
        /* checkL needs to be as long as or longer than the shortest read */
        size_t kkmer = std::min(r2.getLength(), kmer);
        /* Create a map with non-overlapping kmers */
        seqLookup mOne = readOneMap(r1.get_seq(), kkmer, kmerOffset);
        /* returns null if no much
         * r1 and r2 and passed by ref in case only adapter trimming is on */
        SingleEndReadPtr overlapped = getOverlappedReads(r1, r2, mOne, misDensity, mismatch, minOver, checkL, kkmer) ;
        if (!overlapped && swapped){
            std::swap(r1, r2);
        } else if (swapped){
            overlapped->set_read_rc();
        }

        return overlapped;
    }

/*This is the helper class for overlap
 * The idea is in the wet lab, they set up sequences to sequences toward each other
 * Sometimes, these squences can overlap
 * There are two cases in which they can overlap
 * They can either overlap just barely on the ends we call - these mins (medium insert)
 * They can also overlap way to much to the point they have adapters
 * in the read - or a sin (short insert) [adapter are overhangs in an overlap]
 * With a min it is useful to have a higher confidence in the bases in the overlap and longer read
 * With a sin it is useful to have the higher confidence as well as removing the adapters*/
    template <class T, class Impl>
    void do_app(InputReader<T, Impl> &reader, std::shared_ptr<OutputWriter> pe, std::shared_ptr<OutputWriter> se, OverlappingCounters &counters, const po::variables_map &vm) {
        const double misDensity = vm["max-mismatch-errorDensity"].as<double>();
        const size_t mismatch = vm["max-mismatch"].as<size_t>();
        const size_t minOver = vm["min-overlap"].as<size_t>();
        const bool stranded = vm["stranded"].as<bool>();
        const size_t min_length = vm["min-length"].as<size_t>();
        const size_t checkLengths = vm["check-lengths"].as<size_t>();
        const size_t kmer = vm["kmer"].as<size_t>();
        const size_t kmerOffset = vm["kmer-offset"].as<size_t>();
        bool no_orphan = vm["no-orphans"].as<bool>();

        while(reader.has_next()) {
            auto i = reader.next();
            PairedEndRead* per = dynamic_cast<PairedEndRead*>(i.get());
            if (per) {
                counters.input(*per);
                SingleEndReadPtr overlapped = check_read(*per, misDensity, mismatch, minOver, checkLengths, kmer, kmerOffset);
                if (!overlapped) {
                    counters.output(*per, no_orphan);
                    writer_helper(per, pe, se, stranded, no_orphan); //write out as is
                } else if (overlapped) { //if there is an overlap
                    overlapped->checkDiscarded(min_length);
                    unsigned int origLength = std::max(unsigned(per->non_const_read_one().getLength()),unsigned(per->non_const_read_two().getLength()));
                    counters.output(*overlapped, origLength);
                    writer_helper(overlapped.get(), pe, se, stranded);
                }
            } else {
                SingleEndRead* ser = dynamic_cast<SingleEndRead*>(i.get());

                if (ser) {
                    counters.input(*ser);
                    ser->checkDiscarded(min_length);
                    counters.output(*ser, 0);
                    writer_helper(ser, pe, se, stranded);
                } else {
                    throw std::runtime_error("Unknown read type");
                }
            }
        }
    }
};
#endif
