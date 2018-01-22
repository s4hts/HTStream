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

extern template class InputReader<SingleEndRead, SingleEndReadFastqImpl>;
extern template class InputReader<PairedEndRead, PairedEndReadFastqImpl>;
extern template class InputReader<PairedEndRead, InterReadImpl>;
extern template class InputReader<ReadBase, TabReadImpl>;

typedef std::shared_ptr<SingleEndRead> spReadBase;

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

    OverlappingCounters(const std::string &statsFile, bool appendStats, const std::string &program_name, const std::string &notes) : Counters::Counters(statsFile, appendStats, program_name, notes) {

        insertLength.resize(1);

        generic.push_back(std::forward_as_tuple("sins", sins));
        generic.push_back(std::forward_as_tuple("mins", mins));
        generic.push_back(std::forward_as_tuple("lins", lins));
        generic.push_back(std::forward_as_tuple("adapterTrim", Adapter_Trim));
        generic.push_back(std::forward_as_tuple("adapterBpTrim", Adapter_BpTrim));

        se.push_back(std::forward_as_tuple("SE_discarded", SE_Discarded));

        pe.push_back(std::forward_as_tuple("R1_discarded", R1_Discarded));
        pe.push_back(std::forward_as_tuple("R2_discarded", R2_Discarded));
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

        write_labels(generic);
        write_vector("readlength_histogram",iLength);
        write_sublabels("Single_end", se);
        write_sublabels("Paired_end", pe);

        finalize_json();        
    }
};

/* Within the overlap if they are the same bp, then add q scores
 * If they are different bp, subtract q scores and take the larger quality bp*/ 
spReadBase checkIfOverlap(Read &r1, Read &r2, size_t loc1, size_t loc2, const double misDensity, const size_t &mismatch, const size_t minOverlap) {
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
    if (maxLoop <= minOverlap || !threshold_mismatches( i1, i2 , maxLoop, maxMis ) ) {
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

    spReadBase overlap(new SingleEndRead(Read(finalSeq, finalQual, r1.get_id())));
    return overlap;
}

/*Because of the way overlapping works, you only need to check the ends of the shorter read*/
spReadBase getOverlappedReads(Read &r1, Read &r2, const seqLookup &seq1Map,  const double misDensity, const size_t mismatch, const size_t &minOver, const size_t &checkLengths, const size_t kmer) {
    std::string seq2 = r2.get_seq_rc();
    for (size_t bp = 0; bp < checkLengths; ++bp) {
        /*Do a quick check if the shorter read kmer shows up in longer read (read 2)
         * If it does, then try the brute force approach*/
        auto test = seq1Map.equal_range(seq2.substr(bp, kmer));
        for (auto it = test.first; it != test.second; ++it) {
            spReadBase overlapped = checkIfOverlap(r1, r2, it->second, bp, misDensity, mismatch, minOver);
            if (overlapped != nullptr) {
                return overlapped;
            }
        }
    } 
    for (size_t bp = seq2.length() - (checkLengths + kmer); bp < seq2.length() - kmer ; ++bp) {
        /*Do a quick check if the shorter read kmer shows up in longer read (read 2)
         * If it does, then try the brute force approach*/
        auto test = seq1Map.equal_range(seq2.substr(bp, kmer));
        for (auto it = test.first; it != test.second; ++it) {
            spReadBase overlapped = checkIfOverlap(r1, r2, it->second, bp, misDensity, mismatch, minOver);
            if (overlapped != nullptr) {
                return overlapped;
            }
        }
    }
    return nullptr;
}

spReadBase check_read(PairedEndRead &pe, const double misDensity, const size_t &mismatch, const size_t &minOver, const size_t &checkLengths, const size_t kmer, const size_t kmerOffset) {
    
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
    spReadBase overlapped = getOverlappedReads(r1, r2, mOne, misDensity, mismatch, minOver, checkL, kkmer) ;
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
void helper_overlapper(InputReader<T, Impl> &reader, std::shared_ptr<OutputWriter> pe, std::shared_ptr<OutputWriter> se, OverlappingCounters &counters, const double misDensity, const size_t mismatch, const size_t minOver, const bool stranded, const size_t min_length, const size_t checkLengths, const size_t kmer, const size_t kmerOffset, bool no_orphan = false ) {
    while(reader.has_next()) {
        auto i = reader.next();
        PairedEndRead* per = dynamic_cast<PairedEndRead*>(i.get());        
        if (per) {
            counters.input(*per);
            spReadBase overlapped = check_read(*per, misDensity, mismatch, minOver, checkLengths, kmer, kmerOffset);
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

#endif

