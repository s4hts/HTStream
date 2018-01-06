/*The idea of this program is to create a quick kmer lookup table for phix.
 * It accomplishes this by making 2 arrays of 1<<16 size, because we represent
 * each possible kmer in the 2 bit format (A -> 00, T -> 11, C -> 01, G->10) and
 * need foward and Reverse strands. Each read is then parsed into 8mers (again
 * represented by their 2 bit format, and checked against these lookup tables
 * (simply the arrays)*/ 

#ifndef AT_TRIM_H
#define AT_TRIM_H
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

extern template class InputReader<PairedEndRead, InterReadImpl>;

typedef std::unordered_multimap<std::string, std::size_t> seqLookup;
typedef std::shared_ptr<SingleEndRead> spReadBase;

/*Create the quick lookup table
 * Multi map because a single kemr could appear multiple places*/
seqLookup readOneMap(std::string seq1, const size_t kmer, const size_t kmerOffset) {

    seqLookup baseReadMap;
    std::string::iterator it;
    for ( it = seq1.begin() ; it < seq1.end() - ( static_cast<long> ( kmerOffset + kmer ) ) ; it += static_cast<long> ( kmerOffset )   ) {
        baseReadMap.insert(std::make_pair( std::string ( it, it+ static_cast<long> ( kmer )  ) , it - seq1.begin() ));
    }

    if ( seq1.begin() + static_cast<long> (kmer) > seq1.end() ) {
        it = seq1.end() - static_cast<long> ( kmer );
        baseReadMap.insert(std::make_pair(  std::string( it , it + static_cast<long> (kmer)), it - seq1.begin()   ) );
    }

    return baseReadMap;
}
/*If adapater trimming is turned on that means adapter trimming and do not overlap
 * so trim adapter, but don't worry about the overlap.
 * however we still need to change the overlap
 * 
 * Within the overlap if they are the same bp, then add q scores
 * If they are different bp, subtract q scores and take the larger quality bp*/ 
spReadBase checkIfOverlap(Read &r1, Read &r2, size_t loc1, size_t loc2, const double misDensity, size_t minOverlap) {
    size_t minLoc = std::min(loc1, loc2);
    int  loc1_t = loc1 - minLoc;
    int loc2_t = loc2 - minLoc;
    int r1_len = r1.getLength();
    int r2_len = r2.getLength();

    size_t maxLoop = std::min(r1_len - loc1_t, r2_len - loc2_t);
    size_t maxMis = static_cast<size_t>(maxLoop * misDensity);
    
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
spReadBase getOverlappedReads(Read &r1, Read &r2, const seqLookup &seq1Map,  const double misDensity, const size_t &minOver, const size_t &checkLengths, const size_t kmer) {
    std::string seq2 = r2.get_seq_rc();
    for (size_t bp = 0; bp < checkLengths; ++bp) {
        /*Do a quick check if the shorter read kmer shows up in longer read (read 2)
         * If it does, then try the brute force approach*/
        auto test = seq1Map.equal_range(seq2.substr(bp, kmer));
        for (auto it = test.first; it != test.second; ++it) {
            spReadBase overlapped = checkIfOverlap(r1, r2, it->second, bp, misDensity, minOver);
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
            spReadBase overlapped = checkIfOverlap(r1, r2, it->second, bp, misDensity, minOver);
            if (overlapped != nullptr) {
                return overlapped;
            }
        }
    }
    return nullptr;
}

spReadBase check_read(PairedEndRead &pe, const double misDensity, const size_t &minOver, const size_t &checkLengths, const size_t kmer, const size_t kmerOffset) {
    
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
    spReadBase overlapped = getOverlappedReads(r1, r2, mOne, misDensity, minOver, checkL, kkmer) ;
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
 * They can either overlap just barely on the ends we call - these lins (long insert)
 * They can also overlap way to much to the point they have adapters
 * in the read - or a sin (short insert)
 * With a lin it is useful to have a higher confidence in the bases in the overlap and longer read
 * With a sin it is useful to have the higher confidence as well as removing the adapters*/
template <class T, class Impl>
void helper_overlapper(InputReader<T, Impl> &reader, std::shared_ptr<OutputWriter> pe, std::shared_ptr<OutputWriter> se, OverlappingCounters &counters, const double misDensity, const size_t minOver, const bool stranded, const size_t min_length, const size_t checkLengths, const size_t kmer, const size_t kmerOffset, bool no_orphan = false ) {
    while(reader.has_next()) {
        auto i = reader.next();
        PairedEndRead* per = dynamic_cast<PairedEndRead*>(i.get());        
        if (per) {
            counters.input(*per);
            spReadBase overlapped = check_read(*per, misDensity, minOver, checkLengths, kmer, kmerOffset);
            if (!overlapped) {
                writer_helper(per, pe, se, stranded, no_orphan); //write out as is
                counters.output(*per);
            } else if (overlapped) { //if there is an overlap
                overlapped->checkDiscarded(min_length);
                unsigned int origLength = unsigned(per->non_const_read_one().getLength());
                writer_helper(overlapped.get(), pe, se, stranded);
                counters.output(*overlapped, origLength);
            }
        } else {
            SingleEndRead* ser = dynamic_cast<SingleEndRead*>(i.get());
            
            if (ser) {
                counters.input(*ser);
                ser->checkDiscarded(min_length);
                counters.output(*ser);
                writer_helper(ser, pe, se, stranded);
            } else {
                throw std::runtime_error("Unknown read type");
            }
        }
    }
}

#endif

