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

typedef std::unordered_multimap<std::string, std::size_t> seqLookup;
typedef std::shared_ptr<SingleEndRead> spReadBase;
typedef std::shared_ptr<std::vector<unsigned long long int>> histVec;

/*Create the quick lookup table
 * Multi map because a single kemr could appear multiple places*/
seqLookup readOneMap(const std::string &seq1, const size_t kmer, const size_t kmerOffset) {

    seqLookup baseReadMap;
    size_t seqLen = seq1.length() - 1;

    for (size_t bp = 0; bp < seqLen - kmer; bp+=kmerOffset) {
        baseReadMap.insert(std::make_pair(seq1.substr(bp, kmer), bp));
    }
    baseReadMap.insert(std::make_pair(seq1.substr(seqLen - kmer, kmer), seqLen - kmer));

    return baseReadMap;
}
/*If adapater trimming is turned on that means adapter trimming and do not overlap
 * so trim adapter, but don't worry about the overlap.
 * however we still need to change the overlap
 * 
 * Within the overlap if they are the same bp, then add q scores
 * If they are different bp, subtract q scores and take the larger quality bp*/ 
spReadBase checkIfOverlap(Read &r1, Read &r2, size_t loc1, size_t loc2, const double misDensity, size_t minOverlap, bool adapterTrimming) {
    size_t minLoc = std::min(loc1, loc2);
    size_t loc1_t = loc1 - minLoc;
    size_t loc2_t = loc2 - minLoc;
    size_t maxLoop = std::min(r1.getLength() - loc1_t, r2.getLength() - loc2_t);
    size_t maxMis = static_cast<size_t>(maxLoop * misDensity);

    if (maxLoop <= minOverlap) {
        return nullptr;
    }

    size_t read1_bp;
    size_t read2_bp;

    const std::string &seq1 = r1.get_seq();
    const std::string &seq2 = r2.get_seq_rc();
    const std::string &qual1 = r1.get_qual();
    const std::string &qual2 = r2.get_qual_rc();

    std::string finalSeq;
    std::string finalQual;
    size_t misMatches = 0;

    char bp;
    char qual;

    for (size_t i = 0; i < maxLoop; ++i) {
        read1_bp = loc1_t + i;
        read2_bp = loc2_t + i;
        
        if (seq1[read1_bp] == seq2[read2_bp]) {
            bp = seq1[read1_bp];
            qual = static_cast<char>(std::min(qual1[read1_bp] + qual2[read2_bp] - 33, 40 + 33)); //addition of qual (minus just one of the ascii values 
        } else {
            bp = qual1[read1_bp] > qual2[read2_bp] ? seq1[read1_bp] : seq2[read2_bp];
            qual = static_cast<char>(std::max(qual1[read1_bp] - qual2[read2_bp] + 33, 1 + 33));
            
            if (++misMatches > maxMis) {
                /*Not valid match*/
                return nullptr;
            }
        }
        if (!adapterTrimming) {
            finalSeq += bp;
            finalQual += qual;
        } else {
            r1.changeSeq(read1_bp, bp);
            r1.changeQual(read1_bp, qual);
            r2.changeSeq(read2_bp, bp);
            r2.changeQual(read2_bp, qual);
        }
    }

    if ( loc1_t >= loc2_t ) {
        finalSeq = seq1.substr(0, loc1_t) + finalSeq;
        finalQual = qual1.substr(0, loc1_t) + finalQual;
    } else if (adapterTrimming) {
        r1.setLCut(loc1_t);
    }
    
    if (r1.getLength() - loc1_t < r2.getLength()) {
        finalSeq += seq2.substr(maxLoop, r2.getLength() - maxLoop);
        finalQual += qual2.substr(maxLoop, r2.getLength() - maxLoop);
    } else if (adapterTrimming) {
        r2.setRCut(maxLoop);
    }

    spReadBase overlap(new SingleEndRead(Read(finalSeq, finalQual, r1.get_id())));
    return overlap;
}

/*Because of the way overlapping works, you only need to check the ends of the shorter read*/
spReadBase getOverlappedReads(Read &r1, Read &r2, const seqLookup &seq1Map,  const double misDensity, const size_t &minOver, const size_t &checkLengths, const bool &adapterTrimming, const size_t kmer) {
    std::string seq2 = r2.get_seq_rc();
    for (size_t bp = 0; bp < checkLengths; ++bp) {
        /*Do a quick check if the shorter read kmer shows up in longer read (read 2)
         * If it does, then try the brute force approach*/
        auto test = seq1Map.equal_range(seq2.substr(bp, kmer));
        for (auto it = test.first; it != test.second; ++it) {
            spReadBase overlapped = checkIfOverlap(r1, r2, it->second, bp, misDensity, minOver, adapterTrimming);
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
            spReadBase overlapped = checkIfOverlap(r1, r2, it->second, bp, misDensity, minOver, adapterTrimming);
            if (overlapped != nullptr) {
                return overlapped;
            }
        }
    }
    return nullptr;

}

spReadBase check_read(PairedEndRead &pe , Counter &counters, const double misDensity, const size_t &minOver, histVec &insertLength, const size_t &checkLengths, const bool &adapterTrimming, const size_t kmer, const size_t kmerOffset, const size_t minLength) {
    
    Read &r1 = pe.non_const_read_one();
    Read &r2 = pe.non_const_read_two();

    bool swapped = false;
    /*Read1 is always longer than Read 2)*/
    if (r1.getLength() < r2.getLength()) {
        std::swap(r1, r2);
        swapped = true;
    }
    /*Create a map with non-overlapping kmers*/
    seqLookup mOne = readOneMap(r1.get_seq(), kmer, kmerOffset);
    /*returns null if no much
     * r1 and r2 and passed by ref in case only adapter trimming is on*/
    spReadBase overlapped = getOverlappedReads(r1, r2, mOne, misDensity, minOver, checkLengths, adapterTrimming, kmer) ;
    //we need to check it overlapper is greater than min length;

    if (overlapped && overlapped->non_const_read_one().getLength() < minLength) {
        if (insertLength) { /*No overlap*/
            ++(*insertLength)[0];
        }
        return nullptr;
    }

    if (overlapped) {
        std::string s = overlapped->get_read().get_seq();
        size_t len = s.length();

        if (insertLength) { //overlap plus writing out
            /*This is important for the hist file
            * Shows distribution of lins and sins*/
            if (insertLength->size() < len) {
                insertLength->resize(len + 1);
            }
            ++((*insertLength)[len]);
        }
        if (len > r1.getLength() && len > r2.getLength()) {
            ++counters["Lins"];
        } else {
            ++counters["Sins"];
        }
    } else if (insertLength) {
        /*No overlap*/
        ++(*insertLength)[0];
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
void helper_overlapper(InputReader<T, Impl> &reader, std::shared_ptr<OutputWriter> pe, std::shared_ptr<OutputWriter> se, Counter& counters, const double misDensity, const size_t &minOver, histVec &insertLength, const bool &stranded, const size_t &min_length, const size_t &checkLengths, const bool &adapterTrimming, const size_t kmer, const size_t kmerOffset, bool no_orphan = false ) {
    
    while(reader.has_next()) {
        auto i = reader.next();
        //Saves check
        if (insertLength) {
            if (insertLength->size() < 1) {
                insertLength->resize(1);
            }
        }

        ++counters["TotalRecords"];
        PairedEndRead* per = dynamic_cast<PairedEndRead*>(i.get());        
        if (per) {
            spReadBase overlapped = check_read(*per, counters, misDensity, minOver, insertLength, checkLengths, adapterTrimming, kmer, kmerOffset, min_length);
            if (!overlapped) {
                writer_helper(per, pe, se, stranded, counters); //write out as is
            } else if (overlapped) { //if there is an overlap
                if (adapterTrimming) { //just adapter trimming, removes adapters
                    per->checkDiscarded(min_length);    
                    writer_helper(per, pe, se, stranded, counters, no_orphan); 
                } else {
                    overlapped->checkDiscarded(min_length);    
                    writer_helper(overlapped.get(), pe, se, stranded, counters);
                }
            }


        } else {
            SingleEndRead* ser = dynamic_cast<SingleEndRead*>(i.get());
            
            if (ser) {
                ser->checkDiscarded(min_length);
                writer_helper(ser, pe, se, stranded, counters);
            } else {
                throw std::runtime_error("Unknown read type");
            }
        }
    }
}

#endif

