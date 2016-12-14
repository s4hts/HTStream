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

typedef std::unordered_map<std::string, size_t> Counter;
typedef std::multimap<std::string, std::size_t> seqLookup;
typedef std::shared_ptr<SingleEndRead> spReadBase;
typedef std::shared_ptr<std::vector<unsigned long long int>> histVec;

const size_t kmer = 8;
const size_t kmerBits = kmer*2;
const size_t step = 2;

uint8_t getBin(char c) {
    if (c == 'A') 
        return 0;
    else if (c == 'T')
        return 3;
    else if (c == 'C')
        return 1;
    else if (c == 'G')
        return 2;
    return 5;
}

seqLookup readOneMap(const std::string &seq1, size_t step) {

    seqLookup baseReadMap;
    size_t seqLen = seq1.length();

    for (size_t bp = 0; bp < seqLen - kmer; bp+=kmer) {
        std::cout << "Map " << seq1.substr(bp, kmer);
        baseReadMap.insert(std::make_pair(seq1.substr(bp, kmer), bp));
    }
    baseReadMap.insert(std::make_pair(seq1.substr(seqLen - kmer, kmer), seqLen - kmer));

    return baseReadMap;
}

spReadBase check(Read &r1, Read &r2, const size_t &loc1, const size_t &loc2, const size_t &maxMis, const size_t &minOverlap, const bool &adapterTrimming) {
    size_t minLoc = std::min(loc1, loc2);
    size_t loc1_t = loc1 - minLoc;
    size_t loc2_t = loc2 - minLoc;
    size_t maxLoop = std::min(r1.getLength() - loc1_t, r2.getLength() - loc2_t);
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

    if ( loc1_t >= loc2_t ) {
        finalSeq += seq1.substr(0, loc1_t);
        finalQual += qual1.substr(0, loc1_t);
    } else if (adapterTrimming) {
        r1.setLCut(loc1_t);
    }

    char bp;
    char qual;
   
    for (size_t i = 0; i < maxLoop; ++i) {
        read1_bp = loc1_t + i;
        read2_bp = loc2_t + i;
        
        if (seq1[read1_bp] == seq2[read2_bp]) {
            bp = seq1[read1_bp];
            qual = std::min(qual1[read1_bp] + qual2[read2_bp] - 33, 40 + 33); //addition of qual (minus just one of the ascii values 
        } else {
            bp = qual1[read1_bp] > qual2[read2_bp] ? seq1[read1_bp] : seq2[read2_bp];
            qual = std::max(qual1[read1_bp] - qual2[read2_bp] + 33, 1 + 33);
            if (++misMatches > maxMis) {
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
    
    if (r1.getLength() - loc1_t < r2.getLength()) {
        finalSeq += seq2.substr(maxLoop, r2.getLength() - maxLoop);
        finalQual += seq2.substr(maxLoop, r2.getLength() - maxLoop);
    } else if (adapterTrimming) {
        r2.setRCut(maxLoop);
    }

    spReadBase overlap(new SingleEndRead(Read(finalSeq, finalQual, "ID")));
    return overlap;
}

spReadBase getOverlappedReads(Read &r1, Read &r2, const seqLookup &seq1Map, const size_t &maxMis, const size_t &minOver, const size_t &checkLengths, const bool &adapterTrimming) {
    std::string seq2 = r2.get_seq_rc();

    for (size_t bp = 0; bp < checkLengths; ++bp) {
        auto test = seq1Map.equal_range(seq2.substr(bp, kmer));
        for (auto it = test.first; it != test.second; ++it) {
            spReadBase overlapped = check(r1, r2, (*it).second, bp, maxMis, minOver, adapterTrimming);
            if (overlapped != nullptr) {
                return overlapped;
            }
        }
    } 

    for (size_t bp = seq2.length() - (checkLengths + kmer); bp < seq2.length() - kmer ; ++bp) {
        auto test = seq1Map.equal_range(seq2.substr(bp, kmer));
        for (auto it = test.first; it != test.second; ++it) {
            spReadBase overlapped = check(r1, r2, (*it).second, bp, maxMis, minOver, adapterTrimming);
            if (overlapped != nullptr) {
                return overlapped;
            }
        }
    }
    return nullptr;

}

spReadBase check_read(PairedEndRead &pe , const size_t &maxMis, const size_t &minOver, histVec &insertLength, const bool &stranded, const size_t &checkLengths, const bool &adapterTrimming) {
    
    Read &r1 = pe.non_const_read_one();
    Read &r2 = pe.non_const_read_two();

    bool swapped = false;

    if (r1.getLength() < r2.getLength()) {
        std::swap(r1, r2);
        swapped = true;
    }

    seqLookup mOne = readOneMap(r1.get_seq(), step);
    spReadBase overlapped = getOverlappedReads(r1, r2, mOne, maxMis, minOver, checkLengths, adapterTrimming) ;

    if (insertLength && overlapped) { //overlap plus writing out
        std::string s = overlapped->get_read().get_seq();
        size_t len = s.length();
        if (insertLength->size() < len) {
            insertLength->resize(len + 1);
        }
        ++((*insertLength)[len]);

    } else if (insertLength) {
        if (insertLength->size() < 1) {
            insertLength->resize(1);
        }
        ++(*insertLength)[0];
    }
    return overlapped;
}


template <class T, class Impl>
void helper_overlapper(InputReader<T, Impl> &reader, std::shared_ptr<OutputWriter> pe, std::shared_ptr<OutputWriter> se, Counter& counters, const size_t &maxMis, const size_t &minOver, histVec &insertLength, const bool &stranded, const size_t &min_length, const size_t &checkLengths, const bool &adapterTrimming) {
    
    while(reader.has_next()) {
        auto i = reader.next();

        ++counters["TotalRecords"];
        PairedEndRead* per = dynamic_cast<PairedEndRead*>(i.get());        
        if (per) {
            spReadBase overlapped = check_read(*per, maxMis, minOver, insertLength, stranded, checkLengths, adapterTrimming);
            if (!overlapped || adapterTrimming) {
                per->checkDiscarded(min_length);    
                writer_helper(per, pe, se, stranded, counters); 
            } else {
                overlapped->checkDiscarded(min_length);    
                writer_helper(overlapped.get(), pe, se, stranded, counters);
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

