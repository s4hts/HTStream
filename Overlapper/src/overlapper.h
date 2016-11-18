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
   /* uint16_t forwardBits = 0;
    size_t kmerCheck = kmer - 1;
    uint8_t bin;*/


    seqLookup baseReadMap;
    std::cout << seq1 << '\n';
    size_t seqLen = seq1.length();
    /*8 is a random number*/
    for (size_t bp = 0; bp < seqLen - kmer; bp+=kmer) {
        std::cout << seq1.substr(bp, kmer) << '\n';
        baseReadMap.insert(std::make_pair(seq1.substr(bp, kmer), bp));
       /* bin = getBin(seq1[bp]);
        if (bin != 5) {
            forwardBits &= ~(3 << 14); //only check forward strand against both RC and non-rc
            forwardBits <<= 2;
            forwardBits ^= bin;
        } else {
            kmerCheck = kmer - 1;
        }

        if (!kmerCheck) {
        } else {
            --kmerCheck;    
        }*/
    }
    baseReadMap.insert(std::make_pair(seq1.substr(seqLen - kmer, kmer), seqLen - kmer));
    std::cout << seq1.substr(seqLen - kmer, kmer) << '\n';

    return baseReadMap;
}

spReadBase check(Read &r1, Read &r2, const size_t &loc1, const size_t &loc2, const size_t &maxMis, const size_t &minOverlap) {
    size_t minLoc = std::min(loc1, loc2);
    size_t loc1_t = loc1 - minLoc;
    size_t loc2_t = loc2 - minLoc;
    size_t maxLoop = std::min(r1.getLength() - loc1_t, r2.getLength() - loc2_t);
    
    if (maxLoop <= minOverlap) {
        std::cout << "NULL\n";
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
    std::cout << "Loc1 " << loc1_t << " Loc2 " << loc2_t << '\n';

    if ( loc1_t >= loc2_t) {
        std::cout << "SNAG FIRST HALF\n";
        finalSeq += seq1.substr(0, loc1_t);
        finalQual += qual1.substr(0, loc1_t);
    }

    for (size_t i = 0; i < maxLoop; ++i) {
        read1_bp = loc1_t + i;
        read2_bp = loc2_t + i;
        if (seq1[read1_bp] == seq2[read2_bp]) {
            finalSeq += seq1[read1_bp];
            finalQual += std::min(qual1[read1_bp] + qual2[read2_bp] - 33, 40 + 33); //addition of qual (minus just one of the ascii values 
        } else {
            finalSeq += qual1[read1_bp] > qual2[read2_bp] ? seq1[read1_bp] : seq2[read2_bp];
            finalQual += std::max(qual1[read1_bp] - qual2[read2_bp] + 33, 1 + 33);
            ++misMatches;
        }

    }
    
    if (r1.getLength() - loc1_t < r2.getLength()) {
        std::cout << "Snag second half\n";
        finalSeq += seq2.substr(maxLoop, r2.getLength() - maxLoop);
        finalQual += seq2.substr(maxLoop, r2.getLength() - maxLoop);
    }
    std::cout << "Final\n";
    std::cout << finalSeq << '\n';
    std::cout << finalQual << '\n';
    spReadBase overlap(new SingleEndRead(Read(finalSeq, finalQual, "ID")));
    const Read &r = overlap->get_read();
    std::cout << r.get_seq() << '\n';
    return overlap;
}

spReadBase getOverlappedReads(Read &r1, Read &r2, const seqLookup &seq1Map) {
    uint16_t forwardBits = 0;
    size_t kmerCheck = kmer - 1;
    uint8_t bin;
    std::vector< std::pair<size_t, size_t> > loc;
    std::cout << "Start\n";
    std::cout << r2.get_seq_rc() << '\n';
    for (size_t bp = 0; bp < r2.getLength() - kmer; ++bp) {
        auto test = seq1Map.equal_range(r2.get_seq_rc().substr(bp, kmer));
        std::cout << r2.get_seq_rc().substr(bp, kmer) << '\n';
        for (auto it = test.first; it != test.second; ++it) {
            spReadBase overlapped = check(r1, r2, (*it).second, bp, 0, 0);
            if (overlapped != nullptr) {
                std::cout << "Returning Overlapped\n";
                const Read &r = overlapped->get_read();
                std::cout << r.get_seq() << '\n';
                return overlapped;
            }
        }
    } 

    return nullptr;

}

spReadBase check_read(PairedEndRead &pe , const size_t maxMis, const size_t minOver) {
    
    Read &r1 = pe.non_const_read_one();
    Read &r2 = pe.non_const_read_two();

    if (r1.getLength() < r2.getLength()) {
        std::swap(r1, r2);
    }
    seqLookup mOne = readOneMap(r1.get_seq(), step);
    spReadBase overlapped = getOverlappedReads(r1, r2, mOne) ;
    const Read &r = overlapped->get_read();
    std::cout << r.get_seq() << '\n';
    if (overlapped != nullptr) {
       return overlapped;
    } 
}


template <class T, class Impl>
void helper_overlapper(InputReader<T, Impl> &reader, std::shared_ptr<OutputWriter> pe, std::shared_ptr<OutputWriter> se, Counter& counters, size_t hits) {
    
    while(reader.has_next()) {
        auto i = reader.next();
        ++counters["TotalRecords"];
        PairedEndRead* per = dynamic_cast<PairedEndRead*>(i.get());        
        if (per) {
            pe->write(*per);

        } else {
            SingleEndRead* ser = dynamic_cast<SingleEndRead*>(i.get());
            
            if (ser) {
                se->write(*ser);
            } else {
                throw std::runtime_error("Unknown read type");
            }
        }
    }
}

#endif

