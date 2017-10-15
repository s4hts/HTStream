#include "overlapping_utils.h"

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
unsigned int checkIfOverlap(Read &r1, Read &r2, size_t loc1, size_t loc2, const double misDensity, size_t minOverlap) {
    size_t minLoc = std::min(loc1, loc2);
    int loc1_t = loc1 - minLoc;
    int loc2_t = loc2 - minLoc;
    int r1_len = r1.getLength() - 1;
    int r2_len = r2.getLength() - 1;

    size_t maxLoop = std::min(r1.getLength() - loc1_t, r2.getLength() - loc2_t);
    size_t maxMis = static_cast<size_t>(maxLoop * misDensity);
    
    const std::string &seq1 = r1.get_seq();
    const std::string &seq2 = r2.get_seq_rc();

    auto i1 = seq1.begin();
    std::advance(i1, loc1_t);
    auto i2 = seq2.begin();
    std::advance(i2, loc2_t);
    
    if (maxLoop <= minOverlap || !threshold_mismatches(i1, i2, maxLoop, maxMis) ) {
        return 0;
    }

    unsigned int insert_length = maxLoop;
   

    //This does not need to be check for the engulf case
    //R2 will be R1  
    if (loc1_t >= loc2_t) {
        insert_length += loc1_t;
    }

    //LIN
    if (r1_len - loc1_t <= r2_len - loc2_t) {
        insert_length += (r2_len - maxLoop);
    } else {
        //SIN
        //insert length will be the overlap length
    }

   return insert_length;
}

/*Because of the way overlapping works, you only need to check the ends of the shorter read*/
unsigned int getOverlappedReads(Read &r1, Read &r2, const seqLookup &seq1Map,  const double misDensity, const size_t &minOver, const size_t &checkLengths, const size_t kmer) {

    std::string seq2 = r2.get_seq_rc();
    for (size_t bp = seq2.length() - (checkLengths + kmer); bp < seq2.length() - kmer ; ++bp) {
        /*Do a quick check if the shorter read kmer shows up in longer read (read 2)
         * If it does, then try the brute force approach*/
        auto test = seq1Map.equal_range(seq2.substr(bp, kmer));
        for (auto it = test.first; it != test.second; ++it) {
            unsigned int overlapped = checkIfOverlap(r1, r2, it->second, bp, misDensity, minOver);
            if (overlapped) {
                return overlapped + 1;
            }
        }
     }
     for (size_t bp = 0; bp < checkLengths; ++bp) {
        /*Do a quick check if the shorter read kmer shows up in longer read (read 2)
         * If it does, then try the brute force approach*/
        auto test = seq1Map.equal_range(seq2.substr(bp, kmer));
        for (auto it = test.first; it != test.second; ++it) {
            unsigned int overlapped = checkIfOverlap(r1, r2, it->second, bp, misDensity, minOver);
            if (overlapped) {
                return overlapped + 1;
            }
        }
    } 

   return 0;

}

unsigned int getInsertSize(PairedEndRead &pe , const double misDensity, const size_t &minOver, const size_t &checkLengths, const size_t kmer, const size_t kmerOffset, const size_t minLength) {
    
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
    unsigned int overlapped = getOverlappedReads(r1, r2, mOne, misDensity, minOver, checkLengths, kmer) ;
    if (swapped) {
        std::swap(r1, r2);
        r2.set_read_rc();
    }
    //we need to check it overlapper is greater than min length;

    return overlapped;
}

