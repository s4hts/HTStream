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

class AdapterTrimmerCounters : public Counters {

public:
    std::vector<unsigned long long int> insertLength;

    AdapterTrimmerCounters() {
        Common();
        c["lins"] = 0;
        c["sins"] = 0;
        c["SE_Discard"] = 0;
        c["R1_Discard"] = 0;
        c["R2_Discard"] = 0;
        c["R1_Adapter_Trim"] = 0;
        c["R2_Adapter_Trim"] = 0;
        c["Nolins"] = 0;
        insertLength.resize(1);
    }
    
    void output(SingleEndRead &ser) {
        if (ser.non_const_read_one().getDiscard()) {
            ++c["SE_Discard"];
        } else {
            ++c["SE_Out"];
        }

    }

    void output(PairedEndRead &per, unsigned int overlapped) {
        //test lin or sin
        Read &one = per.non_const_read_one();
        Read &two = per.non_const_read_two();

        if (overlapped) {

            if (one.getLengthTrue() > overlapped || two.getLengthTrue() > overlapped ) {
                ++c["sins"]; //adapters must be had (short insert)
            } else {
                ++c["lins"]; //must be a long insert
            }
            if ( overlapped + 1 > insertLength.size() ) {
                insertLength.resize(overlapped + 1);
            }
            ++insertLength[overlapped];
        } else {
            ++c["Nolins"]; //lin
        }

        if (!one.getDiscard() && !two.getDiscard() ) {
            ++c["PE_Out"];
        } else if (one.getDiscard() && two.getDiscard()) {
            ++c["R1_Discard"];
            ++c["R2_Discard"];
        } else if (one.getDiscard()) {
            ++c["R1_Discard"];
        } else if (two.getDiscard()) {
            ++c["R2_Discard"];
        }

    }

    void write_out(std::string &statsFile, bool appendStats, std::string program_name, std::string notes) {

        std::ifstream testEnd(statsFile);
        int end = testEnd.peek();
        testEnd.close();

        std::fstream outStats;
        bool first = true;

        if (appendStats && end != -1) {
            outStats.open(statsFile, std::ios::out|std::ios::in); //overwritte
            outStats.seekp(-1, std::ios::end );
            outStats << "\n,\"" << program_name << "_" << getpid()  << "\": {\n";
        } else {
            outStats.open(statsFile, std::ios::out); //overwritt
            outStats << "{\n \"" << program_name << "_" << getpid() <<  "\": {\n";
        }
        outStats << "\"Notes\" : \"" << notes << "\",\n";
        for (const auto name : c) {
            if (first) {
                first = false;
            } else {
                outStats << ",\n"; //make sure json format is kept
            }
            outStats << "\"" << name.first << "\" : " << name.second; //it will get the comma in conditionals tatement about
        }
        for (int i = 1; i < insertLength.size(); ++i) {
            outStats << ",\n"; //make sure json format is kept
            outStats << '"' << i << '"' << " : "  << insertLength[i];
        }
        outStats << "\n}";
        outStats << "\n}";
        outStats.flush();

    }

};

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

char rc(const char bp) {
    if (bp == 'C') {
        return 'G';
    } else if (bp == 'G') {
        return 'C';
    } else if (bp == 'T') {
        return 'A';
    } else if (bp == 'A') {
        return 'T';
    } else {
        return 'N';
    }
}
/*If adapater trimming is turned on that means adapter trimming and do not overlap
 * so trim adapter, but don't worry about the overlap.
 * however we still need to change the overlap
 * 
 * Within the overlap if they are the same bp, then add q scores
 * If they are different bp, subtract q scores and take the larger quality bp*/ 
unsigned int checkIfOverlap(Read &r1, Read &r2, size_t loc1, size_t loc2, const double misDensity, size_t minOverlap) {
    size_t minLoc = std::min(loc1, loc2);
    size_t loc1_t = loc1 - minLoc;
    size_t loc2_t = loc2 - minLoc;
    size_t maxLoop = std::min(r1.getLength() - loc1_t, r2.getLength() - loc2_t);
    size_t maxMis = static_cast<size_t>(maxLoop * misDensity);

    if (maxLoop <= minOverlap) {
        return 0;
    }

    size_t read1_bp;
    size_t read2_bp;

    const std::string &seq1 = r1.get_seq();
    const std::string &seq2 = r2.get_seq();
    const std::string &qual1 = r1.get_qual();
    const std::string &qual2 = r2.get_qual();

    const size_t seq2_length =  seq2.length() - 1;
    
    std::string finalSeq;
    std::string finalQual;
    size_t misMatches = 0;

    char bp;
    char qual;

    unsigned int insert_length = maxLoop;
    for (size_t i = 0; i < maxLoop; ++i) {
        read1_bp = loc1_t + i;
        read2_bp = seq2_length - (loc2_t + i);
        
        if (seq1[read1_bp] != rc(seq2[read2_bp]) )  {
            bp = qual1[read1_bp] > qual2[read2_bp] ? seq1[read1_bp] : rc(seq2[read2_bp]);
            qual = static_cast<char>(std::max(qual1[read1_bp] - qual2[read2_bp] + 33, 1 + 33));
            
            r1.changeSeq(read1_bp, bp);
            r1.changeQual(read1_bp, qual);
            //Something clever with RC
            r2.changeSeq( read2_bp, rc(bp) );
            r2.changeQual(read2_bp, qual);
            
            if (++misMatches > maxMis) {
                /*Not valid match*/
                return 0;
            }
        }
    }

    if ( loc1_t >= loc2_t ) { //no adapter but need to figure out sides
        insert_length += loc1_t; 
    } else {
        r1.setLCut(loc1_t); //sin
    }
    
    if (r1.getLength() - loc1_t < r2.getLength()) {
        insert_length += r2.getLength() - maxLoop; // lin
    } else {
        r2.setRCut(maxLoop); //sin
    }

    return insert_length;
}

/*Because of the way overlapping works, you only need to check the ends of the shorter read*/
unsigned int getOverlappedReads(Read &r1, Read &r2, const seqLookup &seq1Map,  const double misDensity, const size_t &minOver, const size_t &checkLengths, const size_t kmer) {

    std::string seq2 = r2.get_seq_rc();
    for (size_t bp = 0; bp < checkLengths; ++bp) {
        /*Do a quick check if the shorter read kmer shows up in longer read (read 2)
         * If it does, then try the brute force approach*/
        auto test = seq1Map.equal_range(seq2.substr(bp, kmer));
        for (auto it = test.first; it != test.second; ++it) {
            unsigned int overlapped = checkIfOverlap(r1, r2, it->second, bp, misDensity, minOver);
            if (overlapped) {
                return overlapped;
            }
        }
    } 

    for (size_t bp = seq2.length() - (checkLengths + kmer); bp < seq2.length() - kmer ; ++bp) {
        /*Do a quick check if the shorter read kmer shows up in longer read (read 2)
         * If it does, then try the brute force approach*/
        auto test = seq1Map.equal_range(seq2.substr(bp, kmer));
        for (auto it = test.first; it != test.second; ++it) {
            unsigned int overlapped = checkIfOverlap(r1, r2, it->second, bp, misDensity, minOver);
            if (overlapped) {
                return overlapped;
            }
        }
    }
    return 0;

}

unsigned int check_read(PairedEndRead &pe , const double misDensity, const size_t &minOver, const size_t &checkLengths, const size_t kmer, const size_t kmerOffset, const size_t minLength) {
    
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
    //we need to check it overlapper is greater than min length;

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
void helper_overlapper(InputReader<T, Impl> &reader, std::shared_ptr<OutputWriter> pe, std::shared_ptr<OutputWriter> se, AdapterTrimmerCounters &counter, const double misDensity, const size_t minOver, const bool stranded, const size_t min_length, const size_t checkLengths, const size_t kmer, const size_t kmerOffset, bool no_orphan = false ) {
    
    while(reader.has_next()) {
        auto i = reader.next();
        PairedEndRead* per = dynamic_cast<PairedEndRead*>(i.get());        
        if (per) {
            counter.input(*per);
            unsigned int overlapped = check_read(*per, misDensity, minOver, checkLengths, kmer, kmerOffset, min_length);
            per->checkDiscarded(min_length);    
            writer_helper(per, pe, se, stranded, no_orphan); 
            counter.output(*per, overlapped);
        } else {
            SingleEndRead* ser = dynamic_cast<SingleEndRead*>(i.get());
            
            if (ser) {
                counter.input(*ser);
                ser->checkDiscarded(min_length);
                writer_helper(ser, pe, se, stranded);
                counter.output(*ser);
            } else {
                throw std::runtime_error("Unknown read type");
            }
        }
    }
}

#endif

