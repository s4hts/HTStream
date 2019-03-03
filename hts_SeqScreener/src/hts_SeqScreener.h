/*The idea of this program is to screen out reads that "look like" a reference
 * this is accomplised by making a lookup table (hash table) that uses kmer binaries as their hash, and if enough Kmers
 * "hit" this table, it is considered close enough, and is discarded
  */
#ifndef SEQSCREENER_H
#define SEQSCREENER_H

#define BOOST_DYNAMIC_BITSET_DONT_USE_FRIENDS

#include "utils.h"
#include "ioHandler.h"
#include "counters.h"

#include <map>
#include <unordered_map>
#include <algorithm>
#include <bitset>
#include <unordered_set>
#include <boost/dynamic_bitset.hpp>
#include <boost/functional/hash.hpp>
#include <tuple>

extern template class InputReader<SingleEndRead, SingleEndReadFastqImpl>;
extern template class InputReader<PairedEndRead, PairedEndReadFastqImpl>;
extern template class InputReader<PairedEndRead, InterReadImpl>;
extern template class InputReader<ReadBase, TabReadImpl>;

class SeqScreenerCounters : public Counters {

public:
    std::string screen_file;
    uint64_t screen_bp = 0;
    uint64_t lookup_kmers = 0;
    std::vector<Label> screened_info;

    uint64_t Inverse = 0;
    uint64_t Record = 0;

    uint64_t SE_hits = 0;

    uint64_t PE_hits = 0;

    SeqScreenerCounters(const std::string &statsFile, bool force, bool appendStats, const std::string &program_name, const std::string &notes) : Counters::Counters(statsFile, force, appendStats, program_name, notes) {
        generic.push_back(std::forward_as_tuple("inverse", Inverse));
        generic.push_back(std::forward_as_tuple("record", Record));

        se.push_back(std::forward_as_tuple("SE_hits", SE_hits));

        pe.push_back(std::forward_as_tuple("PE_hits", PE_hits));
        // Need to and in screen_file name, however its a string and the vector can't have mixed types.
        //screened_info.push_back(std::forward_as_tuple("screenFile", screen_file));
        screened_info.push_back(std::forward_as_tuple("screenBP", screen_bp));
        screened_info.push_back(std::forward_as_tuple("lookupKmers", lookup_kmers));

    }

    void set_screeninfo(const std::string &screenFile, uint64_t sbp, uint64_t lkmers) {
        screen_file = screenFile;
        screen_bp = sbp;
        lookup_kmers = lkmers;
    }
    void set_inverse() {
        Inverse = 1;
    }
    void set_record() {
        Record = 1;
    }
    void inc_SE_hits() {
        ++SE_hits;
    }
    void inc_PE_hits() {
        ++PE_hits;
    }
    virtual void write_out() {

        initialize_json();

        write_labels(generic);
        write_sublabels("Screen_info", screened_info);
        write_sublabels("Single_end", se);
        write_sublabels("Paired_end", pe);

        finalize_json();
    }
};

class dbhash {
public:
    std::size_t operator() ( const boost::dynamic_bitset<>& bs) const {
        return boost::hash_value(bs.m_bits);
    }
};

typedef std::unordered_set < boost::dynamic_bitset<>, dbhash> kmerSet;

void setBitsBools(boost::dynamic_bitset<> &bs, size_t loc, bool set1, bool set2) {
    bs[loc] = set1;
    bs[loc + 1] = set2;
}

Read fasta_set_to_one_read(InputReader<SingleEndRead, FastaReadImpl> &faReader     ) {
    std::string s;

    while(faReader.has_next()) {
        auto r = faReader.next();
        s += (r->get_read().get_seq());
    }
    return Read(s, "", "All_Header");
}


std::pair <bool, bool> setBitsChar(char c) {

    switch (std::toupper(c)) {
        case 'A':
            return std::pair<bool, bool> (0, 0);
        case 'T':
            return std::pair<bool, bool> (1, 1);
        case 'C':
            return std::pair<bool, bool> (0, 1);
        case 'G':
            return std::pair<bool, bool> (1, 0);
        default:
            throw std::runtime_error("Unknown base pair in sequence " + c);
    }
}


/*The intent of this is to check  reads against the lookup table and return the number of hits
 * The function arguments are long because I did not want to recalculate those values each time
 * Especially reallocating the bitsets each time was costly.
 * Bitsets are safe, all positions will be overwritten without need to reset
 *
 * Two pairs of critical bitsets are at a work forwardLookup and forwardRest AND reverseLookup and reverseRest
 * the *Lookup varaibles are used as a 'prefix search'.
 *      We go to the prefix search location in the hashmap
 *          If diff (meaning the prefix DOES NOT contains all the kmer in the search) AND that location is not null
 *              check the vector at that location for *Rest (the rest of the data that isn't contained in Lookup)
 *              If the vector contains *Rest ++hit
 *
 *          else if !diff (lookup is the same size as KMER) AND that location is not null
 *              ++hit
 *
 * return number of hits seen */

unsigned int check_read( kmerSet &lookup, const Read &rb, const size_t bitKmer, const size_t lookup_loc, const size_t lookup_loc_rc, boost::dynamic_bitset<> &forwardLookup, boost::dynamic_bitset<> &reverseLookup) {

    std::string seq = rb.get_seq();

    unsigned int hits = 0;
    unsigned int current_added = 0;

    /*The flow of base paires will be
     * Forward add to 0 and 1 location then << 2
     * Once the ForwardLookup is full and start filling ForwardRest
     * ForwardRest will take the highest location of ForwardLookup and put it
     * in the lowest bit of ForwardRest then shift << 2*/

    /*The flow for RC is a bit different
     *The bits go into RestReverse first (if there is a diff is greater than 0)
     *They start of the highest bits in Reverse and flow downward >>= 2
     *Once RestReverse if "full", the lowest bits will be put into the Highest point of LookupRest and shifted down >>=2
     *if diff is 0 - then LookupRest is the only thing that is filled
     *These RC bits are filled iwth the RC( *bp )*/

    /*These the lookups are compared, the larger Lookup is taken and searched for,
     * if there is a "soft hit", initiate a search of the vector with the cooresponding Rest*/
    std::pair <bool, bool> bits;
    for (std::string::iterator bp = seq.begin(); bp < seq.end(); ++bp) { //goes through each bp of the read

        if (std::toupper(*bp) == 'N') { // N resets everythign
            current_added = 0;
        } else {
            reverseLookup >>= 2;
            forwardLookup <<=2; //No special condition needed for Forward Looup
            bits = setBitsChar(*bp);
            forwardLookup[lookup_loc] = bits.first;
            forwardLookup[lookup_loc+1] = bits.second;
            reverseLookup[lookup_loc_rc ] = !bits.first;
            reverseLookup[lookup_loc_rc + 1] = !bits.second;

            current_added += 2;

            if (current_added >= bitKmer) { //initiate hit search
                boost::dynamic_bitset<> &bs = forwardLookup > reverseLookup ? forwardLookup : reverseLookup;
                if (lookup.find(bs) != lookup.end()) {
                    ++hits;
                }
            }

         }
    }

    return hits;
}

template <class T, class Impl>
void helper_discard(InputReader<T, Impl> &reader, std::shared_ptr<OutputWriter> pe, std::shared_ptr<OutputWriter> se, SeqScreenerCounters& c, kmerSet &lookup, double hits, bool checkR2, size_t kmerSize, bool inverse = false, bool record = false) {

    /*These are set here so each read doesn't have to recalcuate these stats*/

    size_t bitKmer = kmerSize * 2;

    size_t lookup_loc = 0; // change bit 0 and 1 the << 2
    size_t lookup_loc_rc = bitKmer - 2; // change bit 31 and 30 then >> 2

    /*Particullary time consuming in each read*/
    boost::dynamic_bitset <> forwardLookup(bitKmer);
    boost::dynamic_bitset <> reverseLookup(bitKmer);

    if (inverse && !record) c.set_inverse();
    if (record) c.set_record();

    while(reader.has_next()) {
        auto i = reader.next();
        PairedEndRead* per = dynamic_cast<PairedEndRead*>(i.get());
        if (per) {
            c.input(*per);
            double val = check_read(lookup, per->get_read_one(), bitKmer, lookup_loc, lookup_loc_rc, forwardLookup, reverseLookup );
            val = val / ( per->get_read_one().getLength() - kmerSize);

            if (checkR2) {
                double val2 = check_read(lookup, per->get_read_two(), bitKmer, lookup_loc, lookup_loc_rc, forwardLookup, reverseLookup );

                val2 = val2 / (per->get_read_one().getLength() - kmerSize);
                val = std::max(val, val2);
            }

            if (val > hits) {
                c.inc_PE_hits();
            }
            if (val <= hits && !inverse && !record) {
                c.output(*per);
                writer_helper(per, pe, se, false);
            } else if (val > hits && inverse && !record) {
                c.output(*per);
                writer_helper(per, pe, se, false);
            } else if (record) {
                c.output(*per);
                writer_helper(per, pe, se, false);
            }

        } else {
            SingleEndRead* ser = dynamic_cast<SingleEndRead*>(i.get());

            if (ser) {
                c.input(*ser);
                double val = check_read(lookup, ser->get_read(), bitKmer, lookup_loc, lookup_loc_rc, forwardLookup, reverseLookup );
                val = val / ( ser->get_read().getLength() - kmerSize);

                if (val > hits) {
                    c.inc_SE_hits();
                }
                if (val <= hits && !inverse && !record) {
                    c.output(*ser);
                    writer_helper(ser, pe, se, false);
                } else if (val > hits && inverse && !record) {
                    c.output(*ser);
                    writer_helper(ser, pe, se, false);
                } else if (record) {
                    c.output(*ser);
                    writer_helper(ser, pe, se, false);
                }
            } else {
                throw std::runtime_error("Unknown read type");
            }
        }
    }
}

    /*The flow of base pairs will be
     * Forward add to 0 and 1 location then << 2
     * Once the ForwardLookup is full and start filling ForwardRest
     * ForwardRest will take the highest location of ForwardLookup and put it
     * in the lowest bit of ForwardRest then shift << 2*/

    /*The flow for RC is a bit different
     *The bits go into RestReverse first (if there is a diff is greater than 0)
     *They start of the highest bits in Reverse and flow downward >>= 2
     *Once RestReverse if "full", the lowest bits will be put into the Highest point of LookupRest and shifted down >>=2
     *if diff is 0 - then LookupRest is the only thing that is filled
     *These RC bits are filled iwth the RC( *bp )*/

    /*These the lookups are compared, the larger Lookup is taken and searched for,
     * if there is a "soft hit", initiate a search of the vector with the cooresponding Rest*/

void setLookup( kmerSet &lookup, Read &rb, size_t kmerSize) {
    /*This function is only called once, these values calculation are minimal cost so they are not in the function args*/
    /*Total size of kmer bits*/
    size_t bitKmer = kmerSize * 2;

    /*lookup locations*/
    size_t lookup_loc = 0;
    size_t lookup_loc_rc = bitKmer - 2;

    size_t current_added = 0;
    boost::dynamic_bitset <> forwardLookup(bitKmer);
    boost::dynamic_bitset <> reverseLookup(bitKmer);

    std::string seq = rb.get_seq();
    std::pair<bool, bool> bits;
    for (std::string::iterator bp = seq.begin(); bp != seq.end(); ++bp) {
        if (std::toupper(*bp) == 'A' || std::toupper(*bp) == 'C' || std::toupper(*bp) == 'G' || std::toupper(*bp) == 'T' ) { // N
            reverseLookup >>= 2;
            forwardLookup <<=2;

            bits = setBitsChar(*bp);
            forwardLookup[lookup_loc] = bits.first;
            forwardLookup[lookup_loc + 1] = bits.second;
            reverseLookup[lookup_loc_rc] = !bits.first;
            reverseLookup[lookup_loc_rc + 1] = !bits.second;
            current_added += 2;
            if (current_added >= bitKmer) {
                boost::dynamic_bitset<> &bs = forwardLookup > reverseLookup ? forwardLookup : reverseLookup;
                lookup.insert(bs);
            }
         } else {
           current_added = 0;
         }
    }

}

#endif