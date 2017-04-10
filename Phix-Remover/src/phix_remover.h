/*The idea of this program is to screen out reads that "look like" a reference
 * this is accomplised by making a lookup table (hash table) that uses kmer binaries as their hash, and if enough Kmers
 * "hit" this table, it is considered close enough, and is discarded
  */
#ifndef PHIX_REMOVER_H
#define PHIX_REMOVER_H

#define BOOST_DYNAMIC_BITSET_DONT_USE_FRIENDS

#include "utils.h"
#include "ioHandler.h"

#include <map>
#include <unordered_map>
#include <algorithm>
#include <bitset>
#include <unordered_set>
#include <boost/dynamic_bitset.hpp>
#include <boost/functional/hash.hpp>

/*Simple class that will be added to the lookup table
 * Each soft hit that have a Rest Bitset will add to the 
 * vector at the location of the soft hit*/
class Lookup {
public:
    Lookup(boost::dynamic_bitset<> bs) {
        vecBits.push_back(bs);
    }
    Lookup() { //empty constructor
    }
    ~Lookup() {
    }
    

    void insert(boost::dynamic_bitset<> bs) {
        std::vector<boost::dynamic_bitset<> >::iterator it= std::lower_bound(vecBits.begin(), vecBits.end(), bs);
        if (*it != bs) {
            vecBits.insert(it, bs);
        }
    }

    void print() {
        for (auto &a : vecBits) {
           std::cout << a << "\n";
        }
    }

    unsigned int check(boost::dynamic_bitset<> &bs) {
        return std::binary_search(vecBits.begin(), vecBits.end(), bs);
    }

    std::vector< boost::dynamic_bitset<> >  vecBits;

};




class dbhash {
public:
    std::size_t operator() ( const boost::dynamic_bitset<>& bs) const {
        return boost::hash_value(bs.m_bits);
    }
};

typedef std::unordered_set < boost::dynamic_bitset<>, dbhash> kmerSet;
typedef std::unique_ptr< std::shared_ptr<Lookup>[] >firstLookup;
typedef std::shared_ptr<Lookup > *firstLookupPointer;

void setBitsBools(boost::dynamic_bitset<> &bs, size_t loc, bool set1, bool set2) {
    bs[loc] = set1;
    bs[loc + 1] = set2;
}

int setBitsChar(boost::dynamic_bitset<> &lookup, size_t loc, char c, bool rc) {

    if (c == 'A') {
        lookup[loc] = (0 ^ rc);
        lookup[loc+1] = (0 ^ rc);
    } else if (c == 'T') {
        lookup[loc] = (1 ^ rc);
        lookup[loc+1] = (1 ^ rc);
    } else if (c == 'C') {
        lookup[loc] = (1 ^ rc);
        lookup[loc+1] = (0 ^ rc);
    } else if (c == 'G') {
        lookup[loc] = (0 ^ rc)  ;
        lookup[loc+1] = (1 ^ rc);
    } else {
        return 0; //N
    }
    return 1;

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


unsigned int check_read( firstLookupPointer lookup, const Read &rb, const size_t kmerSize, const size_t lookupKmer,
        const size_t rest_loc, const size_t rest_loc_rc, const size_t bitKmer, const size_t bitKmerLookupSize, const size_t lookup_loc, const size_t lookup_loc_rc,
        const size_t diff, boost::dynamic_bitset<> &forwardLookup, boost::dynamic_bitset<> &reverseLookup,
        boost::dynamic_bitset<> &forwardRest, boost::dynamic_bitset<> &reverseRest) //This arg list is kind of crazy, but I don't want to recalc each read
         {

    std::string seq = rb.get_seq();

    unsigned long ulLookup = 0;
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

    for (std::string::iterator bp = seq.begin(); bp < seq.end(); ++bp) { //goes through each bp of the read
         
        if (*bp == 'N' ) { // N resets everythign
            current_added = 0;
        } else {
            if (diff) {
                    forwardRest <<= 2; //adds read in the "Normal order"
                    reverseLookup >>= 2; //adds read in the "Reverse Complement" order - in the RC case
                    //Forward follows the progression lookup -> Rest 
                    setBitsBools(forwardRest, rest_loc, forwardLookup[bitKmerLookupSize - 2], forwardLookup[bitKmerLookupSize - 2 + 1]);
                    //Reverse goes Rest  -> Lookup
                    setBitsBools(reverseLookup,lookup_loc_rc, reverseRest[0], reverseRest[1]); //Follows the flow of RC
                    reverseRest >>= 2;
                    setBitsChar(reverseRest, rest_loc_rc, *bp, true);
            } else {
                    reverseLookup >>= 2;
                    setBitsChar(reverseLookup,lookup_loc_rc, *bp, true); //Just add RC to lookup if Rest doesn't exist
            }

            forwardLookup <<=2; //No special condition needed for Forward Looup
            setBitsChar(forwardLookup, lookup_loc, *bp, false);
            current_added += 2;
           
            if (current_added >= bitKmer) { //initiate hit search
                boost::dynamic_bitset<> &bs = forwardLookup > reverseLookup ? forwardRest : reverseRest;
                ulLookup = forwardLookup > reverseLookup ? forwardLookup.to_ulong() : reverseLookup.to_ulong();
                if ( lookup[ulLookup ] ) {
                    if (diff) {
                        hits += lookup[ulLookup]->check(bs); //search vector
                    } else {
                        ++hits;//just check if it isn't null
                    }
                }
            } 
            
         }
    }

    return hits; 
}

template <class T, class Impl>
void helper_discard(InputReader<T, Impl> &reader, std::shared_ptr<OutputWriter> pe, std::shared_ptr<OutputWriter> se, Counter& c, firstLookupPointer lookup, double hits, bool checkR2, size_t kmerSize, size_t kmerLookupSize) {

    /*These are set here so each read doesn't have to recalcuate these stats*/

    size_t bitKmer = kmerSize * 2;
    size_t bitKmerLookupSize = kmerLookupSize * 2;

    size_t lookup_loc = 0; // change bit 0 and 1 the << 2
    size_t lookup_loc_rc = bitKmerLookupSize - 2; // change bit 31 and 30 then >> 2
    size_t rest_loc = 0;
    size_t rest_loc_rc = 0;

    size_t diff = 0;
    /*Particullary time consuming in each read*/
    boost::dynamic_bitset <> forwardLookup(bitKmerLookupSize);
    boost::dynamic_bitset <> reverseLookup(bitKmerLookupSize);

    boost::dynamic_bitset <> forwardRest;
    boost::dynamic_bitset <> reverseRest;

    /*If lookup is less the Rest then we have a special case*/
    if (bitKmer > bitKmerLookupSize) {
        diff = bitKmer - bitKmerLookupSize;
        forwardRest = boost::dynamic_bitset<>(  diff  );
        reverseRest = boost::dynamic_bitset<>(  diff );
        rest_loc = 0;
        rest_loc_rc = diff - 2;
    }
    

    while(reader.has_next()) {
        auto i = reader.next();
        PairedEndRead* per = dynamic_cast<PairedEndRead*>(i.get());
        if (per) {

            double val = check_read(lookup, per->get_read_one(), kmerSize, kmerLookupSize, rest_loc, rest_loc_rc, bitKmer, bitKmerLookupSize, lookup_loc, lookup_loc_rc, diff, forwardLookup, reverseLookup, forwardRest, reverseRest );
            val = val / ( per->get_read_one().getLength() - kmerSize);
           
            if (checkR2) {
                double val2 = check_read(lookup, per->get_read_two(), kmerSize, kmerLookupSize, rest_loc, rest_loc_rc, bitKmer, bitKmerLookupSize, lookup_loc, lookup_loc_rc, diff, forwardLookup, reverseLookup, forwardRest, reverseRest );
                
                val2 = val2 / (per->get_read_one().getLength() - kmerSize);
                val = std::max(val, val2);
            }

            if (val <= hits) {
                writer_helper(per, pe, se, false, c);
            } else {
                ++c["R1_Discarded"];
                ++c["R2_Discarded"];
            }

        } else {
            SingleEndRead* ser = dynamic_cast<SingleEndRead*>(i.get());

            if (ser) {
                double val = check_read(lookup, per->get_read(), kmerSize, kmerLookupSize, rest_loc, rest_loc_rc, bitKmer, bitKmerLookupSize, lookup_loc, lookup_loc_rc, diff, forwardLookup, reverseLookup, forwardRest, reverseRest );
                if (val <= hits) {
                    writer_helper(ser, pe, se, false, c);
                } else {
                    ++c["SE_Discarded"];
                }
            } else {
                throw std::runtime_error("Unknown read type");
            }
        }
    }
}

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

void setLookup( firstLookupPointer lookup, Read &rb, size_t kmerSize, size_t lookupKmer) {
    /*This function is only called once, these values calculation are minimal cost so they are not in the function args*/
    /*Total size of kmer bits*/
    size_t bitKmer = kmerSize * 2;
    size_t bitKmerLookupSize = lookupKmer * 2;
   
    /*lookup locations*/ 
    size_t lookup_loc = 0; 
    size_t lookup_loc_rc = bitKmerLookupSize - 2;
    size_t rest_loc = 0;
    size_t rest_loc_rc = 0;
    
    size_t current_added = 0;
    size_t diff = 0;
    boost::dynamic_bitset <> forwardLookup(bitKmerLookupSize);
    boost::dynamic_bitset <> reverseLookup(bitKmerLookupSize);

    boost::dynamic_bitset <> forwardRest;
    boost::dynamic_bitset <> reverseRest;

    if (bitKmer > bitKmerLookupSize) {
        diff = bitKmer - bitKmerLookupSize;
        forwardRest = boost::dynamic_bitset<>(  diff  );
        reverseRest = boost::dynamic_bitset<>(  diff );
        rest_loc = 0;
        rest_loc_rc = diff - 2;
    }
    
    std::string seq = rb.get_seq();
    unsigned long ulLookup = 0;

    for (std::string::iterator bp = seq.begin(); bp < seq.end(); ++bp) {
        if (*bp == 'N' ) { // N
            current_added = 0;
        } else {
            if (diff) {
                    forwardRest <<= 2;
                    reverseLookup >>= 2;
                    //Forward follows the progression lookup -> Rest 
                    setBitsBools(forwardRest, rest_loc, forwardLookup[bitKmerLookupSize - 2], forwardLookup[bitKmerLookupSize - 2 + 1]);
                    //Reverse goes Rest  -> Lookup
                    setBitsBools(reverseLookup,lookup_loc_rc, reverseRest[0], reverseRest[1]);
                    reverseRest >>= 2;
                    setBitsChar(reverseRest, rest_loc_rc, *bp, true);
                   // std::cout << reverseRest << "<- RRest " << reverseLookup << " <- RLooup\n";
            } else {
                    reverseLookup >>= 2;
                    setBitsChar(reverseLookup,lookup_loc_rc, *bp, true);
            }

            forwardLookup <<=2;
            setBitsChar(forwardLookup, lookup_loc, *bp, false);
            current_added += 2;

            if (current_added >= bitKmer) {
                boost::dynamic_bitset<> &bs = forwardLookup > reverseLookup ? forwardRest : reverseRest;
                ulLookup = forwardLookup > reverseLookup ? forwardLookup.to_ulong() : reverseLookup.to_ulong();
                if ( !lookup[ulLookup ] ) { //null location (create class at that object
                    if (diff) {
                        lookup[ulLookup ] = std::make_shared<Lookup>(bs); //insert rest into that locations vector
                    } else {
                        lookup[ulLookup ] = std::make_shared<Lookup>(); //just create an empty class since there is no rest class
                    }
                } else if (diff) {
                    lookup[ulLookup]->insert(bs);
                }
            } 
            
         }
    }

}


#endif
