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
#include "utils.h"

int dist(int x, int y) {
    return y - x;
}

/*This is a O(N) time algorithm
 * it will search for the longest base pair segment that has
 * no N's within it*/
void trim_n(Read &rb, size_t min_length) {
    
    const std::string &seq = rb.get_seq();
    int len = seq.length() - 1;
    int bestLeft = 0, currentLeft = 0, bestRight = 0;

    for (int i = 0; i <= len; ++i) {
        if (seq[i] == 'N') {
            currentLeft = i + 1;
        } else if (dist(bestLeft, bestRight) <= dist(currentLeft, i)) {
            bestRight = i;
            bestLeft = currentLeft;
        }
    }
    rb.setLCut(bestLeft);
    rb.setRCut(bestRight+1);
}


/*Removes all Ns (ambiguity base) from a read*/
template <class T, class Impl>
void helper_trim(InputReader<T, Impl> &reader, std::shared_ptr<OutputWriter> pe, std::shared_ptr<OutputWriter> se, Counter& counters, bool stranded, size_t min_length) {
    
    while(reader.has_next()) {
        auto i = reader.next();
        ++counters["TotalRecords"];
        PairedEndRead* per = dynamic_cast<PairedEndRead*>(i.get());        
        if (per) {
            trim_n(per->non_const_read_one(), min_length);            
            trim_n(per->non_const_read_two(), min_length); 
            per->checkDiscarded(min_length);
            writer_helper(per, pe, se, stranded, counters);
        } else {
            SingleEndRead* ser = dynamic_cast<SingleEndRead*>(i.get());
            
            if (ser) {
                trim_n(ser->non_const_read_one(), min_length);
                ser->checkDiscarded(min_length);
                writer_helper(ser, pe, se, stranded, counters);
            } else {
                throw std::runtime_error("Unknow read type");
            }
        }
    }

}

#endif
