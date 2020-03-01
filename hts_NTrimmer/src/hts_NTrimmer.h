#ifndef N_TRIM_H
#define N_TRIM_H
//  this is so we can implment hash function for dynamic_bitset
#define BOOST_DYNAMIC_BITSET_DONT_USE_FRIENDS

#include "ioHandler.h"
#include <map>
#include <unordered_map>
#include <boost/dynamic_bitset.hpp>
#include <boost/functional/hash.hpp>
#include <algorithm>
#include "utils.h"

extern template class InputReader<SingleEndRead, SingleEndReadFastqImpl>;
extern template class InputReader<PairedEndRead, PairedEndReadFastqImpl>;
extern template class InputReader<PairedEndRead, InterReadImpl>;
extern template class InputReader<ReadBase, TabReadImpl>;

size_t dist(size_t x, size_t y) {
    assert(x <= y);
    return y - x;
}

/*This is a O(N) time algorithm
 * it will search for the longest base pair segment that has
 * no N's within it*/
void trim_n(Read &rb, bool exclude) {

    std::string seq = rb.get_seq();
    size_t bestLeft = 0, currentLeft = 0, bestRight = 0;
    size_t i = 0;

    for (std::string::iterator it = seq.begin(); it != seq.end(); ++it) {
        i = static_cast<size_t>(it - seq.begin() );

        if (*it == 'N') {
            if (exlude) {
              rb.setLCut(1);
              rb.setRCut(0);
              break;
            } else {
                currentLeft = i + 1;
            }
        } else if (dist(bestLeft, bestRight) < dist(currentLeft, i)) {
            bestRight = i;
            bestLeft = currentLeft;
        }
    }
    rb.setLCut(bestLeft);
    rb.setRCut(bestRight+1);
}


/*Removes all Ns (ambiguity base) from a read*/
template <class T, class Impl>
void helper_trim(InputReader<T, Impl> &reader, std::shared_ptr<OutputWriter> pe, std::shared_ptr<OutputWriter> se, TrimmingCounters& counters, size_t min_length , bool stranded, bool no_orphans, bool exclude) {

    while(reader.has_next()) {
        auto i = reader.next();
        PairedEndRead* per = dynamic_cast<PairedEndRead*>(i.get());
        if (per) {
            counters.input(*per);
            trim_n(per->non_const_read_one(), exclude);
            trim_n(per->non_const_read_two(), exclude);
            per->checkDiscarded(min_length);
            writer_helper(per, pe, se, stranded, no_orphans);
            counters.output(*per, no_orphans);
        } else {
            SingleEndRead* ser = dynamic_cast<SingleEndRead*>(i.get());
            if (ser) {
                counters.input(*ser);
                trim_n(ser->non_const_read_one());
                ser->checkDiscarded(min_length);
                writer_helper(ser, pe, se, false, false);
                counters.output(*ser);
            } else {
                throw std::runtime_error("Unknow read type");
            }
        }
    }

}

#endif
