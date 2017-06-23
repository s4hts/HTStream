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

void cut_trim(Read &r, bool no_left, bool no_right, size_t cut_size) {
    if (!no_left) {
        r.setLCut(cut_size);
    }
    if (!no_right) {
        if (r.getLength() > cut_size) {
            r.setRCut(r.getLength() - cut_size);
        }
    }
}


template <class T, class Impl>
void helper_trim(InputReader<T, Impl> &reader, std::shared_ptr<OutputWriter> pe, std::shared_ptr<OutputWriter> se, TrimmingCounters& counters, size_t min_length, size_t cut_size, bool stranded, bool no_left, bool no_right, bool no_orphans) {
    while(reader.has_next()) {
        auto i = reader.next();
        PairedEndRead* per = dynamic_cast<PairedEndRead*>(i.get());        
        if (per) {
            counters.input(*per);
            cut_trim( per->non_const_read_one(), no_left, no_right, cut_size);
            cut_trim( per->non_const_read_two(), no_left, no_right, cut_size);
            counters.output(*per);
            per->checkDiscarded(min_length);
            writer_helper(per, pe, se, stranded, no_orphans);
        } else {
            SingleEndRead* ser = dynamic_cast<SingleEndRead*>(i.get());
            if (ser) {
                counters.input(*ser);
                cut_trim( ser->non_const_read_one(), no_left, no_right, cut_size);
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
