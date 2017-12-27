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

void cut_trim(Read &r, size_t cut_left, size_t cut_right, size_t max_length) {
    if (max_length > 0) {
        if (r.getLength() > max_length) {
            // set left cut site to new cut_left
            r.setRCut(max_length);
        } else {
            // read is shorter than cut, so set new cut_right to cut_left, effectively 0 length
            r.setRCut(0);
        }
    }
    if (cut_right > 0) {
        if (r.getLength() > cut_right) {
            // set left cut site to new cut_left
            r.setRCut(r.getLength() - cut_right);
        } else {
            // read is shorter than cut, so set new cut_right to cut_left, effectively 0 length
            r.setRCut(0);
        }
    }
    if (cut_left > 0) {
        if (r.getLengthTrue() > cut_left) {
            // set left cut site to new cut_left
            r.setLCut(cut_left);
        } else {
            // read is shorter than cut, so set new cut_left to cut_right, effectively 0 length
            r.setLCut(r.getLengthTrue() + r.getLTrim());
        }
    }
}


template <class T, class Impl>
void helper_trim(InputReader<T, Impl> &reader, std::shared_ptr<OutputWriter> pe, std::shared_ptr<OutputWriter> se, TrimmingCounters& counters, size_t min_length , bool stranded, bool no_orphans, size_t r1_cut_left, size_t r1_cut_right, size_t r2_cut_left, size_t r2_cut_right, size_t max_length ){
    while(reader.has_next()) {
        auto i = reader.next();
        PairedEndRead* per = dynamic_cast<PairedEndRead*>(i.get());        
        if (per) {
            counters.input(*per);
            cut_trim( per->non_const_read_one(), r1_cut_left, r1_cut_right, max_length);
            cut_trim( per->non_const_read_two(), r2_cut_left, r2_cut_right, max_length);
            counters.output(*per);
            per->checkDiscarded(min_length);
            writer_helper(per, pe, se, stranded, no_orphans);
        } else {
            SingleEndRead* ser = dynamic_cast<SingleEndRead*>(i.get());
            if (ser) {
                counters.input(*ser);
                cut_trim( ser->non_const_read_one(), r1_cut_left, r1_cut_right, max_length);
                ser->checkDiscarded(min_length);
                counters.output(*ser);
                writer_helper(ser, pe, se);
            } else {
                throw std::runtime_error("Unknown read type");
            }
        }
    }

}

#endif
