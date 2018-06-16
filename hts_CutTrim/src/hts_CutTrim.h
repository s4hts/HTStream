#ifndef CUT_TRIM_H
#define CUT_TRIM_H
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

/*
cut_trim expected behavior:
cut_trim is expected to be one of the first apps used in a preprocessing pipeline
such as when you want to statically trim off X bases from the end of a read (cut_right), due to
quality, and/or to trim off primer sequence from the beginning of a read (cut_left). As such both
cut_left and cut_right should act on the original length of the read.

max_length then prevents the final cut length from being greater than max_length applying an
additional cut to end of the read if applicable. If one wished to simulate smaller sized reads
can run cut_trim with max_length (no cut to left or right) to effectively reduce the size of reads.
 */
void cut_trim(Read &r, size_t cut_left, size_t cut_right, size_t max_length) {
    if (cut_left) {
        r.setLCut(cut_left);
    }
    if (cut_right) {
        r.setRCut(r.getLength() - cut_right);
    }
    if (max_length && max_length < r.getLengthTrue()) {
        r.setRCut(max_length + r.getLTrim());
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
            per->checkDiscarded(min_length);
            writer_helper(per, pe, se, stranded, no_orphans);
            counters.output(*per, no_orphans);
        } else {
            SingleEndRead* ser = dynamic_cast<SingleEndRead*>(i.get());
            if (ser) {
                counters.input(*ser);
                cut_trim( ser->non_const_read_one(), r1_cut_left, r1_cut_right, max_length);
                ser->checkDiscarded(min_length);
                writer_helper(ser, pe, se, false, false);
                counters.output(*ser);
            } else {
                throw std::runtime_error("Unknown read type");
            }
        }
    }

}

#endif
