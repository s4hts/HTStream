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

typedef std::unordered_map<std::string, size_t> Counter;

template <class T, class Impl>
void helper_trim(InputReader<T, Impl> &reader, std::shared_ptr<OutputWriter> pe, std::shared_ptr<OutputWriter> se, Counter& counters, size_t min_length, size_t cut_size, bool stranded, bool no_left, bool no_right) {
    
    while(reader.has_next()) {
        auto i = reader.next();
        ++counters["TotalRecords"];
        PairedEndRead* per = dynamic_cast<PairedEndRead*>(i.get());        
        if (per) {
            Read &rb1 = per->non_const_read_one();
            Read &rb2 = per->non_const_read_one();
            if (!no_left) {
                rb1.setLCut(cut_size);
                rb2.setLCut(cut_size);
            }
            if (!no_right) {
                rb1.setRCut(rb1.getLength() - cut_size);
                rb2.setRCut(rb2.getLength() - cut_size);
            }
            per->checkDiscarded(min_length);
            writer_helper(per, pe, se, stranded, counters);
        } else {
            SingleEndRead* ser = dynamic_cast<SingleEndRead*>(i.get());
            
            if (ser) {
                Read &rb = per->non_const_read_one();
                if (!no_left) {
                    rb.setLCut(cut_size);
                } 
                if (!no_right) {
                    rb.setRCut((per->non_const_read_two()).getLength() - cut_size);
                }
                ser->checkDiscarded(min_length);
                writer_helper(ser, pe, se, stranded, counters);
            } else {
                throw std::runtime_error("Unknown read type");
            }
        }
    }

}

#endif
