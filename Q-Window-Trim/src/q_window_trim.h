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
#include "counters.h"

extern template class InputReader<SingleEndRead, SingleEndReadFastqImpl>;
extern template class InputReader<PairedEndRead, PairedEndReadFastqImpl>;
extern template class InputReader<PairedEndRead, InterReadImpl>;
extern template class InputReader<ReadBase, TabReadImpl>;

void trim_left(Read &rb, size_t qual_threshold, size_t window_size) {

    std::string qual = rb.get_qual();
    size_t current_sum = 0;
    size_t cut = 0;
    size_t i = 0;
    size_t sum_qual = 0;
    size_t final_qual_threshold = qual_threshold * window_size; //might save a bit of time

    for (std::string::iterator it = qual.begin() ; it != qual.end() ; ++it) {
        i = static_cast<size_t>(it - qual.begin()) + 1;
        current_sum += static_cast<size_t>(*it);
        if (i > window_size) { //once we hit window size, subtract the first value off
            current_sum -= static_cast<size_t>( *(it - static_cast<long>(window_size) ) );
            sum_qual = final_qual_threshold;
        } else {
            sum_qual = i * qual_threshold;
        }

        if (current_sum >= sum_qual) {
            break;
        } else {
            cut = i ;
        }
    }
    
    if (current_sum < sum_qual) {
        rb.setLCut(qual.length() - 1);
    } else {
        rb.setLCut(cut);
    }

}

void trim_right(Read &rb, size_t qual_threshold, size_t window_size) {

    std::string qual = rb.get_qual();
    size_t len = qual.length();
    size_t current_sum = 0;
    size_t cut = 0;
    size_t i = 0;
    size_t sum_qual = 0;
    size_t final_qual_threshold = qual_threshold * window_size; //might save a bit of time

    for (std::string::reverse_iterator it = qual.rbegin() ; it != qual.rend() ; ++it) {
        i = static_cast<size_t>(it - qual.rbegin()) + 1;
        current_sum += static_cast<size_t>(*it);
        if (i > window_size) { //once we hit window size, subtract the first value off
            current_sum -= static_cast<size_t>( *(it - static_cast<long>(window_size)) );
            sum_qual = final_qual_threshold;
        } else {
            sum_qual = i * qual_threshold;
        }
        
        if (current_sum >= sum_qual) {
            break;
        } else {
            cut = i;
        }

    }

    if (current_sum < sum_qual || cut > len ) { // will protect from cut - length causing an underflow
        rb.setRCut(0);
    } else {
        rb.setRCut(len - cut);
    }
 
}

template <class T, class Impl>
void helper_trim(InputReader<T, Impl> &reader, std::shared_ptr<OutputWriter> pe, std::shared_ptr<OutputWriter> se, TrimmingCounters &counters, size_t min_length, size_t qual_threshold, size_t window_size, bool stranded, bool no_left, bool no_right, bool no_orphans) {
    
    while(reader.has_next()) {
        auto i = reader.next();
        PairedEndRead* per = dynamic_cast<PairedEndRead*>(i.get());        
        if (per) {
            counters.input(*per);
            if (!no_left) {
                trim_left(per->non_const_read_one(), qual_threshold, window_size);
                trim_left(per->non_const_read_two(), qual_threshold, window_size);            
            }
            if (!no_right) {
                trim_right(per->non_const_read_one(), qual_threshold, window_size);
                trim_right(per->non_const_read_two(), qual_threshold, window_size);
            }
            per->checkDiscarded(min_length);
            writer_helper(per, pe, se, stranded, no_orphans);
            counters.output(*per);
        } else {
            SingleEndRead* ser = dynamic_cast<SingleEndRead*>(i.get());
            
            if (ser) {
                counters.input(*ser);
                if (!no_left) {
                    trim_left(ser->non_const_read_one(), qual_threshold, window_size);            
                } 
                if (!no_right) {
                    trim_right(ser->non_const_read_one(), qual_threshold, window_size);            
                }
                ser->checkDiscarded(min_length);
                writer_helper(ser, pe, se, stranded);
                counters.output(*ser);
            } else {
                throw std::runtime_error("Unknow read type");
            }
        }
    }

}

#endif
