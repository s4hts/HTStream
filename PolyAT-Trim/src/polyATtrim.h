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

void trim_left(Read &rb, size_t min_trim, size_t max_mismatch) {
    size_t a_mismatch = 0;
    size_t t_mismatch = 0;

    std::string seq = rb.get_seq();
    size_t temp_loc = 0;
    for (size_t i = 0; i < seq.length(); ++i) {
        if (seq[i] != 'A') {
            ++a_mismatch;
        } else if (a_mismatch <= max_mismatch) {
            temp_loc = i;
        }

        if (seq[i] != 'T') {
            ++t_mismatch;
        } else if (t_mismatch <= max_mismatch) {
            temp_loc = i;
        }

        if (t_mismatch > max_mismatch && a_mismatch > max_mismatch) {
            break;    
        }
    }

    if (temp_loc + 1 >= min_trim) {
        rb.setLCut(temp_loc + 1);
    }

}

void trim_right(Read &rb, size_t min_trim, size_t max_mismatch) {
    size_t a_mismatch = 0;
    size_t t_mismatch = 0;

    std::string seq = rb.get_seq();

    size_t len = seq.length() - 1;
    size_t temp_loc = len ;

    for (size_t i = len  ; i > 0; --i) {
        if (seq[i] != 'A') {
            ++a_mismatch;
        } else if (a_mismatch <= max_mismatch) {
            temp_loc = i;
        }

        if (seq[i] != 'T') {
            ++t_mismatch;
        } else if (t_mismatch <= max_mismatch) {
            temp_loc = i;
        }
        if (t_mismatch > max_mismatch && a_mismatch > max_mismatch) {
            break; 
        }
    }
    if ( len - temp_loc  >= min_trim) {
        rb.setRCut(temp_loc );
    }

}

template <class T, class Impl>
void helper_trim(InputReader<T, Impl> &reader, std::shared_ptr<OutputWriter> pe, std::shared_ptr<OutputWriter> se, Counter& counters, size_t min_length, size_t min_trim, size_t max_mismatch, bool stranded, bool no_left, bool no_right, bool no_orphans) {
    
    while(reader.has_next()) {
        auto i = reader.next();
        ++counters["TotalRecords"];
        PairedEndRead* per = dynamic_cast<PairedEndRead*>(i.get());        
        if (per) {
            if (!no_left) {
                trim_left(per->non_const_read_one(), min_trim, max_mismatch);
                trim_left(per->non_const_read_two(), min_trim, max_mismatch);            
            }
            if (!no_right) {
                trim_right(per->non_const_read_one(), min_trim, max_mismatch);
                trim_right(per->non_const_read_two(), min_trim, max_mismatch);            
            }
            per->checkDiscarded(min_length);
            writer_helper(per, pe, se, stranded, counters, no_orphans);
        } else {
            SingleEndRead* ser = dynamic_cast<SingleEndRead*>(i.get());
            
            if (ser) {
                if (!no_left) {
                    trim_left(ser->non_const_read_one(), min_trim, max_mismatch);            
                } 
                if (!no_right) {
                    trim_left(ser->non_const_read_one(), min_trim, max_mismatch);            
                }
                ser->checkDiscarded(min_length);
                writer_helper(ser, pe, se, stranded, counters);
            } else {
                throw std::runtime_error("Unknow read type");
            }
        }
    }

}

#endif
