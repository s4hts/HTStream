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

void trim_left(Read &rb, size_t min_qual) {

    std::string qual = rb.get_qual();

    int i;
    for (i = 0; i < qual.length(); ++i) {
        if (qual[i] > min_qual) {
            break;
        }
    }

    rb.setLCut(i );

}

void trim_right(Read &rb, size_t min_qual) {

    std::string qual = rb.get_qual();
    int len = qual.length() - 1;
    int i;
    for (i = len  ; i >=0; --i) {
        if (qual[i] > min_qual) {
            break;
        }
    }
    
    rb.setRCut(i + 1);

}

template <class T, class Impl>
void helper_trim(InputReader<T, Impl> &reader, std::shared_ptr<OutputWriter> pe, std::shared_ptr<OutputWriter> se, Counter& counters, size_t min_length, size_t min_qual, bool stranded, bool no_left, bool no_right) {
    
    while(reader.has_next()) {
        auto i = reader.next();
        ++counters["TotalRecords"];
        PairedEndRead* per = dynamic_cast<PairedEndRead*>(i.get());        
        if (per) {
            if (!no_left) {
                trim_left(per->non_const_read_one(), min_qual);            
                trim_left(per->non_const_read_two(), min_qual);            
            }
            if (!no_right) {
                trim_right(per->non_const_read_one(), min_qual);
                trim_right(per->non_const_read_two(), min_qual);            
            }
            per->checkDiscarded(min_length);
            writer_helper(per, pe, se, stranded);
        } else {
            SingleEndRead* ser = dynamic_cast<SingleEndRead*>(i.get());
            
            if (ser) {
                if (!no_left) {
                    trim_left(ser->non_const_read_one(), min_qual);            
                } 
                if (!no_right) {
                    trim_left(ser->non_const_read_one(), min_qual);            
                }
                ser->checkDiscarded(min_length);
                writer_helper(ser, pe, se, stranded);
            } else {
                throw std::runtime_error("Unknow read type");
            }
        }
    }

}

#endif
