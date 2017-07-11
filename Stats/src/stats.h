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


class StatsCounters : public Counters {

public:

    StatsCounters() {
        Common();
        c["A"] = 0;
        c["T"] = 0;
        c["C"] = 0;
        c["G"] = 0;
        c["N"] = 0;
        c["R1_Len"] = 0;
        c["R2_Len"] = 0;
        c["SE_Len"] = 0;
        c["R1_Q_Sum"] = 0;
        c["R2_Q_Sum"] = 0;
        c["SE_Q_Sum"] = 0;
    }
   
    void read_stats(Read &r) {
        for (auto bp : r.get_seq()) {
            switch (bp) {
                case 'A':
                    ++c["A"];
                    break;
                case 'T':
                    ++c["T"];
                    break;
                case 'C':
                    ++c["C"];
                    break;
                case 'G':
                    ++c["G"];
                    break;
                case 'N':
                    ++c["N"];
                    break;
                default:
                    throw std::runtime_error("Unknown bp in stats counter");
            }
        }
    }
 
    void q_stats(Read &r, unsigned long long int &val) {
    }
 
    void output(PairedEndRead &per) {
        Counters::output(per);
        Read &one = per.non_const_read_one();
        Read &two = per.non_const_read_two();
        read_stats(one);
        read_stats(two);
        for (auto q : one.get_qual()) {
            c["R1_Q_Sum"] += q - 33;
        }
        for (auto q : two.get_qual()) {
            c["R2_Q_Sum"] += q - 33;
        }
        c["R1_Len"] += one.getLength();
        c["R2_Len"] += two.getLength();
    }
    
    void output(SingleEndRead &ser) {
        Counters::output(ser);
        Read &one = ser.non_const_read_one();
        read_stats(one);
        for (auto q : one.get_qual()) {
            c["SE_Q_Sum"] += q - 33;
        }
        c["SE_Len"] += one.getLength();
    }
 
};




template <class T, class Impl>
void helper_trim(InputReader<T, Impl> &reader, std::shared_ptr<OutputWriter> pe, std::shared_ptr<OutputWriter> se, StatsCounters& counters) {
    while(reader.has_next()) {
        auto i = reader.next();
        PairedEndRead* per = dynamic_cast<PairedEndRead*>(i.get());        
        if (per) {
            counters.input(*per);
            counters.output(*per);
            writer_helper(per, pe, se, false);
        } else {
            SingleEndRead* ser = dynamic_cast<SingleEndRead*>(i.get());
            if (ser) {
                counters.input(*ser);
                counters.output(*ser);
                writer_helper(ser, pe, se, false);
            } else {
                throw std::runtime_error("Unknown read type");
            }
        }
    }

}

#endif
