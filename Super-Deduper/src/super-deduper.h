#ifndef SUPERD_H
#define SUPERD_H
//  this is so we can implment hash function for dynamic_bitset
#define BOOST_DYNAMIC_BITSET_DONT_USE_FRIENDS

#include "utils.h"
#include "ioHandler.h"
#include <map>
#include <unordered_map>
#include <boost/dynamic_bitset.hpp>
#include <boost/functional/hash.hpp>

class dbhash {
public:
    std::size_t operator() (const boost::dynamic_bitset<>& bs) const {
        return boost::hash_value(bs.m_bits);
    }
};

typedef std::unordered_map <boost::dynamic_bitset<>, std::unique_ptr<ReadBase>, dbhash> BitMap;


class SuperDeduperCounters : public Counters {
public:
    SuperDeduperCounters() {
        Common();
        c["Replaced"] = 0;
        c["Ignored"] = 0;
    }

    void increment_replace() {
        ++c["Replaced"];
    }

    void increment_ignored() {
        ++c["Ignored"];
    }
};

template <class T, class Impl>
void load_map(InputReader<T, Impl> &reader, SuperDeduperCounters& counters, BitMap& read_map, std::shared_ptr<OutputWriter> pe, std::shared_ptr<OutputWriter> se, double avg_automatic_write, double discard_qual, size_t start, size_t length) {
    double tmpAvg;

    while(reader.has_next()) {
        auto i = reader.next();
        counters.input(*i);
        //check for existance, store or compare quality and replace:
        if (auto key=i->get_key(start, length)) {
            // find faster than count on some compilers
            if(read_map.find(*key) == read_map.end()) {
                if ((tmpAvg = i->avg_q_score()) > avg_automatic_write) {
                    writer_helper(i.get(), pe, se, false);
                    counters.output(*i);
                    read_map[*key] = nullptr;
                } else if (tmpAvg > discard_qual) {
                    read_map[*key] = std::move(i);
                }
            } else if (read_map[*key] == nullptr) { //key had a q-score of 20 or higher (it was all ready written out)
                counters.increment_replace();
            } else if((tmpAvg = i->avg_q_score()) > read_map[*key]->avg_q_score()){
                if (tmpAvg > avg_automatic_write) {
                    writer_helper(i.get(), pe, se, false);
                    counters.output(*i);
                    read_map[*key] = nullptr;
                } else if (tmpAvg > discard_qual) {
                    read_map[*key] = std::move(i);
                }
                counters.increment_replace();
            }
        } else {  // key had N 
            counters.increment_ignored();
        }
    }
}

#endif
