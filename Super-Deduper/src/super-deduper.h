//  this is so we can implment hash function for dynamic_bitset
#ifndef SUPERD_H
#define SUPERD_H
#define BOOST_DYNAMIC_BITSET_DONT_USE_FRIENDS

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

typedef std::unordered_map<std::string, size_t> Counter;
typedef std::unordered_map <boost::dynamic_bitset<>, std::unique_ptr<ReadBase>, dbhash> BitMap;

template <class T, class Impl>
void load_map(InputReader<T, Impl> &reader, Counter& counters, BitMap& read_map, size_t start, size_t length) {
    while(reader.has_next()) {
        auto i = reader.next();
        ++counters["TotalRecords"];
        //check for existance, store or compare quality and replace:
        if (auto key=i->get_key(start, length)) {
            // find faster than count on some compilers
            if(read_map.find(*key) == read_map.end()) {
                read_map[*key] = std::move(i);
            } else if(i->avg_q_score() > read_map[*key]->avg_q_score()){
                read_map[*key] = std::move(i);
                ++counters["Replaced"];
            }
        } else {  // key had N 
            ++counters["HasN"];
        }
    }
}

#endif
