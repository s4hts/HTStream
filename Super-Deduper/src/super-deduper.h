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


#define AVG_AUTOMATIC_WRITE 20
class dbhash {
public:
    std::size_t operator() (const boost::dynamic_bitset<>& bs) const {
        return boost::hash_value(bs.m_bits);
    }
};

typedef std::unordered_map <boost::dynamic_bitset<>, std::unique_ptr<ReadBase>, dbhash> BitMap;

template <class T, class Impl>
void load_map(InputReader<T, Impl> &reader, Counter& counters, BitMap& read_map, std::shared_ptr<OutputWriter> pe, std::shared_ptr<OutputWriter> se, size_t start, size_t length) {
    double tmpAvg;

    while(reader.has_next()) {
        auto i = reader.next();
        ++counters["TotalRecords"];
        //check for existance, store or compare quality and replace:
        if (auto key=i->get_key(start, length)) {
            // find faster than count on some compilers
            if(read_map.find(*key) == read_map.end()) {
                if (i->avg_q_score() > AVG_AUTOMATIC_WRITE) {
                    writer_helper(i.get(), pe, se, false, counters);
                    read_map[*key] = nullptr;
                } else {
                    read_map[*key] = std::move(i);
                }
            } else if (read_map[*key] == nullptr) { //key had a q-score of 20 or higher (it was all ready written out)
                ++counters["Replaced"];
            } else if((tmpAvg = i->avg_q_score()) > read_map[*key]->avg_q_score()){
                if (tmpAvg > AVG_AUTOMATIC_WRITE) {
                    writer_helper(i.get(), pe, se, false, counters);
                    read_map[*key] = nullptr;
                } else {
                    read_map[*key] = std::move(i);
                }
                ++counters["Replaced"];
            }
        } else {  // key had N 
            ++counters["Ignored"];
        }
    }
}

#endif
