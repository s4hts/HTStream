#ifndef N_TRIM_H
#define N_TRIM_H
//  this is so we can implment hash function for dynamic_bitset
#define BOOST_DYNAMIC_BITSET_DONT_USE_FRIENDS

#include "ioHandler.h"
#include "utils.h"
#include "main_template.h"

#include <map>
#include <unordered_map>
#include <boost/dynamic_bitset.hpp>
#include <boost/functional/hash.hpp>
#include <algorithm>

extern template class InputReader<SingleEndRead, SingleEndReadFastqImpl>;
extern template class InputReader<PairedEndRead, PairedEndReadFastqImpl>;
extern template class InputReader<PairedEndRead, InterReadImpl>;
extern template class InputReader<ReadBase, TabReadImpl>;

class NTrimCounters : public TrimmingCounters {

public:

    uint64_t SE_Discarded = 0;
    uint64_t PE_Discarded = 0;

    NTrimCounters(const std::string &program_name, const po::variables_map &vm) : TrimmingCounters::TrimmingCounters(program_name, vm) {
        se.push_back(std::forward_as_tuple("discarded", SE_Discarded));
        pe.push_back(std::forward_as_tuple("discarded", PE_Discarded));
    }

    void increment_discard_se(){
        SE_Discarded++;
    }

    void increment_discard_pe(){
        PE_Discarded++;
    }
};

class NTrimmer: public MainTemplate<NTrimCounters, NTrimmer> {
public:

    NTrimmer() {
        program_name = "hts_NTrimmer";
        app_description =
            "The hts_NTrimmer application will identify and return the longest\n";
        app_description += "  subsequence that no N characters appear in.\n";
    }

    void add_extra_options(po::options_description &desc) {
        desc.add_options()
            ("exclude,e", po::bool_switch()->default_value(false), "Exclude any sequence with an N character");
    }

    size_t dist(size_t x, size_t y) {
        assert(x <= y);
        return y - x;
    }

/*This is a O(N) time algorithm
 * it will search for the longest base pair segment that has
 * no N's within it*/
    bool trim_n(Read &rb, bool exclude) {

        std::string seq = rb.get_seq();
        size_t bestLeft = 0, currentLeft = 0, bestRight = 0;
        size_t i = 0;

        for (std::string::iterator it = seq.begin(); it != seq.end(); ++it) {
            i = static_cast<size_t>(it - seq.begin() );

            if (*it == 'N') {
                if (exclude) {
                    return false;
                } else {
                    currentLeft = i + 1;
                }
            } else if (dist(bestLeft, bestRight) < dist(currentLeft, i)) {
                bestRight = i;
                bestLeft = currentLeft;
            }
        }
        rb.setLCut(bestLeft);
        rb.setRCut(bestRight+1);
        return true;
    }


/*Removes all Ns (ambiguity base) from a read*/
    template <class T, class Impl>
    void do_app(InputReader<T, Impl> &reader, std::shared_ptr<OutputWriter> pe, std::shared_ptr<OutputWriter> se, NTrimCounters& counter, const po::variables_map &vm) {

        bool exclude = vm["exclude"].as<bool>();
        WriterHelper writer(pe, se);

        auto read_visit = make_read_visitor_func(
            [&](SingleEndRead *ser) {
                if (trim_n(ser->non_const_read_one(), exclude)){
                    writer(*ser);
                    counter.output(*ser);
                } else {
                    counter.increment_discard_se();
                }
            },
            [&](PairedEndRead *per) {
                if (trim_n(per->non_const_read_one(), exclude)){
                    if (trim_n(per->non_const_read_two(), exclude)){
                        writer(*per);
                        counter.output(*per);
                    } else {
                        counter.increment_discard_pe();
                    }
                } else {
                    counter.increment_discard_pe();
                }
            });

        while(reader.has_next()) {
            auto i = reader.next();
            counter.input(*i);
            i->accept(read_visit);
        }
    }
};
#endif
