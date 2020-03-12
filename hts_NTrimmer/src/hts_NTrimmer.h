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

class NTrimmer: public MainTemplate<TrimmingCounters, NTrimmer> {
public:

    NTrimmer() {
        program_name = "hts_NTrimmer";
        app_description =
            "The hts_NTrimmer application will identify and return the longest\n";
        app_description += "  subsequence that no N characters appear in.\n";
    }

    void add_extra_options(po::options_description &desc) {
        setDefaultParamsCutting(desc);
            // no-orphans|n ; stranded|s ; min-length|m

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
    void trim_n(Read &rb, bool exclude) {

        std::string seq = rb.get_seq();
        size_t bestLeft = 0, currentLeft = 0, bestRight = 0;
        size_t i = 0;

        for (std::string::iterator it = seq.begin(); it != seq.end(); ++it) {
            i = static_cast<size_t>(it - seq.begin() );

            if (*it == 'N') {
                if (exclude) {
                    rb.setLCut(1);
                    rb.setRCut(0);
                    return;
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
    }


/*Removes all Ns (ambiguity base) from a read*/
    template <class T, class Impl>
    void do_app(InputReader<T, Impl> &reader, std::shared_ptr<OutputWriter> pe, std::shared_ptr<OutputWriter> se, TrimmingCounters& counter, const po::variables_map &vm) {

        const size_t min_length = vm["min-length"].as<size_t>();
        bool stranded = vm["stranded"].as<bool>();
        bool no_orphans = vm["no-orphans"].as<bool>();
        bool exclude = vm["exclude"].as<bool>();
        
        while(reader.has_next()) {
            auto i = reader.next();
            PairedEndRead* per = dynamic_cast<PairedEndRead*>(i.get());
            if (per) {
                counter.input(*per);
                trim_n(per->non_const_read_one(), exclude);
                trim_n(per->non_const_read_two(), exclude);
                per->checkDiscarded(min_length);
                writer_helper(per, pe, se, stranded, no_orphans);
                counter.output(*per, no_orphans);
            } else {
                SingleEndRead* ser = dynamic_cast<SingleEndRead*>(i.get());
                if (ser) {
                    counter.input(*ser);
                    trim_n(ser->non_const_read_one(), exclude);
                    ser->checkDiscarded(min_length);
                    writer_helper(ser, pe, se, false, false);
                    counter.output(*ser);
                } else {
                    throw std::runtime_error("Unknow read type");
                }
            }
        }

    }
};
#endif
