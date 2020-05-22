#ifndef QWIN_TRIM_H
#define QWIN_TRIM_H
//  this is so we can implment hash function for dynamic_bitset
#define BOOST_DYNAMIC_BITSET_DONT_USE_FRIENDS

#include "ioHandler.h"
#include "utils.h"
#include "counters.h"
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

class QWindowTrim: public MainTemplate<TrimmingCounters, QWindowTrim> {
public:

    QWindowTrim() {
        program_name = "hts_QWindowTrim";
        app_description =
            "hts_QWindowTrim uses a sliding window approach to remove low quality\n";
        app_description += "  bases (5' or 3') from a read. A window will slide from each end of the\n";
        app_description += "  read, moving inwards. Once the window reaches an average quality <avg-qual-score>\n";
        app_description += "  it will stop trimming.";
    }

    void add_extra_options(po::options_description &desc) {
        setDefaultParamsTrim(desc);
        // no-left|l ; no-right|r

        desc.add_options()
            ("window-size,w", po::value<size_t>()->default_value(10)->notifier(boost::bind(&check_range<size_t>, "window-size", _1, 1, 10000)),    "Window size in which to trim (min 1, max 10000)")
            ("avg-qual-score,q", po::value<size_t>()->default_value(20)->notifier(boost::bind(&check_range<size_t>, "avg-qual-score", _1, 1, 10000)),    "Threshold for quality score average in the window (min 1, max 10000)")
            ("qual-offset,o", po::value<size_t>()->default_value(33)->notifier(boost::bind(&check_range<size_t>, "qual-offset", _1, 1, 10000)), "Quality offset for ascii q-score (default is 33) (min 1, max 10000)");
    }

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
    void do_app(InputReader<T, Impl> &reader, std::shared_ptr<OutputWriter> pe, std::shared_ptr<OutputWriter> se, TrimmingCounters &counters, const po::variables_map &vm) {

        size_t qual_threshold = vm["avg-qual-score"].as<size_t>() + vm["qual-offset"].as<size_t>() ;
        size_t window_size = vm["window-size"].as<size_t>();
        bool no_left = vm["no-left"].as<bool>();
        bool no_right = vm["no-right"].as<bool>();
        WriterHelper writer(pe, se);

        auto read_visit = make_read_visitor_func(
            [&](SingleEndRead *ser) {
                if (!no_left) {
                    trim_left(ser->non_const_read_one(), qual_threshold, window_size);
                }
                if (!no_right) {
                    trim_right(ser->non_const_read_one(), qual_threshold, window_size);
                }
            },
            [&](PairedEndRead *per) {
                if (!no_left) {
                    trim_left(per->non_const_read_one(), qual_threshold, window_size);
                    trim_left(per->non_const_read_two(), qual_threshold, window_size);
                }
                if (!no_right) {
                    trim_right(per->non_const_read_one(), qual_threshold, window_size);
                    trim_right(per->non_const_read_two(), qual_threshold, window_size);
                }
            });

        while(reader.has_next()) {
            auto i = reader.next();
            counters.input(*i);

            i->accept(read_visit);

            writer(*i);
            counters.output(*i);
        }
    }
};

#endif
