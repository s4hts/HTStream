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
#include "main_template.h"

extern template class InputReader<SingleEndRead, SingleEndReadFastqImpl>;
extern template class InputReader<PairedEndRead, PairedEndReadFastqImpl>;
extern template class InputReader<PairedEndRead, InterReadImpl>;
extern template class InputReader<ReadBase, TabReadImpl>;

class PolyATTrim: public MainTemplate<TrimmingCounters, PolyATTrim> {
public:

    void add_extra_options(po::options_description &desc) {
        setDefaultParamsCutting(desc);
            // no-orphans|n ; stranded|s ; min-length|m
        setDefaultParamsTrim(desc);
            // no-left|l ; no-right|r

        desc.add_options()
            ("skip_polyA,j", po::bool_switch()->default_value(false), "Skip check for polyA sequence")
            ("skip_polyT,k", po::bool_switch()->default_value(false), "Skip check for polyT sequence")
            ("window-size,w", po::value<size_t>()->default_value(6)->notifier(boost::bind(&check_range<size_t>, "window-size", _1, 1, 10000)),    "Window size in which to trim (min 1, max 10000)")
            ("max-mismatch-errorDensity,e", po::value<double>()->default_value(0.30)->notifier(boost::bind(&check_range<double>, "max-mismatch-errorDensity", _1, 0.0, 1.0)), "Max percent of mismatches allowed in overlapped section (min 0.0, max 1.0)")
            ("perfect-windows,c", po::value<size_t>()->default_value(1)->notifier(boost::bind(&check_range<size_t>, "perfect-windows", _1, 0, 10000)),    "Number perfect match windows needed before a match is reported  (min 1, max 10000)")
            ("min-trim,M", po::value<size_t>()->default_value(5)->notifier(boost::bind(&check_range<size_t>, "min-trim", _1, 1, 10000)), "Min base pairs trim for AT tail (min 1, max 10000)")
            ("max-trim,x", po::value<size_t>()->default_value(30)->notifier(boost::bind(&check_range<size_t>, "max-trim", _1, 0, 10000)), "Max size a polyAT can be (min 0, max 10000)");
    }

    PolyATTrim() {
        program_name = "hts_PolyATTrim";
        app_description =
            "hts_PolyATTrim trims poly A and T sequences from a read.\n";
        app_description += "  The algorithm is borrowed from Fig 2, Bonfert et al. doi: 2017 10.1371/journal.pone.0170914\n";
        app_description += "  A sliding window of <window-size> (=6) is shifted from either end of the read\n";
        app_description += "  (adjustable with --no-left and --no-right) until the <max-mismatch-errorDensity> is\n";
        app_description += "  exceeded. The read is then trimmed as long as the following criteria are met:\n";
        app_description += "  \ta) at least <perfect-windows> (=1) were observed\n";
        app_description += "  \tb) at least <min-trim> (=5) bp will be trimmed\n";
        app_description += "  \tc) no more than <max-trim> (=30) bp will be trimmed\n";
        app_description += "  These settings may need to be adjusted depending on library type.";
    }

    bool trim_left(Read &rb, char ccheck, size_t min_trim, size_t max_trim, size_t window_size, double max_mismatch_errorDensity, size_t perfect_windows){
        size_t i = 0;
        size_t mismatch = 0, pwindows = 0;
        double max_mismatch = static_cast<double> (window_size) * max_mismatch_errorDensity;

        std::string seq = rb.get_seq();
        std::string::iterator current_loc, tmp_loc, trim_loc = seq.begin();
/*
  for ( trim_loc = seq.begin(), i = 0; i < window_size ; ++trim_loc, ++i ){
  if ( *trim_loc != ccheck ) {
  break;
  }
  }
*/
        for ( current_loc = seq.begin() + (window_size-1); current_loc < seq.end() ; ++current_loc ) {
            mismatch=0;
            for (tmp_loc = current_loc, i=window_size; i > 0; --tmp_loc, --i) {
                if ( *tmp_loc != ccheck ) {
                    ++mismatch;
                }
            }
            if (static_cast<double> (mismatch) > max_mismatch){
                break;
            } else if (mismatch == 0) {
                pwindows++;
            }
            trim_loc = current_loc;
            while ( *trim_loc != ccheck ){
                trim_loc--;
            }
        }

        if (pwindows >= perfect_windows && trim_loc - seq.begin() + 1 >= static_cast<long> (min_trim) && trim_loc - seq.begin() + 1 <= static_cast<long> (max_trim) && seq.end() > trim_loc ) {
            rb.setLCut( static_cast<size_t>( (trim_loc) - (seq.begin()) ) + 1 );
            return true;
        }
        return false;
    }

    bool trim_right(Read &rb, char ccheck, size_t min_trim, size_t max_trim, size_t window_size, double max_mismatch_errorDensity, size_t perfect_windows){
        size_t i = 0, pwindows = 0;
        double mismatch = 0;
        double max_mismatch = static_cast<double> (window_size) * max_mismatch_errorDensity;

        std::string seq = rb.get_seq();
        std::string::reverse_iterator current_loc, tmp_loc, trim_loc = seq.rbegin();
/*
  for ( trim_loc = seq.rbegin(), i = 0; i < window_size ; ++trim_loc, ++i ){
  if ( *trim_loc != ccheck ) {
  break;
  }
  }
*/
        for ( current_loc = seq.rbegin() + (window_size-1); current_loc < seq.rend() ; ++current_loc ) {
            mismatch=0;
            for (tmp_loc = current_loc, i=window_size; i > 0; --tmp_loc, --i) {
                if ( *tmp_loc != ccheck ) {
                    mismatch++;
                }
            }
            if (static_cast<double> (mismatch) > max_mismatch){
                break;
            } else if (mismatch == 0) {
                pwindows++;
            }
            trim_loc = current_loc;
            while ( *trim_loc != ccheck ){
                trim_loc--;
            }
        }

        if (pwindows >= perfect_windows && trim_loc - seq.rbegin() + 1 >= static_cast<long> (min_trim) && trim_loc - seq.rbegin() + 1 <= static_cast<long> (max_trim)  && seq.rend() > trim_loc ) {
            rb.setRCut( static_cast<size_t>( seq.rend() - (trim_loc + 1)) );
            return true;
        }
        return false;
    }

/*
  For each read, a sliding window of length window_size is shifted along the left or right end of the read
  sequence and the fraction of A’s (or T’s depending on strandedness of sequencing) is calculated within
  each window. A minimum of perfect_windows mustd have 100 As/Ts and the remaining windows > 0.3 A’s and
  this is used as a candidate poly(A) site.
  https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5279776/
*/
    template <class T, class Impl>
    void do_app(InputReader<T, Impl> &reader, std::shared_ptr<OutputWriter> pe, std::shared_ptr<OutputWriter> se, TrimmingCounters& counters, const po::variables_map &vm) {

        size_t min_length = vm["min-length"].as<std::size_t>();
        bool no_pA = vm["skip_polyA"].as<bool>();
        bool no_pT = vm["skip_polyT"].as<bool>();
        size_t min_trim = vm["min-trim"].as<std::size_t>();
        size_t window_size = vm["window-size"].as<std::size_t>();
        size_t perfect_windows = vm["perfect-windows"].as<std::size_t>();
        size_t max_trim = vm["max-trim"].as<std::size_t>();
        double max_mismatch_errorDensity = vm["max-mismatch-errorDensity"].as<double>();
        bool stranded = vm["stranded"].as<bool>();
        bool no_left = vm["no-left"].as<bool>();
        bool no_right = vm["no-right"].as<bool>();
        bool no_orphans = vm["no-orphans"].as<bool>();

        bool r1fA = false, r1fT = false, r2fA = false, r1rA = false;
        while(reader.has_next()) {
            auto i = reader.next();
            PairedEndRead* per = dynamic_cast<PairedEndRead*>(i.get());
            if (per) {
                counters.input(*per);
                if (!no_left){ // 5' fragment end detection
                    if (!no_pA) r1fA = trim_left(per->non_const_read_one(), 'A', min_trim, max_trim, window_size, max_mismatch_errorDensity, perfect_windows);
                    if (!no_pT && !r1fA) r1fT = trim_left(per->non_const_read_one(), 'T', min_trim, max_trim, window_size, max_mismatch_errorDensity, perfect_windows);
                }
                if (!no_right && !r1fA && !r1fT){ // 3' fragement end detection, can't have 5' and 3' polyAT
                    if (!no_pA) r2fA = trim_left(per->non_const_read_two(), 'A', min_trim, max_trim, window_size, max_mismatch_errorDensity, perfect_windows);
                    if (!no_pT && !r2fA) trim_left(per->non_const_read_two(), 'T', min_trim, max_trim, window_size, max_mismatch_errorDensity, perfect_windows);
                }
                // polyAT must be on one of the two fragment ends (when PE), don't bother checking read ends, reads should be overlapped otherwise.
                per->checkDiscarded(min_length);
                writer_helper(per, pe, se, stranded, no_orphans);
                counters.output(*per, no_orphans);
            } else {
                SingleEndRead* ser = dynamic_cast<SingleEndRead*>(i.get());
                if (ser) {
                    counters.input(*ser);
                    if (!no_left) {
                        if (!no_pA) r1fA = trim_left(ser->non_const_read_one(), 'A', min_trim, max_trim, window_size, max_mismatch_errorDensity, perfect_windows);
                        if (!no_pT && !r1fA) r1fT = trim_left(ser->non_const_read_one(), 'T', min_trim, max_trim, window_size, max_mismatch_errorDensity, perfect_windows);
                    }
                    if (!no_right && !r1fA && !r1fT) {
                        if (!no_pA) r1rA = trim_right(ser->non_const_read_one(), 'A', min_trim, max_trim, window_size, max_mismatch_errorDensity, perfect_windows);
                        if (!no_pT && ! r1rA) trim_right(ser->non_const_read_one(), 'T', min_trim, max_trim, window_size, max_mismatch_errorDensity, perfect_windows);
                    }
                    ser->checkDiscarded(min_length);
                    writer_helper(ser, pe, se, false, false);
                    counters.output(*ser);
                } else {
                    throw std::runtime_error("Unknown read type");
                }
            }
        }

    }

};
#endif
