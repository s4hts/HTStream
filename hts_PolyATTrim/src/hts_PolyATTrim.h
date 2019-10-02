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

extern template class InputReader<SingleEndRead, SingleEndReadFastqImpl>;
extern template class InputReader<PairedEndRead, PairedEndReadFastqImpl>;
extern template class InputReader<PairedEndRead, InterReadImpl>;
extern template class InputReader<ReadBase, TabReadImpl>;

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
void helper_trim(InputReader<T, Impl> &reader, std::shared_ptr<OutputWriter> pe, std::shared_ptr<OutputWriter> se, TrimmingCounters& counters, size_t min_length, bool no_pA, bool no_pT, size_t min_trim, size_t window_size, size_t perfect_windows, size_t max_trim, double max_mismatch_errorDensity, bool stranded, bool no_left, bool no_right, bool no_orphans) {
    bool r1fA, r1fT, r2fA, r2fT, r1rA, r1rT, r2rA, r2rT = false;
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
              if (!no_pT && !r2fA) r2fT = trim_left(per->non_const_read_two(), 'T', min_trim, max_trim, window_size, max_mismatch_errorDensity, perfect_windows);
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
                  if (!no_pT && ! r1rA) r1rT = trim_right(ser->non_const_read_one(), 'T', min_trim, max_trim, window_size, max_mismatch_errorDensity, perfect_windows);
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

#endif
