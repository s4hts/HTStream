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

bool trim_left(Read &rb, char ccheck, size_t min_trim, size_t max_mismatch) {
    size_t mismatch = 0;

    std::string seq = rb.get_seq();
    std::string::iterator current_loc, tmp_loc = seq.begin();

    for (current_loc = seq.begin(); current_loc != seq.end() ; ++current_loc ) {
      if ( *current_loc != ccheck ) {
          ++mismatch;
      } else if (mismatch <= max_mismatch ) {
          tmp_loc = current_loc;
      }
      if (mismatch >= max_mismatch) {
          break;
      }
    }
    if (tmp_loc - seq.begin() + 1 >= static_cast<long> (min_trim) && seq.end() > tmp_loc ) {
        rb.setLCut( static_cast<size_t>( (tmp_loc) - (seq.begin()) ) + 1 );
        return true;
    }
    return false;
}

bool trim_right(Read &rb, char ccheck, size_t min_trim, size_t max_mismatch) {
    size_t mismatch = 0;

    std::string seq = rb.get_seq();
    std::string::reverse_iterator current_loc, tmp_loc = seq.rbegin();

    for (current_loc = seq.rbegin(); current_loc != seq.rend() ; ++current_loc ) {
        if ( *current_loc != ccheck ) {
            ++mismatch;
        } else if (mismatch <= max_mismatch ) {
            tmp_loc = current_loc;
        }
        if (mismatch >= max_mismatch) {
            break;
        }
    }
    if (tmp_loc - seq.rbegin() + 1 >= static_cast<long> (min_trim) && seq.rend() > tmp_loc ) {
        rb.setRCut( static_cast<size_t>( seq.rend() - (tmp_loc + 1)) );
        return true;
    }
    return false;
}

template <class T, class Impl>
void helper_trim(InputReader<T, Impl> &reader, std::shared_ptr<OutputWriter> pe, std::shared_ptr<OutputWriter> se, TrimmingCounters& counters, size_t min_length, bool no_r1, bool no_r2, bool no_pA, bool no_pT, size_t min_trim, size_t max_mismatch, bool stranded, bool no_left, bool no_right, bool no_orphans) {

    while(reader.has_next()) {
        auto i = reader.next();
        PairedEndRead* per = dynamic_cast<PairedEndRead*>(i.get());
        if (per) {
            counters.input(*per);
            if (!no_left) {
                if (!no_r1){
                  if (!no_pA) trim_left(per->non_const_read_one(), 'A', min_trim, max_mismatch);
                  if (!no_pT) trim_left(per->non_const_read_one(), 'T', min_trim, max_mismatch);
                }
                if (!no_r2){
                  if (!no_pA) trim_left(per->non_const_read_two(), 'A', min_trim, max_mismatch);
                  if (!no_pT) trim_left(per->non_const_read_two(), 'T', min_trim, max_mismatch);
                }
            }
            if (!no_right) {
                if (!no_r1){
                  if (!no_pA) trim_right(per->non_const_read_one(), 'A', min_trim, max_mismatch);
                  if (!no_pT) trim_right(per->non_const_read_one(), 'T', min_trim, max_mismatch);
                }
                if (!no_r2){
                  if (!no_pA) trim_right(per->non_const_read_two(), 'A', min_trim, max_mismatch);
                  if (!no_pT) trim_right(per->non_const_read_two(), 'T', min_trim, max_mismatch);
                }
            }
            per->checkDiscarded(min_length);
            writer_helper(per, pe, se, stranded, no_orphans);
            counters.output(*per, no_orphans);
        } else {
            SingleEndRead* ser = dynamic_cast<SingleEndRead*>(i.get());

            if (ser) {
                counters.input(*ser);
                if (!no_left) {
                  if (!no_r1){
                    if (!no_pA) trim_left(per->non_const_read_one(), 'A', min_trim, max_mismatch);
                    if (!no_pT) trim_left(per->non_const_read_one(), 'T', min_trim, max_mismatch);
                  }
                }
                if (!no_right) {
                  if (!no_r1){
                    if (!no_pA) trim_left(per->non_const_read_one(), 'A', min_trim, max_mismatch);
                    if (!no_pT) trim_left(per->non_const_read_one(), 'T', min_trim, max_mismatch);
                  }
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
