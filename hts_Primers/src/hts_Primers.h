edit_distance#ifndef PRIMERS_H
#define PRIMERS_H
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

/*
compute the Levenstein distance between a (primer) and b (sequence read)
pegged to within f bases of the 5' end, bounded by edit distance k and
requiring m perfect matches at end.
*/

unsigned int bounded_edit_distance(std:string primer, Read &read, int float_distance, int edit_distance, int end_matches){
  // float_distance is the pre-primer float value, edit_distance is max error and end_matches is end matches
  size_t x, i, j, l, lastdiag, olddiag, cmin, endmatchcount;
  size_t column[primerlen - m + 1];
  size_t distance = edit_distance+1;

  std::string seq = read.get_seq();
  if (primerlen > read.length()) { // primer should never be greater than the seq
      return 0;
  }
  if (primerlen == 0){
      return 0;
  }
  for (x = 0; x <= float_distance; x++) {
      for (i = 1; i <= primerlen - end_matches; i++)
          column[i] = i;
      for (i = 1; i <= primerlen - end_matches + edit_distance ; i++) { // outer loop is the read
          column[0] = i;
          cmin = i;
          for (j = 1, lastdiag = i-1; j <= primerlen - end_matches ; j++) { // inner loop is the primer
              olddiag = column[j];
              column[j] = MIN3(column[j] + 1, column[j-1] + 1, lastdiag + (primer[j-1] == seq[x+i-1] ? 0 : 1));
              lastdiag = olddiag;
              if (column[j] < cmin) cmin = column[j];
          }
          if (cmin > edit_distance) break; // if the smallest value in the column is > max error break
          if (column[primerlen - end_matches] <= distance ){
              endmatchcount=0;
              for (l = 1; l <= end_matches; l++){
                  if (primer[primerlen - l] != seq[x + i + end_matches - l]){
                      break;
                  }
                  endmatchcount++;
              }
              if (endmatchcount == end_matches){
                  read.add_comment("PF:D", read.get_sub_seq(0,x));
                  read.add_comment("PM:D", primerlen - end_matches);
                  read.setLCut(x + i + end_matches);
                  distance = column[primerlen - end_matches ]; // bottom right node, global alignment
              }
          }
      }
  }
  return distance;
}

static PyObject *
bounded_editdist_distance_list(PyObject *self, PyObject *args)
{
    Tuple r, s;
    Py_ssize_t i;
    char *seq;
    int seqlen, c = 0, f = 0, k = 0, m = 0; // c = current index of best, k = maxdist, m = finalmatch

    PyListObject* primer_list_o = NULL;

    if (!PyArg_ParseTuple(args, "O!s#iii", &PyList_Type, &primer_list_o, &seq, &seqlen, &f, &k, &m)){
        return NULL;
    }
    Py_ssize_t primer_list_o_length = PyList_GET_SIZE(primer_list_o);

    for (i = 0; i < primer_list_o_length; i++) {
        PyObject * primer_o = PyList_GET_ITEM(primer_list_o, i);
        if (PyString_Check(primer_o)) {
            r = bounded_edit_distance(PyString_AS_STRING(primer_o), (int)PyString_GET_SIZE(primer_o), seq, seqlen, f, k, m);
            if ( i == 0){
                c = 0;
                s = r;
            } else if ( r.dist < s.dist ){
                c = (int)i;
                s = r;
            } else if ( s.dist == 0){ // can't get better than a perfect match!
                break;
            }
        } else {
            PyErr_SetString(PyExc_TypeError,
                "first argument must be a list of strings");
            return NULL;
        }
    }

    return Py_BuildValue("iiii", c, s.dist, s.spos, s.epos);
}

size_t dist(size_t x, size_t y) {
    assert(x <= y);
    return y - x;
}


/*This is a O(N) time algorithm
 * it will search for the longest base pair segment that has
 * no N's within it*/
void trim_n(Read &rb) {

    std::string seq = rb.get_seq();
    size_t bestLeft = 0, currentLeft = 0, bestRight = 0;
    size_t i = 0;

    for (std::string::iterator it = seq.begin(); it != seq.end(); ++it) {
        i = static_cast<size_t>(it - seq.begin() );

        if (*it == 'N') {
            currentLeft = i + 1;
        } else if (dist(bestLeft, bestRight) < dist(currentLeft, i)) {
            bestRight = i;
            bestLeft = currentLeft;
        }
    }
    rb.setLCut(bestLeft);
    rb.setRCut(bestRight+1);
}

unsigned int checkIfAdapter(Read &r1, Read &adapter, size_t loc1, size_t loc2, const double misDensity, const size_t &mismatch, const size_t &minOverlap ) {
    size_t minLoc = std::min(loc1, loc2);
    int loc1_t = loc1 - minLoc;
    int loc2_t = loc2 - minLoc;
    int r1_len = r1.getLength();
    int adapter_len = adapter.getLength();

    size_t maxLoop = std::min(r1_len - loc1_t, adapter_len - loc2_t);
    size_t maxMis = std::min(mismatch, static_cast<size_t>(maxLoop * misDensity));

    const std::string &seq1 = r1.get_seq();
    const std::string &seq_adapter = adapter.get_seq();

    auto i1 = seq1.begin();
    std::advance(i1, loc1_t);
    auto i2 = seq_adapter.begin();
    std::advance(i2, loc2_t);
    if (maxLoop < minOverlap || !threshold_mismatches(i1, i2, maxLoop, maxMis) ) {
        // no overlap identified
        return 0;
    }
    // overlap exists thats meet maxMis criteria
    r1.setRCut(loc1_t);

    return 1;
}

void check_read_se(SingleEndRead &se ) {

    Read &r1 = se.non_const_read_one();
    /* if adapter is longer than sequence, set length to same as seq */
    if (adapter_seq.length() < r1.getLength()){
        adapter_seq = adapter_seq.substr(0, r1.getLength());
    }

    Read adapter = Read(adapter_seq, "", "");

   /* checkL needs to be as long as or longer than the shortest read */
    size_t checkL = std::min(adapter.getLength(), checkLengths);
    /* kmer needs to be as long as or longer than the shortest read */
    size_t kkmer = std::min(adapter.getLength(), kmer);
    /* Create a map with non-overlapping kmers*/
    seqLookup mOne = readOneMap(r1.get_seq(), kkmer, kmerOffset);

    /*Do a quick check if the shorter read kmer shows up in longer read (read 2)
     * If it does, then try the brute force approach*/
    for (size_t bp = 0; bp < (checkLengths - kmer); ++bp) {
        // check first checkLength kmers at the beginning of read2 (sins)
        auto test = mOne.equal_range(adapter_seq.substr(bp, kmer));
        for (auto it = test.first; it != test.second; ++it) {
            unsigned int overlapped = checkIfAdapter(r1, adapter, it->second, bp, misDensity, mismatch, minOver);
            if (overlapped) {
                return;
            }
        }
    }
}

/* This is the helper class for Primer
 * The idea is ...
 * */
template <class T, class Impl>
void helper_Primer(InputReader<T, Impl> &reader, std::shared_ptr<OutputWriter> pe, std::shared_ptr<OutputWriter> se, PrimerCounters &counter, po::variables_map vm) {

    while(reader.has_next()) {
        auto i = reader.next();
        PairedEndRead* per = dynamic_cast<PairedEndRead*>(i.get());
        if (per) {
            counter.input(*per);
            check_read_pe(*per, misDensity, mismatch, minOver, checkLengths, kmer, kmerOffset, noFixBases);
            per->checkDiscarded(min_length);
            counter.output(*per, no_orphan);
            writer_helper(per, pe, se, stranded, no_orphan);
        } else {
            SingleEndRead* ser = dynamic_cast<SingleEndRead*>(i.get());
            if (ser) {
                counter.input(*ser);
                check_read_se(*ser, misDensity, mismatch, minOver, checkLengths, kmer, kmerOffset, adapter);
                ser->checkDiscarded(min_length);
                counter.output(*ser);
                writer_helper(ser, pe, se);
            } else {
                throw std::runtime_error("Unknown read type");
            }
        }
    }
}

#endif
