#ifndef PRIMERS_H
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


class PrimerCounters : public Counters {

public:

    PrimerCounters(const std::string &statsFile, bool force, bool appendStats, const std::string &program_name, const std::string &notes) : Counters::Counters(statsFile, force, appendStats, program_name, notes) {
    }

    using Counters::input;
    virtual void input(const ReadBase &read) {
        Counters::input(read);
    }

    virtual void write_out() {

        initialize_json();

        write_labels(generic);
        write_sublabels("Single_end", se);
        write_sublabels("Paired_end", pe);

        finalize_json();
    }
};


/* This is the helper class for Primer
 * The idea is ...
 * */
template <class T, class Impl>
void helper_Primers(InputReader<T, Impl> &reader, std::shared_ptr<OutputWriter> pe, std::shared_ptr<OutputWriter> se, PrimerCounters &counter, po::variables_map vm) {

    while(reader.has_next()) {
        auto i = reader.next();
        PairedEndRead* per = dynamic_cast<PairedEndRead*>(i.get());
        if (per) {
            counter.input(*per);
            counter.output(*per);
            writer_helper(per, pe, se);
        } else {
            SingleEndRead* ser = dynamic_cast<SingleEndRead*>(i.get());
            if (ser) {
                counter.input(*ser);
                counter.output(*ser);
                writer_helper(ser, pe, se);
            } else {
                throw std::runtime_error("Unknown read type");
            }
        }
    }
}

#endif
