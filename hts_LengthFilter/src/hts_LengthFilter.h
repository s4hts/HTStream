#ifndef LENGTH_FILTER_H
#define LENGTH_FILTER_H
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

class LengthFilterCounters : public Counters {

public:

    uint64_t SE_Discarded = 0;
    uint64_t PE_Discarded = 0;
    uint64_t R1_Discarded = 0;
    uint64_t R2_Discarded = 0;
    bool no_orphans = false;

    virtual ~LengthFilterCounters() {}
    LengthFilterCounters(const std::string &program_name, const po::variables_map &vm) : Counters::Counters(program_name, vm) {
        se.push_back(std::forward_as_tuple("discarded", SE_Discarded));
        pe.push_back(std::forward_as_tuple("discarded", PE_Discarded));
        r1.push_back(std::forward_as_tuple("discarded", R1_Discarded));
        r2.push_back(std::forward_as_tuple("discarded", R2_Discarded));
    }
    using Counters::output;

    virtual void output(SingleEndRead &ser) {
        if (!ser.get_read().getDiscard()) {
            ++TotalFragmentsOutput;
            ++SE_Out;
            SE_BpLen_Out += ser.get_read().getLength();
            TotalBasepairsOutput += ser.get_read().getLength();
        } else {
            ++SE_Discarded;
        }
    }

    virtual void output(PairedEndRead &per) {
        if (!per.get_read_one().getDiscard() && !per.get_read_two().getDiscard()) {
            ++TotalFragmentsOutput;
            ++PE_Out;
            R1_BpLen_Out += per.get_read_one().getLengthTrue();
            R2_BpLen_Out += per.get_read_two().getLengthTrue();
            TotalBasepairsOutput += per.get_read_one().getLengthTrue();
            TotalBasepairsOutput += per.get_read_two().getLengthTrue();
        } else if (!per.get_read_one().getDiscard() && !no_orphans) {
            ++TotalFragmentsOutput;
            ++SE_Out;
            SE_BpLen_Out += per.get_read_one().getLengthTrue();
            TotalBasepairsOutput += per.get_read_one().getLengthTrue();

            ++R2_Discarded;
        } else if (!per.get_read_two().getDiscard() && !no_orphans) {
            ++TotalFragmentsOutput;
            ++SE_Out;
            SE_BpLen_Out += per.get_read_two().getLengthTrue();
            TotalBasepairsOutput += per.get_read_two().getLengthTrue();
            ++R1_Discarded;
        } else {
            ++PE_Discarded;
        }
    }
};

class LengthFilter: public MainTemplate<LengthFilterCounters, LengthFilter> {
public:

    LengthFilter() {
        program_name = "hts_LengthFilter";
        app_description =
            "The hts_LengthFilter application trims a set number of bases from the 5'\n";
        app_description += "  and/or 3' end of each read\n";
    }

    void add_extra_options(po::options_description &desc) {

        desc.add_options()
            ("min-length,m", po::value<size_t>()->default_value(0)->notifier(boost::bind(&check_range<size_t>, "min-length", _1, 0, 10000)), "Min length for acceptable output read (min 1, max 10000), default is unset");
        desc.add_options()
            ("max-length,M", po::value<size_t>()->default_value(0)->notifier(boost::bind(&check_range<size_t>, "max-length", _1, 0, 10000)), "Maximum allowed length of read, effectively right trims to max-length (min 1, max 10000), default is unset");
        desc.add_options()
            ("no-orphans,n", po::bool_switch()->default_value(false), "Orphaned SE reads will NOT be written out");
        desc.add_options()
            ("stranded,s", po::bool_switch()->default_value(false),    "If R1 is orphaned, R2 is output in RC (for stranded RNA)");
    }


/* TODO update the comments here
 */
    void length_filter(Read &r, size_t min_length, size_t max_length) {
        if (max_length && max_length < r.getLengthTrue()) {
            r.setDiscard();
        }
        if (min_length && min_length > r.getLengthTrue()) {
            r.setDiscard();
        }
    }


    template <class T, class Impl>
    void do_app(InputReader<T, Impl> &reader, std::shared_ptr<OutputWriter> pe, std::shared_ptr<OutputWriter> se, LengthFilterCounters& counters, const po::variables_map &vm) {

        size_t min_length  = vm["min-length"].as<size_t>();
        bool stranded =  vm["stranded"].as<bool>();
        bool no_orphans = vm["no-orphans"].as<bool>();
        size_t max_length = vm["max-length"].as<size_t>();
        counters.no_orphans = no_orphans;

        while(reader.has_next()) {
            auto i = reader.next();
            std::for_each(i->get_reads_non_const().begin(), i->get_reads_non_const().end(), ([=](const ReadPtr &read) { return length_filter(*read, min_length, max_length); }));
            writer_helper(i.get(), pe, se, stranded, no_orphans);
            counters.output(*i);
        }
    }
};
#endif
