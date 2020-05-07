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


    LengthFilterCounters(const std::string &program_name, const po::variables_map &vm) : Counters::Counters(program_name, vm) {
        se.push_back(std::forward_as_tuple("discarded", SE_Discarded));
        pe.push_back(std::forward_as_tuple("discarded", PE_Discarded));
        r1.push_back(std::forward_as_tuple("discarded", R1_Discarded));
        r2.push_back(std::forward_as_tuple("discarded", R2_Discarded));
    }

    void output(SingleEndRead &ser)  {
        Read &one = ser.non_const_read_one();
        if (!one.getDiscard()) {
            ++TotalFragmentsOutput;
            ++SE_Out;
            SE_BpLen_Out += one.getLength();
            TotalBasepairsOutput += one.getLength();
        } else {
            ++SE_Discarded;
        }
    }

    void output(PairedEndRead &per, bool no_orphans) {
        Read &one = per.non_const_read_one();
        Read &two = per.non_const_read_two();
        if (!one.getDiscard() && !two.getDiscard()) {
            ++TotalFragmentsOutput;
            ++PE_Out;
            R1_BpLen_Out += one.getLengthTrue();
            R2_BpLen_Out += two.getLengthTrue();
            TotalBasepairsOutput += one.getLengthTrue();
            TotalBasepairsOutput += two.getLengthTrue();
        } else if (!one.getDiscard() && !no_orphans) {
            ++TotalFragmentsOutput;
            ++SE_Out;
            SE_BpLen_Out += one.getLengthTrue();
            TotalBasepairsOutput += one.getLengthTrue();
            ++R2_Discarded;
        } else if (!two.getDiscard() && !no_orphans) {
            ++TotalFragmentsOutput;
            ++SE_Out;
            SE_BpLen_Out += two.getLengthTrue();
            TotalBasepairsOutput += two.getLengthTrue();
            ++R1_Discarded;
        } else {
            ++PE_Discarded;
        }
    }

private:
    using Counters::output;  // overload the base class and ignore warnings
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

        while(reader.has_next()) {
            ReadBasePtr i = std::static_pointer_cast<ReadBase>(std::shared_ptr<T>(reader.next()));
            PairedEndReadPtr per = std::dynamic_pointer_cast<PairedEndRead>(i);
            if (per) {
                counters.input(*per);
                length_filter(per->non_const_read_one(), min_length, max_length);
                length_filter(per->non_const_read_two(), min_length, max_length);
                writer_helper(per.get(), pe, se, stranded, no_orphans);
                counters.output(*per, no_orphans);
            } else {
                SingleEndReadPtr ser = std::dynamic_pointer_cast<SingleEndRead>(i);
                if (ser) {
                    counters.input(*ser);
                    length_filter(ser->non_const_read_one(), min_length, max_length);
                    writer_helper(ser.get(), se);
                    counters.output(*ser);
                } else {
                    throw std::runtime_error("Unknown read type");
                }
            }
        }
    }
};
#endif
