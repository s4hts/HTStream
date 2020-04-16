#ifndef CUT_TRIM_H
#define CUT_TRIM_H
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

class CutTrimCounters : public TrimmingCounters {

public:

    uint64_t SE_Discarded = 0;
    uint64_t PE_Discarded = 0;
    uint64_t R1_Discarded = 0;
    uint64_t R2_Discarded = 0;

    CutTrimCounters(const std::string &program_name, const po::variables_map &vm) : TrimmingCounters::TrimmingCounters(program_name, vm) {
        se.push_back(std::forward_as_tuple("discarded", SE_Discarded));
        pe.push_back(std::forward_as_tuple("discarded", PE_Discarded));
        r1.push_back(std::forward_as_tuple("discarded", R1_Discarded));
        r2.push_back(std::forward_as_tuple("discarded", R2_Discarded));
    }
    using TrimmingCounters::SE_stats;
    using TrimmingCounters::R1_stats;
    using TrimmingCounters::R2_stats;

    void output(SingleEndRead &ser)  {
        Read &one = ser.non_const_read_one();
        if (!one.getDiscard()) {
            ++TotalFragmentsOutput;
            ++SE_Out;
            TrimmingCounters::SE_stats(one);
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
            TrimmingCounters::R1_stats(one);
            TrimmingCounters::R2_stats(two);
        } else if (!one.getDiscard() && !no_orphans) {
            ++TotalFragmentsOutput;
            ++SE_Out;
            ++R2_Discarded;
            TrimmingCounters::SE_stats(one);
        } else if (!two.getDiscard() && !no_orphans) {
            ++TotalFragmentsOutput;
            ++SE_Out;
            ++R1_Discarded;
            TrimmingCounters::SE_stats(two);
        } else {
            ++PE_Discarded;
        }
    }
};

class CutTrim: public MainTemplate<CutTrimCounters, CutTrim> {
public:

    CutTrim() {
        program_name = "hts_CutTrim";
        app_description =
            "The hts_CutTrim application trims a set number of bases from the 5'\n";
        app_description += "  and/or 3' end of each read\n";
    }

    void add_extra_options(po::options_description &desc) {
        desc.add_options()
            ("r1-cut-left,a", po::value<size_t>()->default_value(0)->notifier(boost::bind(&check_range<size_t>, "r1-cut-left", _1, 0, 10000)), "Cut length of sequence from read 1 left (5') end (min 0, max 10000)");
        desc.add_options()
            ("r1-cut-right,b", po::value<size_t>()->default_value(0)->notifier(boost::bind(&check_range<size_t>, "r1-cut-right", _1, 0, 10000)), "Cut length of sequence from read 1 right (3') end (min 0, max 10000)");
        desc.add_options()
            ("r2-cut-left,c", po::value<size_t>()->default_value(0)->notifier(boost::bind(&check_range<size_t>, "r2-cut-left", _1, 0, 10000)), "Cut length of sequence from read 2 left (5') end (min 0, max 10000)");
        desc.add_options()
            ("r2-cut-right,d", po::value<size_t>()->default_value(0)->notifier(boost::bind(&check_range<size_t>, "r2-cut-right", _1, 0, 10000)), "Cut length of sequence from read 2 right (3') end (min 0, max 10000)");
        desc.add_options()
            ("min-length,m", po::value<size_t>()->default_value(0)->notifier(boost::bind(&check_range<size_t>, "min-length", _1, 0, 10000)), "Min length for acceptable output read (min 1, max 10000), default is unset");
        desc.add_options()
            ("max-length,M", po::value<size_t>()->default_value(0)->notifier(boost::bind(&check_range<size_t>, "max-length", _1, 0, 10000)), "Maximum allowed length of read, effectively right trims to max-length (min 1, max 10000), default is unset");
        desc.add_options()
            ("no-orphans,n", po::bool_switch()->default_value(false), "Orphaned SE reads will NOT be written out");
        desc.add_options()
            ("stranded,s", po::bool_switch()->default_value(false),    "If R1 is orphaned, R2 is output in RC (for stranded RNA)");
    }


/*
  cut_trim expected behavior:
  cut_trim is expected to be one of the first apps used in a preprocessing pipeline
  such as when you want to statically trim off X bases from the end of a read (cut_right), due to
  quality, and/or to trim off primer sequence from the beginning of a read (cut_left). As such both
  cut_left and cut_right should act on the original length of the read.

  max_length then prevents the final cut length from being greater than max_length applying an
  additional cut to end of the read if applicable. If one wished to simulate smaller sized reads
  can run cut_trim with max_length (no cut to left or right) to effectively reduce the size of reads.
*/
    void cut_trim(Read &r, size_t cut_left, size_t cut_right, size_t min_length, size_t max_length) {
        if (cut_left) {
            r.setLCut(cut_left);
        }
        if (cut_right) {
            r.setRCut(r.getLength() - cut_right);
        }
        if (max_length && max_length < r.getLengthTrue()) {
            std::cout << "max_length:" << max_length << " length:" << r.getLengthTrue() << "\n";
            r.setDiscard();
        }
        if (min_length && min_length > r.getLengthTrue()) {
            std::cout << "min_length:" << min_length << " length:" << r.getLengthTrue() << "\n";
            r.setDiscard();
        }
    }


    template <class T, class Impl>
    void do_app(InputReader<T, Impl> &reader, std::shared_ptr<OutputWriter> pe, std::shared_ptr<OutputWriter> se, CutTrimCounters& counters, const po::variables_map &vm) {

        size_t min_length  = vm["min-length"].as<size_t>();
        bool stranded =  vm["stranded"].as<bool>();
        bool no_orphans = vm["no-orphans"].as<bool>();
        size_t r1_cut_left = vm["r1-cut-left"].as<size_t>();
        size_t r1_cut_right = vm["r1-cut-right"].as<size_t>();
        size_t r2_cut_left = vm["r2-cut-left"].as<size_t>();
        size_t r2_cut_right = vm["r2-cut-right"].as<size_t>();
        size_t max_length = vm["max-length"].as<size_t>();

        while(reader.has_next()) {
            auto i = reader.next();
            PairedEndRead* per = dynamic_cast<PairedEndRead*>(i.get());
            if (per) {
                counters.input(*per);
                cut_trim( per->non_const_read_one(), r1_cut_left, r1_cut_right, min_length, max_length);
                cut_trim( per->non_const_read_two(), r2_cut_left, r2_cut_right, min_length, max_length);
                writer_helper(per, pe, se, stranded, no_orphans);
                counters.output(*per, no_orphans);
            } else {
                SingleEndRead* ser = dynamic_cast<SingleEndRead*>(i.get());
                if (ser) {
                    counters.input(*ser);
                    cut_trim( ser->non_const_read_one(), r1_cut_left, r1_cut_right, min_length, max_length);
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
