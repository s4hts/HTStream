#ifndef SUPERD_H
#define SUPERD_H
//  this is so we can implment hash function for dynamic_bitset
#define BOOST_DYNAMIC_BITSET_DONT_USE_FRIENDS

#include "utils.h"
#include "ioHandler.h"
#include "main_template.h"

#include <map>
#include <unordered_map>
#include <boost/dynamic_bitset.hpp>
#include <boost/functional/hash.hpp>

extern template class InputReader<SingleEndRead, SingleEndReadFastqImpl>;
extern template class InputReader<PairedEndRead, PairedEndReadFastqImpl>;
extern template class InputReader<PairedEndRead, InterReadImpl>;
extern template class InputReader<ReadBase, TabReadImpl>;

class SuperDeduperCounters : public Counters {

public:
    std::vector<Vector> duplicateProportion;

    uint64_t Ignored = 0;
    uint64_t Duplicate = 0;

    SuperDeduperCounters(const std::string &program_name, const po::variables_map &vm) : Counters::Counters(program_name, vm) {
        generic.push_back(std::forward_as_tuple("ignored", Ignored));
        generic.push_back(std::forward_as_tuple("duplicate", Duplicate));
    }

    using Counters::input;
    virtual void input(const ReadBase &read, size_t dup_freq) {
        if (dup_freq > 0 && TotalFragmentsInput % dup_freq == 0 && TotalFragmentsInput != 0){
            duplicateProportion.push_back(std::forward_as_tuple(TotalFragmentsInput, Duplicate));
        }
        Counters::input(read);
    }

    void increment_replace() {
        ++Duplicate;
    }

    void increment_ignored() {
        ++Ignored;
    }

    virtual void write_out() {

        // record final input/dup
        duplicateProportion.push_back(std::forward_as_tuple(TotalFragmentsInput, Duplicate));

        initialize_json();

        write_labels(generic);
        write_vector("duplicate_saturation",duplicateProportion);
        write_sublabels("Single_end", se);
        write_sublabels("Paired_end", pe);

        finalize_json();
    }
};

class dbhash {
public:
    std::size_t operator() (const boost::dynamic_bitset<>& bs) const {
        return boost::hash_value(bs.m_bits);
    }
};

typedef std::unordered_map <boost::dynamic_bitset<>, std::unique_ptr<ReadBase>, dbhash> BitMap;

class SuperDeduper: public MainTemplate<SuperDeduperCounters, SuperDeduper> {
public:

    BitMap read_map;

    SuperDeduper() {
        program_name = "hts_SuperDeduper";
        app_description =
            "hts_SuperDeduper is a reference-free PCR duplicate remover. It uses a subsequence\n";
        app_description += "  within each read as a unique key to detect duplicates in future reads.\n";
        app_description += "  Reads with 'N' character(s) in the key sequence are ignored.\n";
        app_description += "  hts_SuperDeduper is not recommended for single-end reads.\n";
        app_description += "  WARNING: hts_SuperDeduper will only work correctly on untrimmed reads.\n";
    }

    void add_extra_options(po::options_description &desc) {
        desc.add_options()
            ("start,s", po::value<size_t>()->default_value(10)->notifier(boost::bind(&check_range<size_t>, "start", _1, 1, 10000)),  "Start location for unique ID (min 1, max 10000)")
            ("length,l", po::value<size_t>()->default_value(10)->notifier(boost::bind(&check_range<size_t>, "length", _1, 1, 10000)), "Length of unique ID (min 1, max 10000)")
            ("avg-qual-score,q", po::value<double>()->default_value(30)->notifier(boost::bind(&check_range<double>, "avg-qual-score", _1, 1, 10000)), "Avg quality score to have the read written automatically (min 1, max 10000)")
            ("inform-avg-qual-score,a", po::value<double>()->default_value(5)->notifier(boost::bind(&check_range<double>, "inform-avg-qual-score", _1, 1, 10000)), "Avg quality score to consider a read informative (min 1, max 10000)") //I know this says user input is a int, but is actually a double
            ("log_freq,e", po::value<size_t>()->default_value(1000000)->notifier(boost::bind(&check_range<size_t>, "log_freq", _1, 0, 1000000000)), "Frequency in which to log duplicates in reads, can be used to create a saturation plot (0 turns off).");
    }

    template <class T, class Impl>
    void load_map(InputReader<T, Impl> &reader, SuperDeduperCounters& counters, std::shared_ptr<OutputWriter> pe, std::shared_ptr<OutputWriter> se, double avg_automatic_write, double discard_qual, size_t start, size_t length, size_t log_freq ){
        double tmpAvg;

        while(reader.has_next()) {
            auto i = reader.next();
            counters.input(*i, log_freq);
            tmpAvg = i->avg_q_score();
            //check for existance, store or compare quality and replace:
            if ( tmpAvg < discard_qual ){ // averge qual must be less than discard_qual, ignored
                counters.increment_ignored();
            } else if (auto key=i->get_key(start, length)) { // check for duplicate
                // find faster than count on some compilers, new key
                if(read_map.find(*key) == read_map.end()) { // first time the key is seen
                    if ( tmpAvg >= avg_automatic_write ) { // if its greater than avg_automatic_write then write out
                        writer_helper(i.get(), pe, se, false);
                        counters.output(*i);
                        read_map[*key] = nullptr;
                    } else {
                        read_map[*key] = std::move(i);
                    }
                } else if (read_map[*key] == nullptr) { //key already seen and written out, PCR dup
                    counters.increment_replace();
                } else if( tmpAvg > read_map[*key]->avg_q_score()){ // new read is 'better' than old, key not yet read out
                    if (tmpAvg >= avg_automatic_write) { // read qualifies, write out
                        writer_helper(i.get(), pe, se, false);
                        counters.output(*i);
                        read_map[*key] = nullptr;
                    } else if (tmpAvg > discard_qual) {
                        read_map[*key] = std::move(i);
                    }
                    counters.increment_replace();
                } else {
                    counters.increment_replace(); // new read has been seen but is not better then last seen read with same key
                }
            } else {  // key had N, no key obtained
                counters.increment_ignored();
            }
        }
    }

    template <class T, class Impl>
    void do_app(InputReader<T, Impl> &reader, std::shared_ptr<OutputWriter> pe, std::shared_ptr<OutputWriter> se, SuperDeduperCounters &counter, const po::variables_map &vm) {
        const double avg_automatic_write = vm["avg-qual-score"].as<double>();
        const double discard_qual = vm["inform-avg-qual-score"].as<double>();
        const size_t start = vm["start"].as<size_t>() - 1;
        const size_t length = vm["length"].as<size_t>();
        const size_t log_freq = vm["log_freq"].as<size_t>();

        load_map(reader, counter, pe, se, avg_automatic_write, discard_qual, start, length, log_freq);
        for(auto const &i : read_map) {
            if (i.second.get() != nullptr) {
                counter.output(*i.second.get());
                writer_helper(i.second.get(), pe, se, false, false);
            }
        }
    }
};
#endif
