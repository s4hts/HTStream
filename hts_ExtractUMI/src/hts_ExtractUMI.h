#ifndef EXTRACT_UMI_H
#define EXTRACT_UMI_H
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


class ExtractUMI: public MainTemplate<TrimmingCounters, ExtractUMI> {
public:

    std::vector<char> read_options{'F', 'R', 'B'}; // possible paramters for read options

    // bases some parameters and allows passing UMIs between PE reads
    struct UMI {
      std::string seq;
      std::string qual;
      size_t length;
      bool discard;
      size_t qual_threshold;
      size_t avg_qual_threshold;
      size_t qual_offset;
      bool homopolymer;
      bool discard_n;
    };


    ExtractUMI() {
        program_name = "hts_ExtractUMI";
        app_description =
            "The hts_ExtractUMI application trims a set number of bases from the 5'\n";
        app_description += "  end of a read and appends it to the read ID. Program is \n";
        app_description += "  meant to be an HTStream drop-in for umi_tools extract \n";
        app_description += "  with some different features that is compatible with \n";
        app_description += "  HTStream streaming pipelines.\n";
    }

    void add_extra_options(po::options_description &desc) {
        desc.add_options()
            ("read,r", po::value<char>()->default_value('F')->notifier(boost::bind(&check_values<char>, "read", _1, read_options)), "Read from which to extract the UMI (F = Forward, R = Reverse, B = Both), ignored if SE")
            ("umi-length,l", po::value<size_t>()->default_value(6)->notifier(boost::bind(&check_range<size_t>, "umi_length", _1, 1, 36)), "Total length of UMI to extract (1, 36)")
            ("qual-score,q", po::value<size_t>()->default_value(0)->notifier(boost::bind(&check_range<size_t>, "qual-score", _1, 0, 10000)), "Threshold for quality score for any base within a UMI (min 1, max 10000), read pairs are discarded, default is unset")
            ("avg-qual-score,Q", po::value<size_t>()->default_value(0)->notifier(boost::bind(&check_range<size_t>, "avg-qual-score", _1, 0, 10000)), "Threshold for quality score average of UMI (min 1, max 10000), read pairs are discarded, default is unset")
            ("homopolymer,p", po::bool_switch()->default_value(false), "Remove reads with homopolymer UMIs")
            ("discard-n,n", po::bool_switch()->default_value(false), "Remove reads with UMIs containing an N");
    }


    void quality_check(UMI &umi) {

        size_t tmp_qual_threshold; // just to be sure original parameters are not overwritten

        if (umi.qual_threshold != 0) {

            tmp_qual_threshold = umi.qual_threshold + umi.qual_offset;

            for (std::string::iterator it = umi.qual.begin() ; it != umi.qual.end() ; ++it) {
                if (static_cast<size_t>(*it) < tmp_qual_threshold) {
                    umi.discard = true;
                    break;
                }
            }

        } else if (umi.avg_qual_threshold != 0) {

            tmp_qual_threshold = umi.avg_qual_threshold + umi.qual_offset;
            size_t current_sum = 0;
            size_t final_qual_threshold = tmp_qual_threshold * umi.length;

            for (std::string::iterator it = umi.qual.begin() ; it != umi.qual.end() ; ++it) {
                current_sum += static_cast<size_t>(*it);
            }

            if (current_sum < final_qual_threshold) {
                umi.discard = true;
            }

        }

    }


    void homopolymer_check(UMI &umi) { 
        if (umi.seq.find_first_not_of(umi.seq[0]) == std::string::npos) {
            umi.discard = true;
        }
    }

    
    void n_check(UMI &umi) { 
        if (umi.seq.find('N') != std::string::npos) {
            umi.discard = true;
        }
    }


    void extract_umi(Read &r, UMI &umi) {

        if (umi.seq.empty()) {

            umi.seq = r.get_seq().substr(0, umi.length);

            if (umi.qual_threshold + umi.avg_qual_threshold != 0) {
                umi.qual = r.get_qual().substr(0, umi.length);
                quality_check(umi);
            }

            if (umi.homopolymer) {
                homopolymer_check(umi);
            }

            if (umi.discard_n) {
                n_check(umi);
            }

            if (!umi.discard) {
                r.setLCut(umi.length);
            }

        }

        if (!umi.discard) {
            r.set_id_first(r.get_id_first() + "_" + umi.seq);
        } else {
            r.setDiscard();
        }  

    }


    template <class T, class Impl>
    void do_app(InputReader<T, Impl> &reader, std::shared_ptr<OutputWriter> pe, std::shared_ptr<OutputWriter> se, TrimmingCounters& counters, const po::variables_map &vm) {
        
        char read = vm["read"].as<char>();
        size_t umi_length = vm["umi-length"].as<size_t>(); 
        size_t qual_offset = vm["qual-offset"].as<size_t>();
        size_t qual_threshold = vm["qual-score"].as<size_t>();
        size_t avg_qual_threshold = vm["avg-qual-score"].as<size_t>();
        bool homopolymer =  vm["homopolymer"].as<bool>();
        bool discard_n =  vm["discard-n"].as<bool>();


        // init UMI struct
        UMI umi = {
                   "",                     // umi sequence
                   "",                     // qual sequence
                   umi_length,             // umi length
                   false,                  // discard status (init)
                   qual_threshold,         // quality threshold
                   avg_qual_threshold,     // avg quality threshold
                   qual_offset,            // quality offset  
                   homopolymer,            // homopolymer filter
                   discard_n               // discard N containing UMIs
                  };


        WriterHelper writer(pe, se);

        auto read_visit = make_read_visitor_func(
            [&](SingleEndRead *ser) {
                extract_umi( ser->non_const_read_one(), umi);
            },
            [&](PairedEndRead *per) {
                if (read == 'F') {
                    extract_umi( per->non_const_read_one(), umi );
                    extract_umi( per->non_const_read_two(), umi );
                } else if (read == 'R') {
                    extract_umi( per->non_const_read_two(), umi );
                    extract_umi( per->non_const_read_one(), umi );
                } else {
                    extract_umi( per->non_const_read_one(), umi );
                    std::tie(umi.seq, umi.qual) = std::make_tuple("", ""); // reset umi struct
                    extract_umi( per->non_const_read_two(), umi );
                }
            }
            );

        while(reader.has_next()) {
            auto i = reader.next();
            counters.input(*i);
            i->accept(read_visit);
            std::tie(umi.seq, umi.qual, umi.discard) = std::make_tuple("", "", false); // reset umi struct
            writer(*i);
            counters.output(*i);
        }
    }
};
#endif
