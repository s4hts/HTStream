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

    std::vector<char> read_options{'F', 'R', 'B'};

    struct UMI {
      std::string seq = "";
      std::string qual = "";
      int length = 0;
      bool discard = false;
    } umi;

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
            ("read,r", po::value<char>()->default_value('F')->notifier(boost::bind(&check_char<char>, "read", _1, read_options)), "Read from which to extract the UMI (F = Forward, R = Reverse, B = Both), forward if SE")
            ("umi-length,l", po::value<size_t>()->default_value(6)->notifier(boost::bind(&check_range<size_t>, "umi_length", _1, 1, 36)), "Total length of UMI to extract (1, 36)")
            ("qual-score,q", po::value<size_t>()->default_value(0)->notifier(boost::bind(&check_range<size_t>, "qual-score", _1, 0, 10000)), "Threshold for quality score for any base within a UMI (min 1, max 10000), read pairs are discarded, default is unset")
            ("avg-qual-score,Q", po::value<size_t>()->default_value(0)->notifier(boost::bind(&check_range<size_t>, "avg-qual-score", _1, 0, 10000)), "Threshold for quality score average of UMI (min 1, max 10000), read pairs are discarded, default is unset")
            ("homopolymer,p", po::bool_switch()->default_value(false), "Remove reads with homopolymer UMIs")
            ("discard-n,n", po::bool_switch()->default_value(false), "Remove reads with UMIs containing an N");
    }

    void quality_check(size_t qual_threshold, size_t avg_qual_threshold, size_t qual_offset) {

        if (qual_threshold != 0) {

            qual_threshold += qual_offset;

            for (std::string::iterator it = umi.qual.begin() ; it != umi.qual.end() ; ++it) {
                if (static_cast<size_t>(*it) < qual_threshold) {
                    umi.discard = true;
                    break;
                }
            }

        } else if (avg_qual_threshold != 0) {

            avg_qual_threshold += qual_offset;
            size_t current_sum = 0;
            size_t final_qual_threshold = avg_qual_threshold * umi.length;

            for (std::string::iterator it = umi.qual.begin() ; it != umi.qual.end() ; ++it) {
                current_sum += static_cast<size_t>(*it);
            }

            if (current_sum < final_qual_threshold) {
                umi.discard = true;
            }

        }

    }


    void homopolymer_check() { 
        if (umi.seq.find_first_not_of(umi.seq[0]) == std::string::npos) {
            umi.discard = true;
        }
    }
    
    void n_check() { 
        if (umi.seq.find('N') != std::string::npos) {
            umi.discard = true;
        }
    }


    void extract_umi(Read &r, size_t umi_length, size_t qual_threshold, size_t avg_qual_threshold, 
                        size_t qual_offset, bool homopolymer, bool discard_n) {

        if (umi.seq.empty()) {

            umi.seq = r.get_seq().substr(0, umi_length);
            umi.length = umi_length;

            if (qual_threshold + avg_qual_threshold != 0) {
                umi.qual = r.get_qual().substr(0, umi_length);
                quality_check(qual_threshold, avg_qual_threshold, qual_offset);
            }

            if (homopolymer) {
                homopolymer_check();
            }

            if (discard_n) {
                n_check();
            }

            if (!umi.discard) {
                r.setLCut(umi_length);
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
        
        char read = vm["read"].as<size_t>();
        size_t umi_length = vm["umi-length"].as<size_t>(); 
        size_t qual_offset = vm["qual-offset"].as<size_t>();
        size_t qual_threshold = vm["qual-score"].as<size_t>();
        size_t avg_qual_threshold = vm["avg-qual-score"].as<size_t>();
        bool homopolymer =  vm["homopolymer"].as<bool>();
        bool discard_n =  vm["discard-n"].as<bool>();

        WriterHelper writer(pe, se);

        auto read_visit = make_read_visitor_func(
            [&](SingleEndRead *ser) {
                extract_umi( ser->non_const_read_one(), umi_length, qual_threshold, avg_qual_threshold, 
                                                        qual_offset, homopolymer, discard_n );
            },
            [&](PairedEndRead *per) {
                if (read == 'F') {
                    extract_umi( per->non_const_read_one(), umi_length, qual_threshold, avg_qual_threshold, 
                                                            qual_offset, homopolymer, discard_n );
                    extract_umi( per->non_const_read_two(), umi_length, qual_threshold, avg_qual_threshold, 
                                                            qual_offset, homopolymer, discard_n );
                } else {
                    extract_umi( per->non_const_read_two(), umi_length, qual_threshold, avg_qual_threshold, 
                                                            qual_offset, homopolymer, discard_n );
                    extract_umi( per->non_const_read_one(), umi_length, qual_threshold, avg_qual_threshold, 
                                                            qual_offset, homopolymer, discard_n );
                }
            }
            );

        while(reader.has_next()) {
            auto i = reader.next();
            counters.input(*i);
            i->accept(read_visit);
            umi = UMI();
            writer(*i);
            counters.output(*i);
        }
    }
};
#endif
