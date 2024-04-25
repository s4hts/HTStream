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
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/join.hpp>


#include <algorithm>

extern template class InputReader<SingleEndRead, SingleEndReadFastqImpl>;
extern template class InputReader<PairedEndRead, PairedEndReadFastqImpl>;
extern template class InputReader<PairedEndRead, InterReadImpl>;
extern template class InputReader<ReadBase, TabReadImpl>;


class ExtractUMI: public MainTemplate<TrimmingCounters, ExtractUMI> {
public:

    // bases some parameters and allows passing UMIs between PE reads
    struct UMI {
        std::string seq1;
        std::string seq2;
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
        app_description += "   (left) end of a read and appends it to the end of the read ID.\n";
    }

    void add_extra_options(po::options_description &desc) {
        desc.add_options()
        ("read,r", po::value<char>()->default_value('F')->notifier(boost::bind(&check_values<char>, "read", _1, READ_OPTIONS)), "Read from which to extract the UMI (F = Forward, R = Reverse, B = Both), ignored if SE")
        ("umi-length,l", po::value<size_t>()->default_value(6)->notifier(boost::bind(&check_range<size_t>, "umi_length", _1, 1, 36)), "Total length of UMI to extract (1, 36)")
        ("delimiter,d", po::value<char>()->default_value('_')->notifier(boost::bind(&check_values<char>, "delimiter", _1, DEL_OPTIONS)), "Character to separate the UMI sequence from other fields in the Read ID (Possible options: '-', '_', ':')")
        ("qual-score,q", po::value<size_t>()->default_value(0)->notifier(boost::bind(&check_range<size_t>, "qual-score", _1, 0, 10000)), "Threshold for quality score for any base within a UMI (min 1, max 10000), read pairs are discarded, default is unset")
        ("avg-qual-score,Q", po::value<size_t>()->default_value(0)->notifier(boost::bind(&check_range<size_t>, "avg-qual-score", _1, 0, 10000)), "Threshold for quality score average of UMI (min 1, max 10000), read pairs are discarded, default is unset")
        ("homopolymer,p", po::bool_switch()->default_value(false), "Remove reads with homopolymer UMIs")
        ("discard-n,n", po::bool_switch()->default_value(false), "Remove reads with UMIs containing an N")
        ("DRAGEN,D", po::bool_switch()->default_value(false), "Formats UMI addition to Read ID so that it is compatible with Illumina's DRAGEN suite.");
    }


    bool quality_check(UMI &umi) {

        size_t tmp_qual_threshold; // just to be sure original parameters are not overwritten

        if (umi.qual_threshold != 0) {

            tmp_qual_threshold = umi.qual_threshold + umi.qual_offset;

            for (std::string::iterator it = umi.qual.begin() ; it != umi.qual.end() ; ++it) {
                if (static_cast<size_t>(*it) < tmp_qual_threshold) {
                    return true;
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
                return true;
            }
        }

        return false;
    }


    bool homopolymer_check(const std::string &tmp_seq) {
        if (tmp_seq.find_first_not_of(tmp_seq[0]) == std::string::npos) { return true; }
        return false;
    }


    bool n_check(const std::string &tmp_seq) {
        if (tmp_seq.find('N') != std::string::npos) { return true; }
        return false;
    }


    void set_dragen(Read &r, const UMI &umi, const bool &se = false) {

<<<<<<< HEAD
         std::string new_id;
=======
        std::string new_id;
>>>>>>> hts_SuperDeduper
        if (se) {
            new_id = umi.seq1;
        } else {
            new_id = umi.seq1 + "+" + umi.seq2;
        }

        std::vector<std::string> result;
        boost::split(result, r.get_id_first(), boost::is_any_of(":"));
        if (result.size() < 7) {
            throw HtsIOException("DRAGEN read ID format not found. Must have at least 7 fields deliminated by a ':'");               
<<<<<<< HEAD
        } else if (result.size() == 7) {
=======
        } if (result.size() == 7) {
>>>>>>> hts_SuperDeduper
            result.push_back(new_id);
        } else {
            result[7] = new_id;
        }
<<<<<<< HEAD
        r.set_id_first(boost::algorithm::join(result, ":"));
=======
       r.set_id_first(boost::algorithm::join(result, ":"));
>>>>>>> hts_SuperDeduper
    }


    void extract_umi(Read &r, UMI &umi, const char &del, const bool &dragen = false, const bool &both_reads = false) {

        std::string tmp_seq;

        if ( ((umi.seq1.empty()) || dragen) || both_reads ) {

            tmp_seq = r.get_seq().substr(0, umi.length);

            if (umi.qual_threshold + umi.avg_qual_threshold != 0) {
                umi.qual = r.get_qual().substr(0, umi.length);
                umi.discard = quality_check(umi);
            }

            if (umi.homopolymer) {
                umi.discard = homopolymer_check(tmp_seq);
            }

            if (umi.discard_n) {
                umi.discard = n_check(tmp_seq);
            }

            if (!umi.discard) {
                r.setLCut(umi.length);
            }

            if (!umi.seq1.empty() && dragen) {
                umi.seq2 = tmp_seq; // necessary copy?
            } else {
                umi.seq1 = tmp_seq; // necessary copy?
            }

        }

        if ((!umi.discard)) {
            if (!dragen) { 
                r.set_id_first(r.get_id_first() + del + umi.seq1); 
            }
        } else {
            r.setDiscard();
        }

    }

    template <class T, class Impl>
    void do_app(InputReader<T, Impl> &reader, std::shared_ptr<OutputWriter> pe, std::shared_ptr<OutputWriter> se, TrimmingCounters& counters, const po::variables_map &vm) {

        char read = vm["read"].as<char>();
        size_t umi_length = vm["umi-length"].as<size_t>();
        char del = vm["delimiter"].as<char>();
        size_t qual_offset = vm["qual-offset"].as<size_t>();
        size_t qual_threshold = vm["qual-score"].as<size_t>();
        size_t avg_qual_threshold = vm["avg-qual-score"].as<size_t>();
        bool homopolymer = vm["homopolymer"].as<bool>();
        bool discard_n = vm["discard-n"].as<bool>();
        bool dragen = vm["DRAGEN"].as<bool>();


        if (dragen & (del != ':')) {
            throw HtsIOException("Delimiter (--delimiter) must be ':' to be compatible with --DRAGEN parameter");    
        }


        // init UMI struct
        UMI umi = {
            "",                     // umi R1 sequence
            "",                     // umi R2 sequence
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
        [&](SingleEndRead * ser) {
            if (dragen && (umi_length > 15)) {
                throw HtsIOException("UMI length (--umi-length) greater than 15 is not compatible with --DRAGEN parameter for Single End Reads");    
            }
<<<<<<< HEAD
            extract_umi( ser->non_const_read_one(), umi, del, dragen );
            if (dragen) { set_dragen( ser->non_const_read_one(), umi, true); }
=======
            extract_umi( ser->non_const_read_one(), umi, del );
            if (dragen) { set_dragen( ser->non_const_read_one(), umi, true ); }
>>>>>>> hts_SuperDeduper
        },
        [&](PairedEndRead * per) {
            if (dragen && (umi_length > 8)) {
                throw HtsIOException("UMI length (--umi-length) greater than 8 is not compatible with --DRAGEN parameter for PairedEndRead End Reads");    
            }
            if (read == 'F') {
                extract_umi( per->non_const_read_one(), umi, del );
                extract_umi( per->non_const_read_two(), umi, del );
            } else if (read == 'R') {
                extract_umi( per->non_const_read_two(), umi, del );
                extract_umi( per->non_const_read_one(), umi, del );
            } else {
                extract_umi( per->non_const_read_one(), umi, del, dragen, true );
                std::tie(umi.qual) = std::make_tuple(""); // reset umi struct
                extract_umi( per->non_const_read_two(), umi, del, dragen, true );
                if (dragen) {
                    set_dragen( per->non_const_read_one(), umi );
                    set_dragen( per->non_const_read_two(), umi );
                }
            }
        }
                          );

        while (reader.has_next()) {
            auto i = reader.next();
            counters.input(*i);
            i->accept(read_visit);
            std::tie(umi.seq1, umi.seq2, umi.qual, umi.discard) = std::make_tuple("", "", "", false); // reset umi struct
            writer(*i);
            counters.output(*i);
        }
    }
};
#endif
