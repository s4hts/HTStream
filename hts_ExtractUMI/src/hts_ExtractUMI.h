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
        std::string seq;
        std::string qual;
        size_t length;
        bool discard;
        size_t qual_threshold;
        size_t avg_qual_threshold;
        size_t qual_offset;
        bool homopolymer;
        bool discard_n;
        bool add_as_tag;
    };


    ExtractUMI() {
        program_name = "hts_ExtractUMI";
        app_description =
            "The hts_ExtractUMI application trims a set number of bases from the 5'\n";
        app_description += "   (left) end of a read and appends it to the end of the read ID.\n";
    }

    void add_extra_options(po::options_description &desc) {
        desc.add_options()
        ("read,r", po::value<char>()->default_value('F')->notifier(boost::bind(&check_values<char>, "read", _1, READ_OPTIONS)), "Read from which to extract the UMI (F = Forward, R = Reverse, B = Both, P = Paired End), ignored if SE")
        ("umi-length,l", po::value<size_t>()->default_value(6)->notifier(boost::bind(&check_range<size_t>, "umi_length", _1, 1, 36)), "Total length of UMI to extract (1, 36)")
        ("delimiter,d", po::value<char>()->default_value('_')->notifier(boost::bind(&check_values<char>, "delimiter", _1, DEL_OPTIONS)), "Character to separate the UMI sequence from other fields in the Read ID (Possible options: '-', '_', ':'). Ignored if --add-as-tag is used")
        ("separator,s", po::value<char>()->default_value('+')->notifier(boost::bind(&check_values<char>, "separator", _1, DEL_OPTIONS)), "Character to separate the UMI sequence when using Paired End mode (Possible options: '-', '_', ':').")
        ("qual-score,q", po::value<size_t>()->default_value(0)->notifier(boost::bind(&check_range<size_t>, "qual-score", _1, 0, 10000)), "Threshold for quality score for any base within a UMI (min 1, max 10000), read pairs are discarded, default is unset")
        ("avg-qual-score,Q", po::value<size_t>()->default_value(0)->notifier(boost::bind(&check_range<size_t>, "avg-qual-score", _1, 0, 10000)), "Threshold for quality score average of UMI (min 1, max 10000), read pairs are discarded, default is unset")
        ("homopolymer,p", po::bool_switch()->default_value(false), "Remove reads with homopolymer UMIs")
        ("discard-n,n", po::bool_switch()->default_value(false), "Remove reads with UMIs containing an N")
        ("add-as-tag,C", po::bool_switch()->default_value(false), "Appends UMI to tag section of read IDs.")
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


    void set_dragen(Read &r, const UMI &umi) {
        std::vector<std::string> result;
        boost::split(result, r.get_id_first(), boost::is_any_of(":"));
        if (result.size() < 7) {
            throw HtsIOException("DRAGEN read ID format not found. Must have at least 7 fields deliminated by a ':'");               
        } else if (result.size() == 7) {
            result.push_back(umi.seq);
        } else {
            result[7] = umi.seq;
        }

        r.set_id_first(boost::algorithm::join(result, ":"));
    }

    void set_tag(Read &r, const UMI &umi) {
        r.add_comment("RX:Z:" + umi.seq);
    }

    void set_paired(Read &r, const UMI &umi, const char &del) {
        r.set_id_first(r.get_id_first() + del + umi.seq); // append UMI to read ID
    }

    void extract_umi(Read &r, UMI &umi, const char &del, const char &sep, const bool &dragen = false, const bool &paired_reads = false) {

        std::string tmp_seq;

        if ((umi.seq.empty() || dragen) || paired_reads ) {

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

            if ((!umi.seq.empty() && dragen) || (!umi.seq.empty() && paired_reads)) {
                umi.seq = umi.seq + sep + tmp_seq;
            } else {
                umi.seq = tmp_seq;
            }
        }

        if (!umi.discard) {
            if (!dragen && !umi.add_as_tag && !paired_reads) {
                r.set_id_first(r.get_id_first() + del + umi.seq); // append UMI to read ID
            }
        } else {
            r.setDiscard();
        }

    }

    template <class T, class Impl>
    void do_app(InputReader<T, Impl> &reader, std::shared_ptr<OutputWriter> pe, std::shared_ptr<OutputWriter> se, TrimmingCounters& counters, const po::variables_map &vm) {

        char read = vm["read"].as<char>();
        size_t umi_length = vm["umi-length"].as<size_t>();
        bool add_as_tag = vm["add-as-tag"].as<bool>();
        char del = vm["delimiter"].as<char>();
        char sep = vm["separator"].as<char>();
        size_t qual_offset = vm["qual-offset"].as<size_t>();
        size_t qual_threshold = vm["qual-score"].as<size_t>();
        size_t avg_qual_threshold = vm["avg-qual-score"].as<size_t>();
        bool homopolymer = vm["homopolymer"].as<bool>();
        bool discard_n = vm["discard-n"].as<bool>();
        bool dragen = vm["DRAGEN"].as<bool>();


        if (dragen & (del != ':' || sep != '+')) {
            throw HtsIOException("Delimiter (--delimiter) must be ':' and Separator (--separator) was be '+' to be compatible with --DRAGEN parameter");    
        }


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
            discard_n,              // discard N containing UMIs
            add_as_tag              // add as tag
        };


        WriterHelper writer(pe, se);

        auto read_visit = make_read_visitor_func(
        [&](SingleEndRead * ser) {
            if (dragen && (umi_length > 15)) {
                throw HtsIOException("UMI length (--umi-length) greater than 15 is not compatible with --DRAGEN parameter for Single End Reads");    
            }
            extract_umi( ser->non_const_read_one(), umi, del, sep, dragen );
            if (dragen) { set_dragen( ser->non_const_read_one(), umi ); }
            if (add_as_tag) { set_tag( ser->non_const_read_one(), umi ); }
        },
        [&](PairedEndRead * per) {
            if (dragen && (umi_length > 8)) {
                throw HtsIOException("UMI length (--umi-length) greater than 8 is not compatible with --DRAGEN parameter for PairedEndRead End Reads");    
            }
            if (dragen && add_as_tag) {
                throw HtsIOException("DRAGEN parameter (--DRAGEN) is not compatible with --add-as-tag parameter");    
            }
            if (read == 'F') {
                extract_umi( per->non_const_read_one(), umi, del, sep );
                extract_umi( per->non_const_read_two(), umi, del, sep );
            
            } else if (read == 'R') {
                extract_umi( per->non_const_read_two(), umi, del, sep );
                extract_umi( per->non_const_read_one(), umi, del, sep );
            
            } else if (read == 'B') {
                extract_umi( per->non_const_read_one(), umi, del, sep );
                if (dragen) { set_dragen( per->non_const_read_one(), umi ); }
                if (add_as_tag) { set_tag( per->non_const_read_one(), umi ); }
                std::tie(umi.seq) = std::make_tuple(""); // reset umi struct
                extract_umi( per->non_const_read_two(), umi, del, sep );
                if (dragen) { set_dragen( per->non_const_read_two(), umi ); }
                if (add_as_tag) { set_tag( per->non_const_read_two(), umi ); }
                dragen = false;
                add_as_tag = false; 
            
            } else if (read == 'P') {
                extract_umi( per->non_const_read_one(), umi, del, sep, dragen, true );
                extract_umi( per->non_const_read_two(), umi, del, sep, dragen, true );
                if (!dragen && !add_as_tag) {
                    set_paired( per->non_const_read_one(), umi, del ); 
                    set_paired( per->non_const_read_two(), umi, del ); 
                }   
            }
            if (dragen) {
                set_dragen( per->non_const_read_one(), umi );
                set_dragen( per->non_const_read_two(), umi );
            } 
            if (add_as_tag) {
                set_tag( per->non_const_read_one(), umi );
                set_tag( per->non_const_read_two(), umi );
            }
        }
);

        while (reader.has_next()) {
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
