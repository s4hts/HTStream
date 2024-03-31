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

    std::string UMI = "";

    ExtractUMI() {
        program_name = "hts_ExtractUMI";
        app_description =
            "The hts_ExtractUMI application trims a set number of bases from the 5'\n";
        app_description += "  end of a read and appends it to the read ID. Program is \n";
        app_description += "  meant to be an HTStream drop-in for umi_tools extract \n";
        app_description += "  that is compatible with HTStream streaming pipelines.\n";
    }

    void add_extra_options(po::options_description &desc) {
        desc.add_options()
            ("read,r", po::value<size_t>()->default_value(1)->notifier(boost::bind(&check_range<size_t>, "read", _1, 1, 2)), "Read from which to extract the UMI, ignored if SE")
            ("umi_length,l", po::value<size_t>()->default_value(6)->notifier(boost::bind(&check_range<size_t>, "umi_length", _1, 1, 36)), "Total length of UMI to extract (1, 36)")
            ("avg-qual-score,q", po::value<size_t>()->default_value(0)->notifier(boost::bind(&check_range<size_t>, "avg-qual-score", _1, 0, 10000)), "Threshold for quality score average of UMI (min 1, max 10000), read pairs are discarded, default is unset");
    }



    void extract_umi(Read &r, size_t umi_length, size_t qual_threshold, size_t qual_offset) {

        qual_threshold += qual_offset;
        std::string qual = r.get_qual();
        size_t current_sum = 0;
        size_t final_qual_threshold = qual_threshold * umi_length;
        size_t i = 0;

        /* 
            This condition is used to transfer UMI or discard status
                between PE reads.
        */
        if (UMI.empty()) {

            if (qual_threshold != qual_offset) { // Threshold of 0 = no quality filtering

                for (std::string::iterator it = qual.begin() ; it != qual.end() ; ++it) {
                    current_sum += static_cast<size_t>(*it);
                    i = static_cast<size_t>(it - qual.begin()) + 1;
                    if (i >= umi_length) {
                        break;
                    }
                }
            }


            if (qual_threshold == qual_offset || current_sum >= final_qual_threshold) { // if threshold is met or unset
                UMI = r.get_seq().substr(0, umi_length);
                r.setLCut(umi_length);
            
            } else {
                UMI = "discard";
            }
        
        }


        if (UMI.compare("discard") == 0) {
            r.setDiscard();
        } else {
            r.set_id_first(r.get_id_first() + "_" + UMI);   
        }   
        
    }


    template <class T, class Impl>
    void do_app(InputReader<T, Impl> &reader, std::shared_ptr<OutputWriter> pe, std::shared_ptr<OutputWriter> se, TrimmingCounters& counters, const po::variables_map &vm) {
        
        size_t read = vm["read"].as<size_t>();
        size_t umi_length = vm["umi_length"].as<size_t>(); 
        size_t qual_offset = vm["qual-offset"].as<size_t>();
        size_t qual_threshold = vm["avg-qual-score"].as<size_t>();

        WriterHelper writer(pe, se);

        auto read_visit = make_read_visitor_func(
            [&](SingleEndRead *ser) {
                extract_umi( ser->non_const_read_one(), umi_length, qual_threshold, qual_offset );
            },
            [&](PairedEndRead *per) {
                if (read == 1) {
                    extract_umi( per->non_const_read_one(), umi_length, qual_threshold, qual_offset );
                    extract_umi( per->non_const_read_two(), umi_length, qual_threshold, qual_offset );
                } else {
                    extract_umi( per->non_const_read_two(), umi_length, qual_threshold, qual_offset );
                    extract_umi( per->non_const_read_one(), umi_length, qual_threshold, qual_offset );
                }
            }
            );

        while(reader.has_next()) {
            auto i = reader.next();
            counters.input(*i);
            i->accept(read_visit);
            UMI = "";
            writer(*i);
            counters.output(*i);
        }
    }
};
#endif
