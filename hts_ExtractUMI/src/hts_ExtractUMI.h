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
#include <iostream>

extern template class InputReader<SingleEndRead, SingleEndReadFastqImpl>;
extern template class InputReader<PairedEndRead, PairedEndReadFastqImpl>;
extern template class InputReader<PairedEndRead, InterReadImpl>;
extern template class InputReader<ReadBase, TabReadImpl>;


class ExtractUMI: public MainTemplate<TrimmingCounters, ExtractUMI> {
public:

    ExtractUMI() {
        program_name = "hts_ExtractUMI";
        app_description =
            "The hts_ExtractUMI application trims a set number of bases from the 5'\n";
        app_description += "  end of a read and appends it to the read ID.\n";
    }

    void add_extra_options(po::options_description &desc) {
        desc.add_options()
            ("read,r", po::value<size_t>()->default_value(1)->notifier(boost::bind(&check_range<size_t>, "read", _1, 1, 2)), "Read from which to extract the UMI.");
        desc.add_options()
            ("umi_length,l", po::value<size_t>()->default_value(6)->notifier(boost::bind(&check_range<size_t>, "umi_length", _1, 1, 36)), "Total length of UMI to extract (1, 36)");
    }



    std::string extract_umi(Read &r, std::string umi = "", size_t umi_length = 6) {
        if (umi == "") {
            umi = r.get_seq().substr(0, umi_length);
            r.setLCut(umi_length);
        }
        r.set_id_first(r.get_id_first() + "_" + umi);       
        return umi;
    }


    template <class T, class Impl>
    void do_app(InputReader<T, Impl> &reader, std::shared_ptr<OutputWriter> pe, std::shared_ptr<OutputWriter> se, TrimmingCounters& counters, const po::variables_map &vm) {
        
        size_t read = vm["read"].as<size_t>();
        size_t umi_length = vm["umi_length"].as<size_t>(); 
        std::string UMI = "";


        WriterHelper writer(pe, se);

        auto read_visit = make_read_visitor_func(
            [&](SingleEndRead *ser) {
                UMI = extract_umi( ser->non_const_read_one(), UMI, umi_length );
            },
            [&](PairedEndRead *per) {
                if (read == 1) {
                    UMI = extract_umi( per->non_const_read_one(), UMI, umi_length );
                    UMI = extract_umi( per->non_const_read_two(), UMI, umi_length );
                } else {
                    UMI = extract_umi( per->non_const_read_two(), UMI, umi_length );
                    UMI = extract_umi( per->non_const_read_one(), UMI, umi_length );
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
