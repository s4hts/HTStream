#include <iostream>
#include <string>
#include <boost/program_options.hpp>
#include <vector>
#include <fstream>
#include "ioHandler.h"
#include <map>
#include <unordered_map>
#include <boost/dynamic_bitset.hpp>

namespace
{
    const size_t SUCCESS = 0;
    const size_t ERROR_IN_COMMAND_LINE = 1;
    const size_t ERROR_UNHANDLED_EXCEPTION = 2;

} // namespace

typedef std::unordered_map<std::string, size_t> Counter;
typedef std::map <boost::dynamic_bitset<>, std::unique_ptr<ReadBase>> BitMap;

template <class T, class Impl>
void load_map(InputReader<T, Impl> &reader, Counter& counters, BitMap& read_map, size_t start, size_t length) {
    while(reader.has_next()) {
        auto i = reader.next();
        ++counters["TotalRecords"];
        //std::cout << "read id: " << i->get_read().get_id() << std::endl;
        //check for existance, store or compare quality and replace:
        if (auto key=i->get_key(start, length)) {
            if(!read_map.count(*key)) {
                read_map[*key] = std::move(i);
            } else if(i->avg_q_score() > read_map[*key]->avg_q_score()){
                read_map[*key] = std::move(i);
                ++counters["Replaced"];
            }
        } else {  // key had N 
            ++counters["HasN"];
        }
    }
}


int main(int argc, char** argv)
{
    BitMap read_map;
    Counter counters;
    counters["TotalRecords"] = 0;
    counters["Replaced"] = 0;
    counters["HasN"] = 0;
    size_t start, length = 0;
    std::string prefix;
    std::vector<std::string> default_pe = {"PE1", "PE2"};
    bool fastq_out;
    bool tab_out;
    bool std_out;
    
    try
    {
        /** Define and parse the program options
         */
        namespace po = boost::program_options;
        po::options_description desc("Options");
        desc.add_options()
            ("version,v",                  "Version print")
            ("read1-input,1", po::value< std::vector<std::string> >(),
                                           "Read 1 input <comma sep for multiple files>") 
            ("read2-input,2", po::value< std::vector<std::string> >(), 
                                           "Read 2 input <comma sep for multiple files>")
            ("singleend-input,U", po::value< std::vector<std::string> >(),
                                           "Single end read input <comma sep for multiple files>")
            ("tab-input,T", po::value< std::vector<std::string> >(),
                                           "Tab input <comma sep for multiple files>")
            ("interleaved-input,I", po::value< std::vector<std::string> >(),
                                           "Interleaved input I <comma sep for multiple files>")
            ("stdin-input,S",              "STDIN input <MUST BE TAB DELIMITED INPUT>")
            ("start,s", po::value<size_t>(&start)->default_value(10),  "Start location for unique ID <int>")
            ("length,l", po::value<size_t>(&length)->default_value(10), "Length of unique ID <int>")
            ("quality-check-off,q",        "Quality Checking Off First Duplicate seen will be kept")
            ("gzip-output,g",              "Output gzipped")
            ("interleaved-output, i",      "Output to interleaved")
            ("fastq-output,f", po::bool_switch(&fastq_out)->default_value(false), "Fastq format output")
            ("force,F", po::bool_switch()->default_value(true),         "Forces overwrite of files")
            ("tab-output,t", po::bool_switch(&tab_out)->default_value(false),   "Tab-delimited output")
            ("to-stdout,O", po::bool_switch(&std_out)->default_value(false),    "Prints to STDOUT in Tab Delimited")
            ("prefix,p", po::value<std::string>(&prefix)->default_value("output_nodup_"),
                                           "Prefix for outputted files")
            ("log-file,L",                 "Output-Logfile")
            ("no-log,N",                   "No logfile <outputs to stderr>")
            ("help,h",                     "Prints help.");

        po::variables_map vm;
        try
        {
            po::store(po::parse_command_line(argc, argv, desc),
                      vm); // can throw

            /** --help option
             */
            if ( vm.count("help")  || vm.size() == 0)
            {
                std::cout << "Super-Deduper" << std::endl
                          << desc << std::endl;
                return SUCCESS;
            }

            po::notify(vm); // throws on error, so do after help in case
            // there are any problems
            if(vm.count("read1-input")) {
                if (!vm.count("read2-input")) {
                    throw std::runtime_error("must specify both read1 and read2 input files.");
                } else if (vm.count("read2-input") != vm.count("read1-input")) {
                    throw std::runtime_error("must have same number of input files for read1 and read2");
                }
                auto read1_files = vm["read1-input"].as<std::vector<std::string> >();
                auto read2_files = vm["read2-input"].as<std::vector<std::string> >();

                for(size_t i = 0; i < read1_files.size(); ++i) {
                    // todo: check file exists etc
                    std::ifstream read1(read1_files[i], std::ifstream::in);
                    std::ifstream read2(read2_files[i], std::ifstream::in);
                    InputReader<PairedEndRead, PairedEndReadFastqImpl> ifp(read1, read2);
                    load_map(ifp, counters, read_map, start, length);
                }
            }
            else if(vm.count("singleend-input")) {
                auto read_files = vm["singleend-input"].as<std::vector<std::string> >();
                for (auto file : read_files) {
                    std::ifstream read1(file, std::ifstream::in);
                    InputReader<SingleEndRead, SingleEndReadFastqImpl> ifs(read1);
                    load_map(ifs, counters, read_map, start, length);
                }
            }

                if (fastq_out || (! std_out && ! tab_out) ) {
                    default_pe[0] += ".fastq";
                    default_pe[1] += ".fastq";
                    default_pe[0] = prefix + default_pe[0];
                    default_pe[1] = prefix + default_pe[1];
                      
                    std::ofstream out1(default_pe[0], std::ofstream::out);
                    std::ofstream out2(default_pe[1], std::ofstream::out);
                    OutputWriter<PairedEndRead, PairedEndReadOutFastq> ofs(out1, out2);
                    for(auto const &i : read_map) {
                        ofs.write(*dynamic_cast<PairedEndRead*>(i.second.get()));
                    }
                }
            
        }
        catch(po::error& e)
        {
            std::cerr << "ERROR: " << e.what() << std::endl << std::endl;
            std::cerr << desc << std::endl;
            return ERROR_IN_COMMAND_LINE;
        }

        // application code here //

    }
    catch(std::exception& e)
    {
        std::cerr << "Unhandled Exception reached the top of main: "
                  << e.what() << ", application will now exit" << std::endl;
        return ERROR_UNHANDLED_EXCEPTION;

    }

    std::cerr << "TotalRecords:" << counters["TotalRecords"] << "\tReplaced:" << counters["Replaced"]
              << "\tKept:" << read_map.size() << "\tRemoved:" << counters["TotalRecords"] - read_map.size()
              << "\tHasN:" << counters["HasN"] << std::endl;
    return SUCCESS;

}
