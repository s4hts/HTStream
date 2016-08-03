#include <iostream>
#include <string>
#include <boost/program_options.hpp>
#include <vector>
#include <fstream>
#include "ioHandler.h"
#include <map>
#include <boost/dynamic_bitset.hpp>

namespace
{
    const size_t SUCCESS = 0;
    const size_t ERROR_IN_COMMAND_LINE = 1;
    const size_t ERROR_UNHANDLED_EXCEPTION = 2;

} // namespace

typedef std::map<std::string, int> Counter;
typedef std::map <boost::dynamic_bitset<>, std::unique_ptr<ReadBase>> BitMap;

template <class T, class Impl>
void load_map(InputReader<T, Impl> &reader, Counter& counters, BitMap& read_map, size_t start, size_t length) {
    while(reader.has_next()) {
        auto i = reader.next();
        counters["TotalRecords"]++;
        //std::cout << "read id: " << i->get_read().get_id() << std::endl;
        //check for existance, store or compare quality and replace:
        auto key=i->get_key(start, length);
        if(!read_map.count(key)) {
            read_map[key] = std::move(i);
        } else if(i->avg_q_score() > read_map[key]->avg_q_score()){
            read_map[key] = std::move(i);
            counters["Replaced"]++;
        }
    }
}


int main(int argc, char** argv)
{
    BitMap read_map;
    Counter counters;
    counters["TotalRecords"] = 0;
    counters["Replaced"] = 0;
    size_t start, length = 0;

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
            ("start,s", po::value<size_t>(&start),  "Start location for unique ID <int>")
            ("length,l", po::value<size_t>(&length), "Length of unique ID <int>")
            ("quality-check-off,q",        "Quality Checking Off First Duplicate seen will be kept")
            ("gzip-output,g",              "Output gzipped")
            ("interleaved-output, i",      "Output to interleaved")
            ("fastq-output,f",             "Fastq format outputted <R1 and R2>")
            ("force,F",                    "Forces overwrite of files")
            ("tab-output,t",               "Tab-delimited output")
            ("to-stdout,O",                "Prints to STDOUT in Tab Delimited")
            ("prefix,p",                   "Prefix for outputted files")
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
                auto read2_files = vm["read1-input"].as<std::vector<std::string> >();
                for(size_t i = 0; i < read1_files.size(); ++i) {
                    // todo: check file exists etc
                    std::ifstream read1(read1_files[i], std::ifstream::in);
                    std::ifstream read2(read2_files[i], std::ifstream::in);
                    InputReader<PairedEndRead, PairedEndReadImpl> ifp(read1, read2);
                    load_map(ifp, counters, read_map, start, length);
                }
            }
            else if(vm.count("singleend-input")) {
                auto read_files = vm["singleend-input"].as<std::vector<std::string> >();
                for (auto file : read_files) {
                    std::ifstream read1(file, std::ifstream::in);
                    InputReader<SingleEndRead, SingleEndReadImpl> ifs(read1);
                    load_map(ifs, counters, read_map, start, length);
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
              << "\tKept:" << read_map.size() << std::endl;
    return SUCCESS;

}
