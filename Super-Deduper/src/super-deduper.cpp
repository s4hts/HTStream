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

int main(int argc, char** argv)
{
    std::map <boost::dynamic_bitset<>, std::unique_ptr<ReadBase>> read_map;
    std::map <std::string, int> counters;
    counters["TotalRecords"] = 0;
    counters["Replaced"] = 0;

    try
    {
        /** Define and parse the program options
         */
        namespace po = boost::program_options;
        po::options_description desc("Options");
        desc.add_options()
            ("help,h", "Prints help.")
            ("read1-input,1", po::value< std::vector<std::string> >(),
             "Read 1 input <comma sep for multiple files>") 
            ("read2-input,2", po::value< std::vector<std::string> >(), 
             "Read 2 input <comma sep for multiple files>")
            ("singleend-input,U","Single end read input <comma sep for multiple files>")
            ("tab-input,T", "Tab input <comma sep for multiple files>")
            ("interleaved-input,I",  "Interleaved input I <comma sep for multiple files>")
            ("stdin-input,S", "STDIN input <MUST BE TAB DELIMITED INPUT>")
            ("start,s", po::value<int>(),  "Start location for unique ID <int>")
            ("length,l", po::value<int>(), "Length of unique ID <int>");

        po::variables_map vm;
        try
        {
            po::store(po::parse_command_line(argc, argv, desc),
                      vm); // can throw

            /** --help option
             */
            if ( vm.count("help")  )
            {
                std::cout << "Super-Deduper" << std::endl
                          << desc << std::endl;
                return SUCCESS;
            }

            po::notify(vm); // throws on error, so do after help in case
            // there are any problems
            if(vm.count("read1-input")) {
                auto read1_files = vm["read1-input"].as<std::vector<std::string> >();
                std::ifstream read1(read1_files[0], std::ifstream::in);
                InputReader<SingleEndRead, SingleEndReadImpl> ifs(read1);

                while(ifs.has_next()) {
                    auto i = ifs.next();
                    counters["TotalRecords"]++;
                    //std::cout << "read id: " << i->get_read().get_id() << std::endl;
                    //check for existance, store or compare quality and replace:
                    auto key=i->getKey(10, 10);
                    if(!read_map.count(key)) {
                        read_map[key] = std::move(i);
                    } else if(i->avg_q_score() > read_map[key]->avg_q_score()){
                        read_map[key] = std::move(i);
                        counters["Replaced"]++;
                      }
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

    std::cout << "TotalRecords:" << counters["TotalRecords"] << "\tReplaced:" << counters["Replaced"]
              << "\tKept:" << read_map.size() << std::endl;
    return SUCCESS;

}
