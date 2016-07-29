#include <iostream> 
#include <string> 
#include <boost/program_options.hpp>
#include <vector>
#include <fstream>
#include "ioHandler.h"

namespace 
{ 
    const size_t SUCCESS = 0; 
    const size_t ERROR_IN_COMMAND_LINE = 1; 
    const size_t ERROR_UNHANDLED_EXCEPTION = 2; 
 
} // namespace 
 
int main(int argc, char** argv) 
{ 
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
                inputFastqSingle ifs(read1);
                
                for (auto i = ifs.begin(); i != ifs.end(); ++i) {
                    std::cout << "read id: " << i->get_read().get_id() << std::endl;
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
 
    return SUCCESS; 

}
