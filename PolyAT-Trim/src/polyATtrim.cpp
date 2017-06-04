#include <iostream>
#include <string>
#include <boost/program_options.hpp>
#include <vector>
#include <fstream>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/device/file_descriptor.hpp>
#include <boost/iostreams/stream.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/device/mapped_file.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem.hpp>
#include <map>
#include <unordered_map>
#include <boost/functional/hash.hpp>

#include "polyATtrim.h"

namespace
{
    const size_t SUCCESS = 0;
    const size_t ERROR_IN_COMMAND_LINE = 1;
    const size_t ERROR_UNHANDLED_EXCEPTION = 2;

} // namespace

namespace bi = boost::iostreams;

int main(int argc, char** argv)
{

    const std::string program_name = "AT_Trim";

    TrimmingCounters counters;

    try
    {
        /** Define and parse the program options
         */
        namespace po = boost::program_options;
        po::options_description desc("Options");
        setDefaultParams(desc, program_name);
        setDefaultParamsCutting(desc);
        setDefaultParamsTrim(desc);

        desc.add_options()
            ("max-mismatch,x", po::value<size_t>()->default_value(3),    "Max amount of mismatches allowed in trimmed area")
            ("min-trim,M", po::value<size_t>()->default_value(5),    "Min base pairs trim for AT tail");

        po::variables_map vm;
        try
        {
            po::store(po::parse_command_line(argc, argv, desc),
                      vm); // can throw

            version_or_help( program_name, desc, vm);
            
            /** --help option
             */
            po::notify(vm); // throws on error, so do after help in case
            

            std::string statsFile(vm["stats-file"].as<std::string>());
            std::string prefix(vm["prefix"].as<std::string>());

            std::shared_ptr<OutputWriter> pe = nullptr;
            std::shared_ptr<OutputWriter> se = nullptr;

            outputWriters(pe, se, vm["fastq-output"].as<bool>(), vm["tab-output"].as<bool>(), vm["interleaved-output"].as<bool>(), vm["unmapped-output"].as<bool>(), vm["force"].as<bool>(), vm["gzip-output"].as<bool>(), vm["to-stdout"].as<bool>(), prefix );

            if(vm.count("read1-input")) {
                if (!vm.count("read2-input")) {
                    throw std::runtime_error("must specify both read1 and read2 input files.");
                } else if (vm.count("read2-input") != vm.count("read1-input")) {
                    throw std::runtime_error("must have same number of input files for read1 and read2");
                }
                auto read1_files = vm["read1-input"].as<std::vector<std::string> >();
                auto read2_files = vm["read2-input"].as<std::vector<std::string> >();

                for(size_t i = 0; i < read1_files.size(); ++i) {
                    bi::stream<bi::file_descriptor_source> is1{check_open_r(read1_files[i]), bi::close_handle};
                    bi::stream<bi::file_descriptor_source> is2{check_open_r(read2_files[i]), bi::close_handle};
                   
                    InputReader<PairedEndRead, PairedEndReadFastqImpl> ifp(is1, is2);
                    helper_trim(ifp, pe, se, counters, vm["min-length"].as<std::size_t>(), vm["min-trim"].as<std::size_t>(), vm["max-mismatch"].as<std::size_t>(), vm["stranded"].as<bool>(), vm["no-left"].as<bool>(), vm["no-right"].as<bool>(), vm["no-orphans"].as<bool>() );
                }
            }

            if(vm.count("singleend-input")) {
                auto read_files = vm["singleend-input"].as<std::vector<std::string> >();
                for (auto file : read_files) {
                    bi::stream<bi::file_descriptor_source> sef{ check_open_r(file), bi::close_handle};
                    InputReader<SingleEndRead, SingleEndReadFastqImpl> ifs(sef);
                    helper_trim(ifs, pe, se, counters, vm["min-length"].as<std::size_t>(), vm["min-trim"].as<std::size_t>(), vm["max-mismatch"].as<std::size_t>(), vm["stranded"].as<bool>(), vm["no-left"].as<bool>(), vm["no-right"].as<bool>(), vm["no-orphans"].as<bool>() );
                }
            }
            
            if(vm.count("tab-input")) {
                auto read_files = vm["tab-input"].as<std::vector<std::string> > ();
                for (auto file : read_files) {
                    bi::stream<bi::file_descriptor_source> tabin{ check_open_r(file), bi::close_handle};
                    InputReader<ReadBase, TabReadImpl> ift(tabin);
                    helper_trim(ift, pe, se, counters, vm["min-length"].as<std::size_t>(), vm["min-trim"].as<std::size_t>(), vm["max-mismatch"].as<std::size_t>(), vm["stranded"].as<bool>(), vm["no-left"].as<bool>(), vm["no-right"].as<bool>(), vm["no-orphans"].as<bool>() );
                }
            }
            
            if (vm.count("interleaved-input")) {
                auto read_files = vm["interleaved-input"].as<std::vector<std::string > >();
                for (auto file : read_files) {
                    bi::stream<bi::file_descriptor_source> inter{ check_open_r(file), bi::close_handle};
                    InputReader<PairedEndRead, InterReadImpl> ifp(inter);
                    helper_trim(ifp, pe, se, counters, vm["min-length"].as<std::size_t>(), vm["min-trim"].as<std::size_t>(), vm["max-mismatch"].as<std::size_t>(), vm["stranded"].as<bool>(), vm["no-left"].as<bool>(), vm["no-right"].as<bool>(), vm["no-orphans"].as<bool>() );
                }
            }
           
            if (vm.count("std-input")) {
                bi::stream<bi::file_descriptor_source> tabin {fileno(stdin), bi::close_handle};
                InputReader<ReadBase, TabReadImpl> ift(tabin);
                helper_trim(ift, pe, se, counters, vm["min-length"].as<std::size_t>(), vm["min-trim"].as<std::size_t>(), vm["max-mismatch"].as<std::size_t>(), vm["stranded"].as<bool>(), vm["no-left"].as<bool>(), vm["no-right"].as<bool>(), vm["no-orphans"].as<bool>() );
            }  

        }
        catch(po::error& e)
        {
            std::cerr << "ERROR: " << e.what() << std::endl << std::endl;
            std::cerr << desc << std::endl;
            return ERROR_IN_COMMAND_LINE;
        }

    }
    catch(std::exception& e)
    {
        std::cerr << "\n\tUnhandled Exception: "
                  << e.what() << std::endl;
        return ERROR_UNHANDLED_EXCEPTION;

    }

    return SUCCESS;

}
