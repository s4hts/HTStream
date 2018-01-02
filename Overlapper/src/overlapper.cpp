#include <iostream>
#include <string>
#include <boost/program_options.hpp>
#include <vector>
#include <fstream>
#include "ioHandler.h"
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

#include "overlapper.h"

namespace
{
    const size_t SUCCESS = 0;
    const size_t ERROR_IN_COMMAND_LINE = 1;
    const size_t ERROR_UNHANDLED_EXCEPTION = 2;

} // namespace

namespace bi = boost::iostreams;

int main(int argc, char** argv)
{
    const std::string program_name = "overlapper";
    std::string app_description = 
                       "The overlapper application overlaps paired end reads producing\n";

    app_description += "  a longer overlapped single-end read (lins), or in the case of\n";
    app_description += "  complete short overlapping (sins) also removes adapters.\n";
    app_description += "  Reads left not overlapped are lables as 'nins'.\n";

    OverlappingCounters counters;

    try
    {
        /** Define and parse the program options
         */
        namespace po = boost::program_options;
        po::options_description standard = setStandardOptions();
            // version|v ; help|h ; notes|N ; stats-file|L ; append-stats-file|A
        po::options_description input = setInputOptions();
            // read1-input|1 ; read2-input|2 ; singleend-input|U
            // tab-input|T ; interleaved-input|I ; from-stdin|S
        po::options_description output = setOutputOptions(program_name);
            // force|F ; prefix|p ; gzip-output,g ; fastq-output|f
            // tab-output|t ; interleaved-output|i ; unmapped-output|u ; to-stdout,O

        po::options_description desc("Application Specific Options");

        setDefaultParamsCutting(desc);
            // no-orphans|n ; stranded|s ; min-length|m

        setDefaultParamsOverlapping(desc);
            // kmer|k ; kmer-offset|r ; max-mismatch-errorDensity|x
            // check-lengths|c ; min-overlap|o
        
        po::options_description cmdline_options;
        cmdline_options.add(standard).add(input).add(output).add(desc);

        po::variables_map vm;
        try
        {
            po::store(po::parse_command_line(argc, argv, cmdline_options),
                      vm); // can throw
            
            version_or_help(program_name, app_description, cmdline_options, vm);

            po::notify(vm); // throws on error, so do after help in case
         
            if (vm["max-mismatch-errorDensity"].as<double>() < 0.0 ||  vm["max-mismatch-errorDensity"].as<double>() > 1.0) {
                throw std::runtime_error("Woah, there human. It seems you have entered a wacky, zany number that isn't between 0.0 and 1.0 for max mismatch errorDensity.\nThat is why there is this error message. (>^_^)>\n\n");
            }

            if (vm["kmer-offset"].as<size_t>() == 0) {
                throw std::runtime_error("Human - you cannot have a kmer offset of zero! Madness would ensue. Please, try again, but this time make sure the kmer-offset,r flag is set to above 0.\n\n");
            }
        
            std::string statsFile(vm["stats-file"].as<std::string>());
            std::string prefix(vm["prefix"].as<std::string>());
   
            std::shared_ptr<OutputWriter> pe = nullptr;
            std::shared_ptr<OutputWriter> se = nullptr;
            
            outputWriters(pe, se, vm["fastq-output"].as<bool>(), vm["tab-output"].as<bool>(), vm["interleaved-output"].as<bool>(), vm["unmapped-output"].as<bool>(), vm["force"].as<bool>(), vm["gzip-output"].as<bool>(), vm["to-stdout"].as<bool>(), prefix );

            if(vm.count("read1-input")) {
                if (!vm.count("read2-input")) {
                    throw std::runtime_error("must specify both read1 and read2 input files.");
                }
                auto read1_files = vm["read1-input"].as<std::vector<std::string> >();
                auto read2_files = vm["read2-input"].as<std::vector<std::string> >();
                
                if (read1_files.size() != read2_files.size()) {
                    throw std::runtime_error("must have same number of input files for read1 and read2");
                }

                for(size_t i = 0; i < read1_files.size(); ++i) {
                    bi::stream<bi::file_descriptor_source> is1{check_open_r(read1_files[i]), bi::close_handle};
                    bi::stream<bi::file_descriptor_source> is2{check_open_r(read2_files[i]), bi::close_handle};

                    InputReader<PairedEndRead, PairedEndReadFastqImpl> ifp(is1, is2);
                    helper_overlapper(ifp, pe, se, counters, vm["max-mismatch-errorDensity"].as<double>(),  vm["min-overlap"].as<size_t>(), vm["stranded"].as<bool>(), vm["min-length"].as<size_t>(), vm["check-lengths"].as<size_t>(),   vm["kmer"].as<size_t>(), vm["kmer-offset"].as<size_t>(), vm["no-orphans"].as<bool>() ); }
            }

            if(vm.count("singleend-input")) {
                auto read_files = vm["singleend-input"].as<std::vector<std::string> >();
                for (auto file : read_files) {
                    bi::stream<bi::file_descriptor_source> sef{ check_open_r(file), bi::close_handle};
                    InputReader<SingleEndRead, SingleEndReadFastqImpl> ifs(sef);
                    helper_overlapper(ifs, pe, se, counters, vm["max-mismatch-errorDensity"].as<double>(),  vm["min-overlap"].as<size_t>(), vm["stranded"].as<bool>(), vm["min-length"].as<size_t>(), vm["check-lengths"].as<size_t>(),   vm["kmer"].as<size_t>(), vm["kmer-offset"].as<size_t>(), vm["no-orphans"].as<bool>() );
                }
            }
            
            if(vm.count("tab-input")) {
                auto read_files = vm["tab-input"].as<std::vector<std::string> > ();
                for (auto file : read_files) {
                    bi::stream<bi::file_descriptor_source> tabin{ check_open_r(file), bi::close_handle};
                    InputReader<ReadBase, TabReadImpl> ift(tabin);
                    helper_overlapper(ift, pe, se, counters, vm["max-mismatch-errorDensity"].as<double>(),  vm["min-overlap"].as<size_t>(), vm["stranded"].as<bool>(), vm["min-length"].as<size_t>(), vm["check-lengths"].as<size_t>(),   vm["kmer"].as<size_t>(), vm["kmer-offset"].as<size_t>(), vm["no-orphans"].as<bool>() );
                }
            }
            
            if (vm.count("interleaved-input")) {
                auto read_files = vm["interleaved-input"].as<std::vector<std::string > >();
                for (auto file : read_files) {
                    bi::stream<bi::file_descriptor_source> inter{ check_open_r(file), bi::close_handle};
                    InputReader<PairedEndRead, InterReadImpl> ifp(inter);
                    helper_overlapper(ifp, pe, se, counters, vm["max-mismatch-errorDensity"].as<double>(),  vm["min-overlap"].as<size_t>(), vm["stranded"].as<bool>(), vm["min-length"].as<size_t>(), vm["check-lengths"].as<size_t>(),   vm["kmer"].as<size_t>(), vm["kmer-offset"].as<size_t>(), vm["no-orphans"].as<bool>() );
                }
            }
           
            if (vm.count("from-stdin")) {
                bi::stream<bi::file_descriptor_source> tabin {fileno(stdin), bi::close_handle};
                InputReader<ReadBase, TabReadImpl> ift(tabin);
                helper_overlapper(ift, pe, se, counters, vm["max-mismatch-errorDensity"].as<double>(),  vm["min-overlap"].as<size_t>(), vm["stranded"].as<bool>(), vm["min-length"].as<size_t>(), vm["check-lengths"].as<size_t>(),   vm["kmer"].as<size_t>(), vm["kmer-offset"].as<size_t>(), vm["no-orphans"].as<bool>() );
            }  
            counters.write_out(statsFile, vm["append-stats-file"].as<bool>(), program_name, vm["notes"].as<std::string>());
        }
        catch(po::error& e)
        {
            std::cerr << "ERROR: " << e.what() << std::endl << std::endl;
            version_or_help(program_name, app_description, cmdline_options, vm, true);
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

