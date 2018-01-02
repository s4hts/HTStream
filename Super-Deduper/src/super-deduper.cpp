#include "super-deduper.h"
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

namespace
{
    const size_t SUCCESS = 0;
    const size_t ERROR_IN_COMMAND_LINE = 1;
    const size_t ERROR_UNHANDLED_EXCEPTION = 2;

} // namespace

namespace bi = boost::iostreams;

int main(int argc, char** argv)
{
    const std::string program_name = "super-deduper";
    std::string app_description =
                       "Super Deduper is a PCR duplicate remover from sequence data. It uses a\n";
    app_description += "  subsequence within each read to detect duplicates.\n";
    app_description += "  Reads with 'N' character(s) isn the key sequence are ignored.\n";    
    app_description += "  Not recommended for single-end reads.";


    BitMap read_map;
    
    SuperDeduperCounters counters;
    
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

        desc.add_options()
            ("start,s", po::value<size_t>()->default_value(10)->notifier(boost::bind(&check_range<size_t>, "start", _1, 1, 10000)),  "Start location for unique ID (min 1, max 10000)")
            ("length,l", po::value<size_t>()->default_value(10)->notifier(boost::bind(&check_range<size_t>, "length", _1, 1, 10000)), "Length of unique ID (min 1, max 10000)")
            ("avg-qual-score,q", po::value<double>()->default_value(30)->notifier(boost::bind(&check_range<size_t>, "avg-qual-score", _1, 1, 10000)), "Avg quality score to have the read written automatically (min 1, max 10000)")
            ("inform-avg-qual-score,a", po::value<double>()->default_value(5)->notifier(boost::bind(&check_range<size_t>, "inform-avg-qual-score", _1, 1, 10000)), "Avg quality score to considered a read informative (min 1, max 10000)");

        po::options_description cmdline_options;
        cmdline_options.add(standard).add(input).add(output).add(desc);

        po::variables_map vm;
        try
        {
            po::store(po::parse_command_line(argc, argv, cmdline_options),
                      vm); // can throw

            /** --help option
             */

            po::notify(vm); // throws on error, so do after help in case

            version_or_help( program_name, app_description, cmdline_options, vm);
            
            std::string statsFile(vm["stats-file"].as<std::string>());
            std::string prefix(vm["prefix"].as<std::string>());

            std::shared_ptr<OutputWriter> pe = nullptr;
            std::shared_ptr<OutputWriter> se = nullptr;

            outputWriters(pe, se, vm["fastq-output"].as<bool>(), vm["tab-output"].as<bool>(), vm["interleaved-output"].as<bool>(), vm["unmapped-output"].as<bool>(), vm["force"].as<bool>(), vm["gzip-output"].as<bool>(), vm["to-stdout"].as<bool>(), prefix );

            // there are any problems
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
                    load_map(ifp, counters, read_map, pe, se,vm["avg-qual-score"].as<double>(), vm["inform-avg-qual-score"].as<double>(), vm["start"].as<size_t>() - 1, vm["length"].as<size_t>());
                }
            }

            if(vm.count("singleend-input")) {
                auto read_files = vm["singleend-input"].as<std::vector<std::string> >();
                for (auto file : read_files) {
                    bi::stream<bi::file_descriptor_source> ser{ check_open_r(file), bi::close_handle};
                    InputReader<SingleEndRead, SingleEndReadFastqImpl> ifs(ser);
                    load_map(ifs, counters, read_map, pe, se,vm["avg-qual-score"].as<double>(), vm["inform-avg-qual-score"].as<double>(), vm["start"].as<size_t>() - 1, vm["length"].as<size_t>());
                }
            }
            
            if(vm.count("tab-input")) {
                auto read_files = vm["tab-input"].as<std::vector<std::string> > ();
                for (auto file : read_files) {
                    bi::stream<bi::file_descriptor_source> tabin{ check_open_r(file), bi::close_handle};
                    InputReader<ReadBase, TabReadImpl> ift(tabin);
                    load_map(ift, counters, read_map, pe, se,vm["avg-qual-score"].as<double>(),  vm["inform-avg-qual-score"].as<double>(),  vm["start"].as<size_t>() - 1, vm["length"].as<size_t>());
                }
            }
            
            if (vm.count("interleaved-input")) {
                auto read_files = vm["interleaved-input"].as<std::vector<std::string > >();
                for (auto file : read_files) {
                    bi::stream<bi::file_descriptor_source> inter{ check_open_r(file), bi::close_handle};
                    InputReader<PairedEndRead, InterReadImpl> ifp(inter);
                    load_map(ifp, counters, read_map, pe, se,vm["avg-qual-score"].as<double>(),  vm["inform-avg-qual-score"].as<double>(), vm["start"].as<size_t>() - 1, vm["length"].as<size_t>());
                }
            }
            
            if (vm.count("from-stdin")) {
                bi::stream<bi::file_descriptor_source> tabin {fileno(stdin), bi::close_handle};
                InputReader<ReadBase, TabReadImpl> ift(tabin);
                load_map(ift, counters, read_map, pe, se,vm["avg-qual-score"].as<double>(),  vm["inform-avg-qual-score"].as<double>(),  vm["start"].as<size_t>() - 1, vm["length"].as<size_t>());
            }
            
            for(auto const &i : read_map) {
                if (i.second.get() != nullptr) {
                    counters.output(*i.second.get()); 
                    writer_helper(i.second.get(), pe, se, false);
                }
            }
            counters.write_out(statsFile, vm["append-stats-file"].as<bool>() , program_name, vm["notes"].as<std::string>());
        } catch(po::error& e) {
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
