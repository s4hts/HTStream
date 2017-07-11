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
    const std::string program_name = "Super-Deduper";

    BitMap read_map;
    
    SuperDeduperCounters counters;
    
    try
    {
        /** Define and parse the program options
         */
        po::options_description desc("Options");
        setDefaultParams(desc, program_name);

        desc.add_options()
            ("start,s", po::value<size_t>()->default_value(10),  "Start location for unique ID <int>")
            ("length,l", po::value<size_t>()->default_value(10), "Length of unique ID <int>")
            ("avg-qual-score,q", po::value<double>()->default_value(20), "Avg quality score to have the read written automatically <int>")
            ("inform-avg-qual-score,a", po::value<double>()->default_value(5), "Avg quality score to considered a read informative <int> "); //I know this says user input is a int, but is actually a double
                   po::variables_map vm;
        try
        {
            po::store(po::parse_command_line(argc, argv, desc),
                      vm); // can throw

            /** --help option
             */

            po::notify(vm); // throws on error, so do after help in case

            version_or_help( program_name, desc, vm);
            
            std::string statsFile(vm["stats-file"].as<std::string>());
            std::string prefix(vm["prefix"].as<std::string>());

            std::shared_ptr<OutputWriter> pe = nullptr;
            std::shared_ptr<OutputWriter> se = nullptr;

            outputWriters(pe, se, vm["fastq-output"].as<bool>(), vm["tab-output"].as<bool>(), vm["interleaved-output"].as<bool>(), vm["unmapped-output"].as<bool>(), vm["force"].as<bool>(), vm["gzip-output"].as<bool>(), vm["to-stdout"].as<bool>(), prefix );

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
            
            if (vm.count("std-input")) {
                bi::stream<bi::file_descriptor_source> tabin {fileno(stdin), bi::close_handle};
                InputReader<ReadBase, TabReadImpl> ift(tabin);
                load_map(ift, counters, read_map, pe, se,vm["avg-qual-score"].as<double>(),  vm["inform-avg-qual-score"].as<double>(),  vm["start"].as<size_t>() - 1, vm["length"].as<size_t>());
            }
            
            for(auto const &i : read_map) {
                if (i.second.get() != nullptr) {
                    writer_helper(i.second.get(), pe, se, false);
                }
            }
            counters.write_out(statsFile, vm["append-stats-file"].as<bool>() , program_name, vm["notes"].as<std::string>());
        } catch(po::error& e) {
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
