#ifndef MAIN_TEMPLATE
#define MAIN_TEMPLATE
//  this is so we can implment hash function for dynamic_bitset
#define BOOST_DYNAMIC_BITSET_DONT_USE_FRIENDS

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
#include <algorithm>
#include "utils.h"

extern template class InputReader<SingleEndRead, SingleEndReadFastqImpl>;
extern template class InputReader<PairedEndRead, PairedEndReadFastqImpl>;
extern template class InputReader<PairedEndRead, InterReadImpl>;
extern template class InputReader<ReadBase, TabReadImpl>;
namespace po = boost::program_options;
namespace bi = boost::iostreams;

template<typename CounterType, typename DerivedType>
class MainTemplate {
public:
    const size_t SUCCESS = 0;
    const size_t ERROR_IN_COMMAND_LINE = 1;
    const size_t ERROR_UNHANDLED_EXCEPTION = 2;
    const size_t ERROR_HTS_EXCEPTION = 3;

    std::string program_name;
    std::string app_description;

    int main_func(int argc, char** argv) {
        try
        {
            /** Define and parse the program options
             */
            po::options_description standard = setStandardOptions();
            // version|v ; help|h ; notes|N ; stats-file|L ; append-stats-file|A
            po::options_description input = setInputOptions();
            // read1-input|1 ; read2-input|2 ; singleend-input|U
            // tab-input|T ; interleaved-input|I
            po::options_description output = setOutputOptions(program_name);
            // force|F ; uncompressed|u ; fastq-output|f
            // tab-output|t ; interleaved-output|i ; unmapped-output|z

            po::options_description desc("Application Specific Options");

            static_cast<DerivedType*>(this)->add_extra_options(desc);

            po::options_description cmdline_options;
            cmdline_options.add(standard).add(input).add(output).add(desc);
            po::variables_map vm;

            try
            {
                po::store(po::parse_command_line(argc, argv, cmdline_options), vm); // can throw

                /** --help option
                 */
                version_or_help(program_name, app_description, cmdline_options, vm);
                po::notify(vm); // throws on error, so do after help in case

                std::shared_ptr<OutputWriter> pe = nullptr;
                std::shared_ptr<OutputWriter> se = nullptr;
                outputWriters(pe, se, vm);

                CounterType counters(program_name, vm);

                if(vm.count("read1-input")) {
                    if (!vm.count("read2-input")) {
                        throw HtsRuntimeException("must specify both read1 and read2 input files.");
                    }
                    auto read1_files = vm["read1-input"].as<std::vector<std::string> >();
                    auto read2_files = vm["read2-input"].as<std::vector<std::string> >();
                    if (read1_files.size() != read2_files.size()) {
                        throw HtsRuntimeException("must have same number of input files for read1 and read2");
                    }
                    for(size_t i = 0; i < read1_files.size(); ++i) {
                        bi::stream<bi::file_descriptor_source> is1{check_open_r(read1_files[i]), bi::close_handle};
                        bi::stream<bi::file_descriptor_source> is2{check_open_r(read2_files[i]), bi::close_handle};
                        InputReader<PairedEndRead, PairedEndReadFastqImpl> ifp(is1, is2);
                        static_cast<DerivedType*>(this)->do_app(ifp, pe, se, counters, vm);
                    }
                }
                if (vm.count("interleaved-input")) {
                    auto read_files = vm["interleaved-input"].as<std::vector<std::string > >();
                    for (auto file : read_files) {
                        bi::stream<bi::file_descriptor_source> inter{ check_open_r(file), bi::close_handle};
                        InputReader<PairedEndRead, InterReadImpl> ifi(inter);
                        static_cast<DerivedType*>(this)->do_app(ifi, pe, se, counters, vm);
                    }
                }
                if(vm.count("singleend-input")) {
                    auto read_files = vm["singleend-input"].as<std::vector<std::string> >();
                    for (auto file : read_files) {
                        bi::stream<bi::file_descriptor_source> sef{ check_open_r(file), bi::close_handle};
                        InputReader<SingleEndRead, SingleEndReadFastqImpl> ifs(sef);
                        static_cast<DerivedType*>(this)->do_app(ifs, pe, se, counters, vm);
                    }
                }
                if(vm.count("tab-input")) {
                    auto read_files = vm["tab-input"].as<std::vector<std::string> > ();
                    for (auto file : read_files) {
                        bi::stream<bi::file_descriptor_source> tabin{ check_open_r(file), bi::close_handle};
                        InputReader<ReadBase, TabReadImpl> ift(tabin);
                        static_cast<DerivedType*>(this)->do_app(ift, pe, se, counters, vm);
                    }
                }
                if (!isatty(fileno(stdin))) {
                    bi::stream<bi::file_descriptor_source> tabin {fileno(stdin), bi::close_handle};
                    InputReader<ReadBase, TabReadImpl> ift(tabin);
                    static_cast<DerivedType*>(this)->do_app(ift, pe, se, counters, vm);
                }

                counters.write_out();
            }
            catch(po::error& e)
            {
                std::cerr << "ERROR: " << e.what() << std::endl << std::endl;
                return ERROR_IN_COMMAND_LINE;
            }

        }
        catch(HtsException& e)
        {
            std::cerr << "ERROR: " << e.what() << std::endl;
            return ERROR_HTS_EXCEPTION;

        }
        catch(std::exception& e)
        {
            std::cerr << "ERROR: Unhandled Exception: "
                      << e.what() << std::endl;
            return ERROR_UNHANDLED_EXCEPTION;

        }
        return SUCCESS;
    }
};

#endif
