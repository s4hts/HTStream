//  this is so we can implment hash function for dynamic_bitset
#define BOOST_DYNAMIC_BITSET_DONT_USE_FRIENDS

#include <iostream>
#include <string>
#include <boost/program_options.hpp>
#include <vector>
#include <fstream>
#include "ioHandler.h"
#include <map>
#include <unordered_map>
#include <boost/dynamic_bitset.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/device/file_descriptor.hpp>
#include <boost/iostreams/stream.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/device/mapped_file.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem.hpp>
#include <boost/functional/hash.hpp>

namespace
{
    const size_t SUCCESS = 0;
    const size_t ERROR_IN_COMMAND_LINE = 1;
    const size_t ERROR_UNHANDLED_EXCEPTION = 2;

} // namespace

class dbhash {
public:
    std::size_t operator() (const boost::dynamic_bitset<>& bs) const {
        return boost::hash_value(bs.m_bits);
    }
};

typedef std::unordered_map<std::string, size_t> Counter;
typedef std::unordered_map <boost::dynamic_bitset<>, std::unique_ptr<ReadBase>, dbhash> BitMap;

template <class T, class Impl>
void load_map(InputReader<T, Impl> &reader, Counter& counters, BitMap& read_map, size_t start, size_t length) {
    while(reader.has_next()) {
        auto i = reader.next();
        ++counters["TotalRecords"];
        //check for existance, store or compare quality and replace:
        if (auto key=i->get_key(start, length)) {
            // find faster than count on some compilers
            if(read_map.find(*key) == read_map.end()) {
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

template <class T>
void output_read_map_tab(const BitMap& read_map, T& out1) {
    OutputWriter<ReadBase, ReadBaseOutTab> tabs(out1);
    for(auto const &i : read_map) {
        ReadBase* rb = dynamic_cast<ReadBase*>(i.second.get());
        tabs.write(*rb);
    }
}

template <class T>
void output_read_map_inter(const BitMap& read_map, T& out1, T& single) {
    OutputWriter<PairedEndRead, PairedEndReadOutInter> pofs(out1);
    OutputWriter<SingleEndRead, SingleEndReadOutFastq> sofs(single);
    for(auto const &i : read_map) {
        PairedEndRead* per = dynamic_cast<PairedEndRead*>(i.second.get());
        if (per) {
            pofs.write(*per);
        } else {
            SingleEndRead* ser = dynamic_cast<SingleEndRead*>(i.second.get());
            if(ser) {
                sofs.write(*ser);
            }
            else {
                throw std::runtime_error("Unkown read found");
            }
        }
    }
}

template <class T>
void output_read_map_fastq(const BitMap& read_map, T& out1, T& out2, T& single) {

    OutputWriter<PairedEndRead, PairedEndReadOutFastq> pofs(out1, out2);
    OutputWriter<SingleEndRead, SingleEndReadOutFastq> sofs(single);
    for(auto const &i : read_map) {
        PairedEndRead* per = dynamic_cast<PairedEndRead*>(i.second.get());
        if (per) {
            pofs.write(*per);
        } else {
            SingleEndRead* ser = dynamic_cast<SingleEndRead*>(i.second.get());
            if(ser) {
                sofs.write(*ser);
            }
            else {
                throw std::runtime_error("Unkown read found");
            }
        }
    }
}

namespace bi = boost::iostreams;
namespace bf = boost::filesystem;

int check_open_r(const std::string& filename) {
    bf::path p(filename);
    if (!bf::exists(p)) {
        throw std::runtime_error("File " + filename + " was not found.");
    }
    
    if (p.extension() == ".gz") {
        return fileno(popen(("gunzip -c " + filename).c_str(), "r"));
    } else {
        return fileno(fopen(filename.c_str(), "r"));
    }
}

int main(int argc, char** argv)
{
    BitMap read_map;
    Counter counters;
    counters["TotalRecords"] = 0;
    counters["Replaced"] = 0;
    counters["HasN"] = 0;
    size_t start = 0, length = 0;
    std::string prefix;
    std::vector<std::string> default_outfiles = {"PE1", "PE2", "SE"};

    bool fastq_out;
    bool tab_out;
    bool std_out;
    bool gzip_out;
    bool interleaved_out;
    
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
            ("gzip-output,g", po::bool_switch(&gzip_out)->default_value(false),  "Output gzipped")
            ("interleaved-output, i", po::bool_switch(&interleaved_out)->default_value(false),     "Output to interleaved")
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
                    bi::stream<bi::file_descriptor_source> is1{check_open_r(read1_files[i]), bi::close_handle};
                    bi::stream<bi::file_descriptor_source> is2{check_open_r(read2_files[i]), bi::close_handle};
                    
                    InputReader<PairedEndRead, PairedEndReadFastqImpl> ifp(is1, is2);
                    load_map(ifp, counters, read_map, start, length);
                }
            }
            if(vm.count("singleend-input")) {
                auto read_files = vm["singleend-input"].as<std::vector<std::string> >();
                for (auto file : read_files) {
                    std::ifstream read1(file, std::ifstream::in);
                    InputReader<SingleEndRead, SingleEndReadFastqImpl> ifs(read1);
                    load_map(ifs, counters, read_map, start, length);
                }
            }
            
            if(vm.count("tab-input")) {
                auto read_files = vm["tab-input"].as<std::vector<std::string> > ();
                for (auto file : read_files) {
                    std::ifstream tabReads(file, std::ifstream::in);
                    InputReader<ReadBase, TabReadImpl> ift(tabReads);
                    load_map(ift, counters, read_map, start, length);
                }
            }
            
            if (vm.count("interleaved-input")) {
                auto read_files = vm["interleaved-input"].as<std::vector<std::string > >();
                for (auto file : read_files) {
                    std::ifstream interReads(file, std::ifstream::in);
                    InputReader<PairedEndRead, InterReadImpl> ifp(interReads);
                    load_map(ifp, counters, read_map, start, length);
                }
            }

            if (fastq_out || (! std_out && ! tab_out) ) {
                for (auto& outfile: default_outfiles) {
                    outfile = prefix + outfile + ".fastq";
                }
                
                if (gzip_out) {
                    bi::stream<bi::file_descriptor_sink> out1{fileno(popen(("gzip > " + default_outfiles[0] + ".gz").c_str(), "w")), bi::close_handle};
                    bi::stream<bi::file_descriptor_sink> out2{fileno(popen(("gzip > " + default_outfiles[1] + ".gz").c_str(), "w")), bi::close_handle};
                    bi::stream<bi::file_descriptor_sink> out3{fileno(popen(("gzip > " + default_outfiles[2] + ".gz").c_str(), "w")), bi::close_handle};
                    output_read_map_fastq(read_map, out1, out2, out3);
                } else {
                    // note: mapped file is faster but uses more memory
                    std::ofstream out1(default_outfiles[0], std::ofstream::out);
                    std::ofstream out2(default_outfiles[1], std::ofstream::out);
                    std::ofstream out3(default_outfiles[2], std::ofstream::out);
                    //bi::stream<bi::mapped_file_sink> out1{default_outfiles[0].c_str()};
                    //bi::stream<bi::mapped_file_sink> out2{default_outfiles[1].c_str()};
                    //bi::stream<bi::mapped_file_sink> out3{default_outfiles[2].c_str()};
                    output_read_map_fastq(read_map, out1, out2, out3);
                }
            } else if (interleaved_out)  {
                for (auto& outfile: default_outfiles) {
                    outfile = prefix + "INTER" + ".fastq";
                }

                if (gzip_out) {
                    bi::stream<bi::file_descriptor_sink> out1{fileno(popen(("gzip > " + default_outfiles[0] + ".gz").c_str(), "w")), bi::close_handle};
                    bi::stream<bi::file_descriptor_sink> out3{fileno(popen(("gzip > " + default_outfiles[2] + ".gz").c_str(), "w")), bi::close_handle};
                    output_read_map_inter(read_map, out1, out3);
                } else {
                    std::ofstream out1(default_outfiles[0], std::ofstream::out);
                    std::ofstream out3(default_outfiles[2], std::ofstream::out);
                    output_read_map_inter(read_map, out1, out3);
                }
            } else if (tab_out) {
                for (auto& outfile: default_outfiles) {
                    outfile = prefix + "tab" + ".tastq";
                }
                
                if (gzip_out) {
                    bi::stream<bi::file_descriptor_sink> out1{fileno(popen(("gzip > " + default_outfiles[0] + ".gz").c_str(), "w")), bi::close_handle};
                    output_read_map_tab(read_map, out1);
                } else {
                    std::ofstream out1(default_outfiles[0], std::ofstream::out);
                    output_read_map_tab(read_map, out1);
                }
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

    std::cerr << "TotalRecords:" << counters["TotalRecords"] << "\tReplaced:" << counters["Replaced"]
              << "\tKept:" << read_map.size() << "\tRemoved:" << counters["TotalRecords"] - read_map.size()
              << "\tHasN:" << counters["HasN"] << std::endl;
    return SUCCESS;

}
