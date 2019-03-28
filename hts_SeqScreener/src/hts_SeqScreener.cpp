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
#include <cstdlib>

#include "hts_SeqScreener.h"

const std::string phixSeq_True = "GAGTTTTATCGCTTCCATGACGCAGAAGTTAACACTTTCGGATATTTCTGATGAGTCGAAAAATTATCTTGATAAAGCAGGAATTACTACTGCTTGTTTACGAATTAAATCGAAGTGGACTGCTGGCGGAAAATGAGAAAATTCGACCTATCCTTGCGCAGCTCGAGAAGCTCTTACTTTGCGACCTTTCGCCATCAACTAACGATTCTGTCAAAAACTGACGCGTTGGATGAGGAGAAGTGGCTTAATATGCTTGGCACGTTCGTCAAGGACTGGTTTAGATATGAGTCACATTTTGTTCATGGTAGAGATTCTCTTGTTGACATTTTAAAAGAGCGTGGATTACTATCTGAGTCCGATGCTGTTCAACCACTAATAGGTAAGAAATCATGAGTCAAGTTACTGAACAATCCGTACGTTTCCAGACCGCTTTGGCCTCTATTAAGCTCATTCAGGCTTCTGCCGTTTTGGATTTAACCGAAGATGATTTCGATTTTCTGACGAGTAACAAAGTTTGGATTGCTACTGACCGCTCTCGTGCTCGTCGCTGCGTTGAGGCTTGCGTTTATGGTACGCTGGACTTTGTGGGATACCCTCGCTTTCCTGCTCCTGTTGAGTTTATTGCTGCCGTCATTGCTTATTATGTTCATCCCGTCAACATTCAAACGGCCTGTCTCATCATGGAAGGCGCTGAATTTACGGAAAACATTATTAATGGCGTCGAGCGTCCGGTTAAAGCCGCTGAATTGTTCGCGTTTACCTTGCGTGTACGCGCAGGAAACACTGACGTTCTTACTGACGCAGAAGAAAACGTGCGTCAAAAATTACGTGCGGAAGGAGTGATGTAATGTCTAAAGGTAAAAAACGTTCTGGCGCTCGCCCTGGTCGTCCGCAGCCGTTGCGAGGTACTAAAGGCAAGCGTAAAGGCGCTCGTCTTTGGTATGTAGGTGGTCAACAATTTTAATTGCAGGGGCTTCGGCCCCTTACTTGAGGATAAATTATGTCTAATATTCAAACTGGCGCCGAGCGTATGCCGCATGACCTTTCCCATCTTGGCTTCCTTGCTGGTCAGATTGGTCGTCTTATTACCATTTCAACTACTCCGGTTATCGCTGGCGACTCCTTCGAGATGGACGCCGTTGGCGCTCTCCGTCTTTCTCCATTGCGTCGTGGCCTTGCTATTGACTCTACTGTAGACATTTTTACTTTTTATGTCCCTCATCGTCACGTTTATGGTGAACAGTGGATTAAGTTCATGAAGGATGGTGTTAATGCCACTCCTCTCCCGACTGTTAACACTACTGGTTATATTGACCATGCCGCTTTTCTTGGCACGATTAACCCTGATACCAATAAAATCCCTAAGCATTTGTTTCAGGGTTATTTGAATATCTATAACAACTATTTTAAAGCGCCGTGGATGCCTGACCGTACCGAGGCTAACCCTAATGAGCTTAATCAAGATGATGCTCGTTATGGTTTCCGTTGCTGCCATCTCAAAAACATTTGGACTGCTCCGCTTCCTCCTGAGACTGAGCTTTCTCGCCAAATGACGACTTCTACCACATCTATTGACATTATGGGTCTGCAAGCTGCTTATGCTAATTTGCATACTGACCAAGAACGTGATTACTTCATGCAGCGTTACCATGATGTTATTTCTTCATTTGGAGGTAAAACCTCTTATGACGCTGACAACCGTCCTTTACTTGTCATGCGCTCTAATCTCTGGGCATCTGGCTATGATGTTGATGGAACTGACCAAACGTCGTTAGGCCAGTTTTCTGGTCGTGTTCAACAGACCTATAAACATTCTGTGCCGCGTTTCTTTGTTCCTGAGCATGGCACTATGTTTACTCTTGCGCTTGTTCGTTTTCCGCCTACTGCGACTAAAGAGATTCAGTACCTTAACGCTAAAGGTGCTTTGACTTATACCGATATTGCTGGCGACCCTGTTTTGTATGGCAACTTGCCGCCGCGTGAAATTTCTATGAAGGATGTTTTCCGTTCTGGTGATTCGTCTAAGAAGTTTAAGATTGCTGAGGGTCAGTGGTATCGTTATGCGCCTTCGTATGTTTCTCCTGCTTATCACCTTCTTGAAGGCTTCCCATTCATTCAGGAACCGCCTTCTGGTGATTTGCAAGAACGCGTACTTATTCGCCACCATGATTATGACCAGTGTTTCCAGTCCGTTCAGTTGTTGCAGTGGAATAGTCAGGTTAAATTTAATGTGACCGTTTATCGCAATCTGCCGACCACTCGCGATTCAATCATGACTTCGTGATAAAAGATTGAGTGTGAGGTTATAACGCCGAAGCGGTAAAAATTTTAATTTTTGCCGCTGAGGGGTTGACCAAGCGAAGCGCGGTAGGTTTTCTGCTTAGGAGTTTAATCATGTTTCAGACTTTTATTTCTCGCCATAATTCAAACTTTTTTTCTGATAAGCTGGTTCTCACTTCTGTTACTCCAGCTTCTTCGGCACCTGTTTTACAGACACCTAAAGCTACATCGTCAACGTTATATTTTGATAGTTTGACGGTTAATGCTGGTAATGGTGGTTTTCTTCATTGCATTCAGATGGATACATCTGTCAACGCCGCTAATCAGGTTGTTTCTGTTGGTGCTGATATTGCTTTTGATGCCGACCCTAAATTTTTTGCCTGTTTGGTTCGCTTTGAGTCTTCTTCGGTTCCGACTACCCTCCCGACTGCCTATGATGTTTATCCTTTGAATGGTCGCCATGATGGTGGTTATTATACCGTCAAGGACTGTGTGACTATTGACGTCCTTCCCCGTACGCCGGGCAATAACGTTTATGTTGGTTTCATGGTTTGGTCTAACTTTACCGCTACTAAATGCCGCGGATTGGTTTCGCTGAATCAGGTTATTAAAGAGATTATTTGTCTCCAGCCACTTAAGTGAGGTGATTTATGTTTGGTGCTATTGCTGGCGGTATTGCTTCTGCTCTTGCTGGTGGCGCCATGTCTAAATTGTTTGGAGGCGGTCAAAAAGCCGCCTCCGGTGGCATTCAAGGTGATGTGCTTGCTACCGATAACAATACTGTAGGCATGGGTGATGCTGGTATTAAATCTGCCATTCAAGGCTCTAATGTTCCTAACCCTGATGAGGCCGCCCCTAGTTTTGTTTCTGGTGCTATGGCTAAAGCTGGTAAAGGACTTCTTGAAGGTACGTTGCAGGCTGGCACTTCTGCCGTTTCTGATAAGTTGCTTGATTTGGTTGGACTTGGTGGCAAGTCTGCCGCTGATAAAGGAAAGGATACTCGTGATTATCTTGCTGCTGCATTTCCTGAGCTTAATGCTTGGGAGCGTGCTGGTGCTGATGCTTCCTCTGCTGGTATGGTTGACGCCGGATTTGAGAATCAAAAAGAGCTTACTAAAATGCAACTGGACAATCAGAAAGAGATTGCCGAGATGCAAAATGAGACTCAAAAAGAGATTGCTGGCATTCAGTCGGCGACTTCACGCCAGAATACGAAAGACCAGGTATATGCACAAAATGAGATGCTTGCTTATCAACAGAAGGAGTCTACTGCTCGCGTTGCGTCTATTATGGAAAACACCAATCTTTCCAAGCAACAGCAGGTTTCCGAGATTATGCGCCAAATGCTTACTCAAGCTCAAACGGCTGGTCAGTATTTTACCAATGACCAAATCAAAGAAATGACTCGCAAGGTTAGTGCTGAGGTTGACTTAGTTCATCAGCAAACGCAGAATCAGCGGTATGGCTCTTCTCATATTGGCGCTACTGCAAAGGATATTTCTAATGTCGTCACTGATGCTGCTTCTGGTGTGGTTGATATTTTTCATGGTATTGATAAAGCTGTTGCCGATACTTGGAACAATTTCTGGAAAGACGGTAAAGCTGATGGTATTGGCTCTAATTTGTCTAGGAAATAACCGTCAGGATTGACACCCTCCCAATTGTATGTTTTCATGCCTCCAAATCTTGGAGGCTTTTTTATGGTTCGTTCTTATTACCCTTCTGAATGTCACGCTGATTATTTTGACTTTGAGCGTATCGAGGCTCTTAAACCTGCTATTGAGGCTTGTGGCATTTCTACTCTTTCTCAATCCCCAATGCTTGGCTTCCATAAGCAGATGGATAACCGCATCAAGCTCTTGGAAGAGATTCTGTCTTTTCGTATGCAGGGCGTTGAGTTCGATAATGGTGATATGTATGTTGACGGCCATAAGGCTGCTTCTGACGTTCGTGATGAGTTTGTATCTGTTACTGAGAAGTTAATGGATGAATTGGCACAATGCTACAATGTGCTCCCCCAACTTGATATTAATAACACTATAGACCACCGCCCCGAAGGGGACGAAAAATGGTTTTTAGAGAACGAGAAGACGGTTACGCAGTTTTGCCGCAAGCTGGCTGCTGAACGCCCTCTTAAGGATATTCGCGATGAGTATAATTACCCCAAAAAGAAAGGTATTAAGGATGAGTGTTCAAGATTGCTGGAGGCCTCCACTATGAAATCGCGTAGAGGCTTTGCTATTCAGCGTTTGATGAATGCAATGCGACAGGCTCATGCTGATGGTTGGTTTATCGTTTTTGACACTCTCACGTTGGCTGACGACCGATTAGAGGCGTTTTATGATAATCCCAATGCTTTGCGTGACTATTTTCGTGATATTGGTCGTATGGTTCTTGCTGCCGAGGGTCGCAAGGCTAATGATTCACACGCCGACTGCTATCAGTATTTTTGTGTGCCTGAGTATGGTACAGCTAATGGCCGTCTTCATTTCCATGCGGTGCACTTTATGCGGACACTTCCTACAGGTAGCGTTGACCCTAATTTTGGTCGTCGGGTACGCAATCGCCGCCAGTTAAATAGCTTGCAAAATACGTGGCCTTATGGTTACAGTATGCCCATCGCAGTTCGCTACACGCAGGACGCTTTTTCACGTTCTGGTTGGTTGTGGCCTGTTGATGCTAAAGGTGAGCCGCTTAAAGCTACCAGTTATATGGCTGTTGGTTTCTATGTGGCTAAATACGTTAACAAAAAGTCAGATATGGACCTTGCTGCTAAAGGTCTAGGAGCTAAAGAATGGAACAACTCACTAAAAACCAAGCTGTCGCTACTTCCCAAGAAGCTGTTCAGAATCAGAATGAGCCGCAACTTCGGGATGAAAATGCTCACAATGACAAATCTGTCCACGGAGTGCTTAATCCAACTTACCAAGCTGGGTTACGACGCGACGCCGTTCAACCAGATATTGAAGCAGAACGCAAAAAGAGAGATGAGATTGAGGCTGGGAAAAGTTACTGTAGCCGACGTTTTGGCGGCGCAACCTGTGACGACAAATCTGCTCAAATTTATGCGCGCTTCGATAAAAATGATTGGCGTATCCAACCTGCA";

namespace
{
    const size_t SUCCESS = 0;
    const size_t ERROR_IN_COMMAND_LINE = 1;
    const size_t ERROR_UNHANDLED_EXCEPTION = 2;
    const size_t ERROR_NOLOOKUP_EXCEPTION = 3;

} // namespace

namespace bi = boost::iostreams;

int main(int argc, char** argv)
{
    const std::string program_name = "hts_SeqScreener";
    std::string app_description =
                       "hts_SeqScreener identifies and removes any reads which appear to have originated\n";
    app_description += "  from a contaminant DNA source. Because bacteriophage Phi-X is common spiked\n";
    app_description += "  into Illumina runs for QC purposes, sequences originating from Phi-X are removed\n";
    app_description += "  by default. If other contaminants are suspected their sequence can be supplied\n";
    app_description += "  as a fasta file <seq>, however the algorithm has been tuned for short contaminant\n";
    app_description += "  sequences, and may not work well with sequences significantly longer than Phi-X (5Kb).\n";

    try
    {
        /** Define and parse the program options
         */
        namespace po = boost::program_options;
        po::options_description standard = setStandardOptions();
            // version|v ; help|h ; notes|N ; stats-file|L ; append-stats-file|A
        po::options_description input = setInputOptions();
            // read1-input|1 ; read2-input|2 ; singleend-input|U
            // tab-input|T ; interleaved-input|I
        po::options_description output = setOutputOptions(program_name);
          // force|F ; uncompressed|u ; fastq-output|f
          // tab-output|t ; interleaved-output|i ; unmapped-output|z

        po::options_description desc("Application Specific Options");

        desc.add_options()
            ("seq,s", po::value<std::string>(), "Please supply a fasta file - default - Phix Sequence - default https://www.ncbi.nlm.nih.gov/nuccore/9626372")
            ("check-read-2,C", po::bool_switch()->default_value(false), "Check R2 as well as R1 (pe)")
            ("kmer,k", po::value<size_t>()->default_value(12)->notifier(boost::bind(&check_range<size_t>, "kmer", _1, 5, 256)), "Kmer size of the lookup table (min 5, max 256)")
            ("percentage-hits,x", po::value<double>()->default_value(.25)->notifier(boost::bind(&check_range<double>, "percentage-hits", _1, 0.0, 1.0)), "Proportion of kmer percentage-hits to sequence need to happen to discard (min 0.0, max 1.0)")
            ("inverse,n", po::bool_switch()->default_value(false), "Output reads that are ABOVE the kmer hit threshold")
            ("record,r", po::bool_switch()->default_value(false), "Only record the reads that pass the kmer hit threshold, output all reads");

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

            std::string statsFile(vm["stats-file"].as<std::string>());
            SeqScreenerCounters counters(statsFile, vm["force"].as<bool>(), vm["append-stats-file"].as<bool>(), program_name, vm["notes"].as<std::string>());

            //sets read information
            //Phix isn't set to default since it makes help a PITA to read
            //sets kmer lookup arrays
            kmerSet lookup;
            uint64_t screen_len;
            std::string lookup_file;
            if (vm.count("seq")) {
                lookup_file = vm["seq"].as<std::string>();
                bi::stream <bi::file_descriptor_source> fa{check_open_r(lookup_file), bi::close_handle};
                InputReader<SingleEndRead, FastaReadImpl> faReader(fa);
                screen_len = setLookup_fasta(lookup, faReader, vm["kmer"].as<size_t>());
            } else {
                Read readSeq;
                lookup_file = "PhiX";
                readSeq = Read(phixSeq_True, "", "");
                screen_len = readSeq.getLength();
                setLookup_read(lookup, readSeq, vm["kmer"].as<size_t>());
            }
            if (lookup.size() == 0){
                std::cerr << "\n\tException lookup table contains no kmers" << std::endl;
                return ERROR_NOLOOKUP_EXCEPTION;
            }

            counters.set_screeninfo(lookup_file, screen_len, lookup.size());

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
                    helper_discard(ifp, pe, se, counters, lookup, vm["percentage-hits"].as<double>(), vm["check-read-2"].as<bool>(),vm["kmer"].as<size_t>(), vm["inverse"].as<bool>(), vm["record"].as<bool>());
                }
                if(vm.count("singleend-input")) { // can have paired-end reads and/or single-end reads
                    auto read_files = vm["singleend-input"].as<std::vector<std::string> >();
                    for (auto file : read_files) {
                        bi::stream<bi::file_descriptor_source> sef{ check_open_r(file), bi::close_handle};
                        InputReader<SingleEndRead, SingleEndReadFastqImpl> ifs(sef);
                        helper_discard(ifs, pe, se, counters, lookup, vm["percentage-hits"].as<double>(), vm["check-read-2"].as<bool>(),vm["kmer"].as<size_t>(), vm["inverse"].as<bool>(), vm["record"].as<bool>());
                    }
                }
            } else if (vm.count("interleaved-input")) {
                auto read_files = vm["interleaved-input"].as<std::vector<std::string > >();
                for (auto file : read_files) {
                    bi::stream<bi::file_descriptor_source> inter{ check_open_r(file), bi::close_handle};
                    InputReader<PairedEndRead, InterReadImpl> ifp(inter);
                    helper_discard(ifp, pe, se, counters, lookup, vm["percentage-hits"].as<double>(), vm["check-read-2"].as<bool>(),vm["kmer"].as<size_t>(), vm["inverse"].as<bool>(), vm["record"].as<bool>());
                }
                if(vm.count("singleend-input")) { // can have interleaved paired-end reads and/or single-end reads
                    auto read_files = vm["singleend-input"].as<std::vector<std::string> >();
                    for (auto file : read_files) {
                        bi::stream<bi::file_descriptor_source> sef{ check_open_r(file), bi::close_handle};
                        InputReader<SingleEndRead, SingleEndReadFastqImpl> ifs(sef);
                        helper_discard(ifs, pe, se, counters, lookup, vm["percentage-hits"].as<double>(), vm["check-read-2"].as<bool>(),vm["kmer"].as<size_t>(), vm["inverse"].as<bool>(), vm["record"].as<bool>());
                    }
                }
            } else if(vm.count("singleend-input")) { // can also just have single-end reads
                auto read_files = vm["singleend-input"].as<std::vector<std::string> >();
                for (auto file : read_files) {
                    bi::stream<bi::file_descriptor_source> sef{ check_open_r(file), bi::close_handle};
                    InputReader<SingleEndRead, SingleEndReadFastqImpl> ifs(sef);
                    helper_discard(ifs, pe, se, counters, lookup, vm["percentage-hits"].as<double>(), vm["check-read-2"].as<bool>(),vm["kmer"].as<size_t>(), vm["inverse"].as<bool>(), vm["record"].as<bool>());
                }
            } else if(vm.count("tab-input")) {
                auto read_files = vm["tab-input"].as<std::vector<std::string> > ();
                for (auto file : read_files) {
                    bi::stream<bi::file_descriptor_source> tabin{ check_open_r(file), bi::close_handle};
                    InputReader<ReadBase, TabReadImpl> ift(tabin);
                    helper_discard(ift, pe, se, counters, lookup, vm["percentage-hits"].as<double>(), vm["check-read-2"].as<bool>(),vm["kmer"].as<size_t>(), vm["inverse"].as<bool>(), vm["record"].as<bool>());
                }
            }

            if (vm.count("interleaved-input")) {
                auto read_files = vm["interleaved-input"].as<std::vector<std::string > >();
                for (auto file : read_files) {
                    bi::stream<bi::file_descriptor_source> inter{ check_open_r(file), bi::close_handle};
                    InputReader<PairedEndRead, InterReadImpl> ifp(inter);
                    helper_discard(ifp, pe, se, counters, lookup, vm["percentage-hits"].as<double>(), vm["check-read-2"].as<bool>(),vm["kmer"].as<size_t>(), vm["inverse"].as<bool>(), vm["record"].as<bool>());
                }
            }

            if (!isatty(fileno(stdin))) {
                bi::stream<bi::file_descriptor_source> tabin {fileno(stdin), bi::close_handle};
                InputReader<ReadBase, TabReadImpl> ift(tabin);
                helper_discard(ift, pe, se, counters, lookup, vm["percentage-hits"].as<double>(), vm["check-read-2"].as<bool>(),vm["kmer"].as<size_t>(), vm["inverse"].as<bool>(), vm["record"].as<bool>());
            } else {
              std::cerr << "ERROR: " << "Input file type absent from command line" << std::endl << std::endl;
              version_or_help(program_name, app_description, cmdline_options, vm, true);
              exit(ERROR_IN_COMMAND_LINE); //success
            }
            counters.write_out();
        }
        catch(po::error& e)
        {
            std::cerr << "ERROR: " << e.what() << std::endl << std::endl;
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
