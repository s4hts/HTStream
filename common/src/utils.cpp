#include "utils.h"

bool threshold_mismatches(std::string::const_iterator r1, std::string::const_iterator r2, size_t length, size_t max) {
    for (size_t i = 0; i < length; ++i) {
        if (*r1 != *r2) {
            --max;
        }
        if (!max) {
            return false;
        }
        ++r1;
        ++r2;
    }
    return true;
}


/* Create the quick lookup table
 * Multi map because a single kemr could appear multiple places */
seqLookup readOneMap(std::string seq1, const size_t kmer, const size_t kmerOffset) {

    seqLookup baseReadMap;
    std::string::iterator it;
    for ( it = seq1.begin() ; it <= seq1.end() - ( static_cast<long> ( kmer ) ) ; it += static_cast<long> ( kmerOffset )   ) {
        baseReadMap.insert(std::make_pair( std::string ( it, it+ static_cast<long> ( kmer )  ) , it - seq1.begin() ));
    }

    return baseReadMap;
}

void outputWriters(std::shared_ptr<OutputWriter> &pe, std::shared_ptr<OutputWriter> &se, bool fastq_out, bool tab_out, bool interleaved_out, bool unmapped_out,  bool force, bool gzip_out, std::string &prefix) {

    std::vector<std::string> default_outfiles = {"_R1", "_R2", "_SE", "_INTERLEAVED"};

    std::shared_ptr<HtsOfstream> out_1 = nullptr;
    std::shared_ptr<HtsOfstream> out_2 = nullptr;
    std::shared_ptr<HtsOfstream> out_3 = nullptr;

    if (interleaved_out) {
        for (auto& outfile: default_outfiles) {
            outfile = prefix + outfile + ".fastq";
        }

        out_1= std::make_shared<HtsOfstream>(default_outfiles[3], force, gzip_out, false);
        out_3= std::make_shared<HtsOfstream>(default_outfiles[2], force, gzip_out, false);

        pe= std::make_shared<PairedEndReadOutInter>(out_1);
        se= std::make_shared<SingleEndReadOutFastq>(out_3);
    } else if (unmapped_out) {
        for (auto& outfile: default_outfiles) {
            outfile = prefix + ".sam";
        }
        out_1= std::make_shared<HtsOfstream>(default_outfiles[0], force, gzip_out, false);

        pe= std::make_shared<ReadBaseOutUnmapped>(out_1);
        se= std::make_shared<ReadBaseOutUnmapped>(out_1);
    } else if (tab_out) {
        for (auto& outfile: default_outfiles) {
            outfile = prefix + ".tab6";
        }
        out_1= std::make_shared<HtsOfstream>(default_outfiles[0], force, gzip_out, false);

        pe= std::make_shared<ReadBaseOutTab>(out_1);
        se= std::make_shared<ReadBaseOutTab>(out_1);
    } else if (fastq_out) {
        for (auto& outfile: default_outfiles) {
            outfile = prefix + outfile + ".fastq";
        }
        out_1= std::make_shared<HtsOfstream>(default_outfiles[0], force, gzip_out, false);
        out_2= std::make_shared<HtsOfstream>(default_outfiles[1], force, gzip_out, false);
        out_3= std::make_shared<HtsOfstream>(default_outfiles[2], force, gzip_out, false);

        pe= std::make_shared<PairedEndReadOutFastq>(out_1, out_2);
        se= std::make_shared<SingleEndReadOutFastq>(out_3);
    } else {
        for (auto& outfile: default_outfiles) {
            outfile = prefix + ".tab6";
        }
        out_1= std::make_shared<HtsOfstream>(default_outfiles[0], force, gzip_out, true);

        pe= std::make_shared<ReadBaseOutTab>(out_1);
        se= std::make_shared<ReadBaseOutTab>(out_1);
    }
}

po::options_description setInputOptions(){

    po::options_description input("Input Options");
    input.add_options()
            //input options
            ("read1-input,1", po::value< std::vector<std::string> >()->multitoken(),
                                           "Read 1 paired end fastq input <space seperated for multiple files>")
            ("read2-input,2", po::value< std::vector<std::string> >()->multitoken(),
                                           "Read 2 paired end fastq input <space seperated for multiple files>")
            ("singleend-input,U", po::value< std::vector<std::string> >()->multitoken(),
                                           "Single end read fastq input <space seperated for multiple files>")
            ("tab-input,T", po::value< std::vector<std::string> >()->multitoken(),
                                           "Tab input <space seperated for multiple files>")
            ("interleaved-input,I", po::value< std::vector<std::string> >()->multitoken(),
                                           "Interleaved fastq input <space seperated for multiple files>")
            ("from-stdin,S", "STDIN input <MUST BE TAB DELIMITED INPUT>, DEFAULT AND DEPRECATED PARAMETER");
    return input;
}

po::options_description setOutputOptions(std::string program_name){
    po::options_description output("Output Options");
    output.add_options()
            //output options
            ("force,F", po::bool_switch()->default_value(false),         "Forces overwrite of files")
            ("prefix,p", po::value<std::string>()->default_value(program_name),
                                           "Prefix for output files")
            ("gzip-output,g", po::bool_switch()->default_value(false),  "Output gzipped files")
            ("fastq-output,f", po::bool_switch()->default_value(false), "Output to Fastq format <PE AND/OR SE files>")
            ("tab-output,t", po::bool_switch()->default_value(false),   "Output to tab-delimited file format")
            ("interleaved-output,i", po::bool_switch()->default_value(false),     "Output to interleaved fastq file <PE ONLY>")
            ("unmapped-output,u", po::bool_switch()->default_value(false),   "Output to unmapped sam file format")
            ("to-stdout,O", po::bool_switch()->default_value(false),    "Output to STDOUT in tab-delimited file format, DEFAULT AND DEPRECATED PARAMETER");
    return output;
}

po::options_description setStandardOptions(){
    po::options_description standard("Standard Options");
    standard.add_options()
            // version, help, notes
            ("version,v", "Version print")
            ("help,h",  "Prints help documentation")
            ("notes,N", po::value<std::string>()->default_value(""),  "Notes for the stats JSON")
            //stats file
            ("stats-file,L", po::value<std::string>()->default_value("stats.log") , "String for output stats file name")
            ("append-stats-file,A", po::bool_switch()->default_value(false),  "Append to stats file");
    return standard;
}

void setDefaultParamsTrim(po::options_description &desc) {
    desc.add_options()
            ("no-left,l", po::bool_switch()->default_value(false),    "Turns off trimming of the left side of the read")
            ("no-right,r", po::bool_switch()->default_value(false),    "Turns off trimming of the right side of the read");
}

void setDefaultParamsCutting(po::options_description &desc) {

    desc.add_options()
            ("no-orphans,n", po::bool_switch()->default_value(false), "Orphaned SE reads will NOT be written out")
            ("stranded,s", po::bool_switch()->default_value(false),    "If R1 is orphaned, R2 is output in RC (for stranded RNA)")
            ("min-length,m", po::value<size_t>()->default_value(1)->notifier(boost::bind(&check_range<size_t>, "min-length", _1, 1, 10000)), "Min length for acceptable output read (min 1, max 10000)");
}

void setDefaultParamsOverlapping(po::options_description &desc) {
    desc.add_options()
            ("kmer,k", po::value<size_t>()->default_value(8)->notifier(boost::bind(&check_range<size_t>, "kmer", _1, 5, 64)), "Kmer size of the lookup table for the longer read (min 5, max 64)")
            ("kmer-offset,r", po::value<size_t>()->default_value(1)->notifier(boost::bind(&check_range<size_t>, "kmer-offset", _1, 1, 10000)), "Offset of kmers. Offset of 1, would be perfect overlapping kmers. An offset of kmer would be non-overlapping kmers that are right next to each other. Must be greater than 0.")
            ("max-mismatch-errorDensity,e", po::value<double>()->default_value(.25)->notifier(boost::bind(&check_range<double>, "max-mismatch-errorDensity", _1, 0.0, 1.0)), "Max percent of mismatches allowed in overlapped section (min 0.0, max 1.0)")
            ("max-mismatch,x", po::value<size_t>()->default_value(100)->notifier(boost::bind(&check_range<size_t>, "max-mismatch", _1, 0, 10000)), "Max number of total mismatches allowed in overlapped section (min 0, max 10000)")
            ("check-lengths,c", po::value<size_t>()->default_value(20)->notifier(boost::bind(&check_range<size_t>, "check-lengths", _1, 5, 10000)), "Check lengths of the ends (min 5, max 10000)")
            ("min-overlap,o", po::value<size_t>()->default_value(8)->notifier(boost::bind(&check_range<size_t>, "min-length", _1, 5, 10000)), "Min overlap required to merge two reads (min 5, max 10000)");
}

void version_or_help(std::string program_name, std::string app_description, po::options_description &desc, po::variables_map vm, bool error) {

    std::string prolog="HTStream <https://github.com/ibest/HTStream> application: " + program_name;
    std::string epilog="Please report any issues, request for enhancement, or comments to <https://github.com/ibest/HTStream/issues>";
    int SUCCESS = 0;
    if (vm.count("version")) {
        std::cout << VERSION << std::endl;
        exit(SUCCESS); //success
    } else if ( vm.count("help")  || error || vm.size() == 0) {
        std::cout << prolog << std::endl;
        std::cout << "Version: " << VERSION << std::endl;
        std::cout << app_description << std::endl;
        std::cout << desc << std::endl;
        std::cout << std::endl << epilog << std::endl;
        exit(SUCCESS); //success
    } /* else if ( !vm.count("read1-input") & !vm.count("singleend-input") & !vm.count("tab-input") & !vm.count("interleaved-input") & !vm.count("from-stdin") ){
        std::cout << prolog << std::endl;
        std::cout << "Version: " << VERSION << std::endl;
        std::cout << app_description << std::endl;
        std::cout << desc << std::endl;
        std::cout << std::endl << epilog << std::endl;
        exit(SUCCESS); //success
    } */
}

char rc (const char &bp) {
    switch (bp) {
        case 'C':
            return 'G';
        case 'G':
            return 'C';
        case 'A':
            return 'T';
        case 'T':
            return 'A';
        case 'N':
            return 'N';
        default:
            throw std::runtime_error("Unknown base provided to rc");
    }
}
