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

void outputWriters(std::shared_ptr<OutputWriter> &pe, std::shared_ptr<OutputWriter> &se, bool fastq_out, bool tab_out, bool interleaved_out, bool unmapped_out,  bool force, bool gzip_out, bool std_out, std::string &prefix) {

    std::vector<std::string> default_outfiles = {"_R1", "_R2", "_SE"};

    std::shared_ptr<HtsOfstream> out_1 = nullptr;
    std::shared_ptr<HtsOfstream> out_2 = nullptr;
    std::shared_ptr<HtsOfstream> out_3 = nullptr;
    
    if (interleaved_out)  {
        for (auto& outfile: default_outfiles) {
            outfile = prefix + "INTER" + ".fastq";
        }

        out_1= std::make_shared<HtsOfstream>(default_outfiles[0], force, gzip_out, false);
        out_3= std::make_shared<HtsOfstream>(default_outfiles[1], force, gzip_out, false);

        pe= std::make_shared<PairedEndReadOutInter>(out_1);
        se= std::make_shared<SingleEndReadOutFastq>(out_3);
    } else if (unmapped_out) {
        for (auto& outfile: default_outfiles) {
            outfile = prefix + ".sam";
        }
        out_1= std::make_shared<HtsOfstream>(default_outfiles[0], force, gzip_out, std_out);

        pe= std::make_shared<ReadBaseOutUnmapped>(out_1);
        se= std::make_shared<ReadBaseOutUnmapped>(out_1);

    } else if (tab_out || std_out) {
        for (auto& outfile: default_outfiles) {
            outfile = prefix + "tab" + ".tastq";
        }
        out_1= std::make_shared<HtsOfstream>(default_outfiles[0], force, gzip_out, std_out);

        pe= std::make_shared<ReadBaseOutTab>(out_1);
        se= std::make_shared<ReadBaseOutTab>(out_1);
    } else if (fastq_out || (! std_out && ! tab_out) ) {
        for (auto& outfile: default_outfiles) {
            outfile = prefix + outfile + ".fastq";
        }

        out_1= std::make_shared<HtsOfstream>(default_outfiles[0], force, gzip_out, false);
        out_2= std::make_shared<HtsOfstream>(default_outfiles[1], force, gzip_out, false);
        out_3= std::make_shared<HtsOfstream>(default_outfiles[2], force, gzip_out, false);

        pe= std::make_shared<PairedEndReadOutFastq>(out_1, out_2);
        se= std::make_shared<SingleEndReadOutFastq>(out_3);
    }
}

void setDefaultParams(po::options_description &desc, std::string program_name) {

    desc.add_options()
            ("version,v", "Version print")
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
            ("std-input,S", "STDIN input <MUST BE TAB DELIMITED INPUT>")
            ("gzip-output,g", po::bool_switch()->default_value(false),  "Output gzipped")
            ("interleaved-output,i", po::bool_switch()->default_value(false),     "Output to interleaved")
            ("fastq-output,f", po::bool_switch()->default_value(true), "Fastq format output")
            ("force,F", po::bool_switch()->default_value(false),         "Forces overwrite of files")
            ("tab-output,t", po::bool_switch()->default_value(false),   "Tab-delimited output")
            ("unmapped-output,u", po::bool_switch()->default_value(false),   "Unmapped sam output")
            ("to-stdout,O", po::bool_switch()->default_value(false),    "Prints to STDOUT in Tab Delimited")
            ("stats-file,L", po::value<std::string>()->default_value("stats.log") , "String for output stats file name")
            ("append-stats-file,A", po::bool_switch()->default_value(false),  "Append Stats file.")
            ("notes,N", po::value<std::string>()->default_value(""),  "Notes for the JSON.")
            ("prefix,p", po::value<std::string>()->default_value(program_name),
                                           "Prefix for outputted files") 


            ("help,h",                     "Prints help.");

}

void setDefaultParamsTrim(po::options_description &desc) {
    desc.add_options()
        ("no-left,l", po::bool_switch()->default_value(false),    "Turns of trimming of the left side of the read")
        ("no-right,r", po::bool_switch()->default_value(false),    "Turns of trimming of the right side of the read");
 
}

void setDefaultParamsCutting(po::options_description &desc) {

    desc.add_options()
            ("no-orphans,n", po::bool_switch()->default_value(false), "SE reads will be NOT be written out")
            ("stranded,s", po::bool_switch()->default_value(false),    "If R1 is orphaned, R2 is RC (for stranded RNA)")
            ("min-length,m", po::value<size_t>()->default_value(50),    "Min length for acceptable outputted read");

}

void version_or_help(std::string program_name, po::options_description &desc, po::variables_map vm) {

    int SUCCESS = 0;
    if (vm.count("version")) {
        std::cout << program_name << std::endl;
        std::cout << "Version " << VERSION << std::endl;
        exit(SUCCESS); //success
    } else if ( vm.count("help")  || vm.size() == 0) {
        std::cout << program_name << std::endl
                  << "Version " << VERSION << std::endl
                  << desc << std::endl;
        exit(SUCCESS); //success
    } 
}

char rc (const char bp) {
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
            throw std::runtime_error("Unknown base alled in rc");
    }
}
