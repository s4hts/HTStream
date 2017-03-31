#include "utils.h"

void setupCounter(Counter &c) {

    c["TotalReadsInput"] = 0;
    c["PE_In"] = 0;
    c["PE_Out"] = 0;
    c["SE_In"] = 0;
    c["SE_Out"] = 0;
    c["R1_Length"] = 0;
    c["R2_Length"] = 0;
    c["SE_Length"] = 0;
    c["R1_Discarded"] = 0;
    c["R2_Discarded"] = 0;
    c["SE_Discarded"] = 0;
    c["R1_Left_Trim"] = 0;
    c["R1_Right_Trim"] = 0;
    c["R2_Left_Trim"] = 0;
    c["R2_Right_Trim"] = 0;
    c["SE_Right_Trim"] = 0;
    c["SE_Left_Trim"] = 0;
    c["Overlap_BPs"] = 0;
    c["Sins"] = 0;
    c["Lins"] = 0;
    c["Replaced"] = 0;
    c["Ignored"] = 0;

}

void write_stats(const std::string &statsFile, const bool &appendStats, const Counter &c, const std::string &program_name) {

    std::ifstream testEnd(statsFile);
    int end = testEnd.peek();
    testEnd.close();

    std::ofstream outStats;

    if (appendStats) {
        outStats.open(statsFile, std::ofstream::out | std::ofstream::app); //overwritte
    } else {
        outStats.open(statsFile, std::ofstream::out); //overwritte
    }

    if (end == -1 || !appendStats) {
        std::string header("Program\t");
        for (const auto name : c) {
            header += name.first + '\t';
        }
        header.replace(header.length()-1, 1, "\n");
        outStats << header;
    }

    std::string info;

    outStats << program_name << '\t';
    for (const auto name : c) {
        info += std::to_string(name.second) + '\t';
    }
    info.replace(info.length()-1, 1, "\n");
    outStats << info;
}

void outputWriters(std::shared_ptr<OutputWriter> &pe, std::shared_ptr<OutputWriter> &se, bool fastq_out, bool tab_out, bool interleaved_out, bool unmapped_out,  bool force, bool gzip_out, bool std_out, std::string &prefix) {

    std::vector<std::string> default_outfiles = {"PE1", "PE2", "SE"};

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
            ("stdin-input,S", "STDIN input <MUST BE TAB DELIMITED INPUT>")
            ("gzip-output,g", po::bool_switch()->default_value(false),  "Output gzipped")
            ("interleaved-output,i", po::bool_switch()->default_value(false),     "Output to interleaved")
            ("fastq-output,f", po::bool_switch()->default_value(true), "Fastq format output")
            ("force,F", po::bool_switch()->default_value(false),         "Forces overwrite of files")
            ("tab-output,t", po::bool_switch()->default_value(false),   "Tab-delimited output")
            ("unmapped-output,u", po::bool_switch()->default_value(false),   "Unmapped sam output")
            ("to-stdout,O", po::bool_switch()->default_value(false),    "Prints to STDOUT in Tab Delimited")
            ("stats-file,L", po::value<std::string>()->default_value("stats.log") , "String for output stats file name")
            ("append-stats-file,A", po::bool_switch()->default_value(false),  "Append Stats file.")
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
