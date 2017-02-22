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

    if (fastq_out || (! std_out && ! tab_out) ) {
        for (auto& outfile: default_outfiles) {
            outfile = prefix + outfile + ".fastq";
        }

        out_1.reset(new HtsOfstream(default_outfiles[0], force, gzip_out, false));
        out_2.reset(new HtsOfstream(default_outfiles[1], force, gzip_out, false));
        out_3.reset(new HtsOfstream(default_outfiles[2], force, gzip_out, false));

        pe.reset(new PairedEndReadOutFastq(out_1, out_2));
        se.reset(new SingleEndReadOutFastq(out_3));
    } else if (interleaved_out)  {
        for (auto& outfile: default_outfiles) {
            outfile = prefix + "INTER" + ".fastq";
        }

        out_1.reset(new HtsOfstream(default_outfiles[0], force, gzip_out, false));
        out_3.reset(new HtsOfstream(default_outfiles[1], force, gzip_out, false));

        pe.reset(new PairedEndReadOutInter(out_1));
        se.reset(new SingleEndReadOutFastq(out_3));
    } else if (tab_out || std_out) {
        for (auto& outfile: default_outfiles) {
            outfile = prefix + "tab" + ".tastq";
        }
        out_1.reset(new HtsOfstream(default_outfiles[0], force, gzip_out, std_out));

        pe.reset(new ReadBaseOutTab(out_1));
        se.reset(new ReadBaseOutTab(out_1));
    }

}

void setDefaultParams(boost::program_options::options_description &descp) {
    descp.add_options() ("test,t", "Test");
    /*namespace po = boost::program_options;
    po::options_description test;
    test.add_options()
        ("version,v", "Version Print");*/
/*
        desc.add_options()
        ("version,v",                  "Version print");*/
}

