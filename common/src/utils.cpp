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
    outStats << "Program" << '\t';
    if (end == -1 || !appendStats) {
        std::string header;
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

