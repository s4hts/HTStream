#ifndef COUNTERS_H
#define COUNTERS_H

#include <unordered_map>
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include "read.h"
#include "typedefs.h"
#include <unistd.h>
class Counters {
    
public:
    Counter c;
   
    void Common() {
        c["TotalReadsInput"] = 0;
        c["PE_In"] = 0;
        c["SE_In"] = 0;
        c["TotalReadsOutput"] = 0;
        c["PE_Out"] = 0;
        c["SE_Out"] = 0;
    }

    Counters() {
        Common();
    }
    
    virtual void input(const PairedEndRead &read) {
        ++c["TotalReadsInput"];
        ++c["PE_In"]; 
    }



    virtual void input(const SingleEndRead &read) {
        ++c["TotalReadsInput"];
        ++c["SE_In"]; 
    }

    virtual void input(const ReadBase &read) {
        const PairedEndRead *per = dynamic_cast<const PairedEndRead *>(&read);
        if (per) {
            input(*per);
        } else {
            const SingleEndRead *ser = dynamic_cast<const SingleEndRead *>(&read);
            if (ser) {
                input(*ser); 
            } else {
                throw std::runtime_error("In utils.h output: read type not valid");
            }
        } 
    }

    virtual void output(PairedEndRead &read, bool no_orhpans = false) {
        ++c["TotalReadsOutput"];
        ++c["PE_Out"];
    }


    virtual void output(SingleEndRead &read) {
        ++c["TotalReadsOutput"];
        ++c["SE_Out"];
    }

    virtual void output(ReadBase &read) {
        PairedEndRead *per = dynamic_cast<PairedEndRead *>(&read);
        if (per) {
            output(*per);
        } else {
            SingleEndRead *ser = dynamic_cast<SingleEndRead *>(&read);
            if (ser) {
                output(*ser); 
            } else {
                throw std::runtime_error("In utils.h output: read type not valid");
            }
        } 
    }

    virtual void write_out(std::string &statsFile, bool appendStats, std::string program_name, std::string notes) {
        
        std::ifstream testEnd(statsFile);
        int end = testEnd.peek();
        testEnd.close();
        
        std::fstream outStats;
        bool first = true;
        
        if (appendStats && end != -1) {
            //outStats.open(statsFile, std::ofstream::out | std::ofstream::app); //overwritte
            outStats.open(statsFile, std::ios::out|std::ios::in); //overwritte
            outStats.seekp(-1, std::ios::end );
            outStats << "\n,\"" << program_name << "_" << getpid()  << "\": {\n";
        } else {
            //outStats.open(statsFile, std::ofstream::out); //overwritte
            outStats.open(statsFile, std::ios::out); //overwritt
            outStats << "{\n \"" << program_name << "_" << getpid() <<  "\": {\n";
        }
        outStats << "\"Notes\" : \"" << notes << "\",\n";
        for (const auto name : c) {
            if (first) {
                first = false;
            } else {
                outStats << ",\n"; //make sure json format is kept
            }
            outStats << "\"" << name.first << "\" : " << name.second; //it will get the comma in conditionals tatement about
        }
        outStats << "\n}";
        outStats << "\n}";
        outStats.flush();
         
    }
};


class TrimmingCounters : public Counters {

public: 
    TrimmingCounters() {
        Common();
        c["R1_Left_Trim"] = 0;
        c["R1_Right_Trim"] = 0;
        c["R2_Left_Trim"] = 0;
        c["R2_Right_Trim"] = 0;
        c["SE_Right_Trim"] = 0;
        c["SE_Left_Trim"] = 0;

        c["R1_Discarded"] = 0;
        c["R2_Discarded"] = 0;
        c["SE_Discarded"] = 0;

    }

    void R1_stats(Read &one) {
        c["R1_Length"] += one.getLengthTrue();
        c["R1_Left_Trim"] += one.getLTrim();
        c["R1_Right_Trim"] += one.getRTrim();
    }

    void R2_stats(Read &two) {
        c["R2_Length"] += two.getLengthTrue();
        c["R2_Left_Trim"] += two.getLTrim();
        c["R2_Right_Trim"] += two.getRTrim();
    }

    void SE_stats(Read &se) {
        c["SE_Length"] += se.getLengthTrue();
        c["SE_Left_Trim"] += se.getLTrim();
        c["SE_Right_Trim"] += se.getRTrim();
    }

    void output(PairedEndRead &per, bool no_orphans = false) {
        Read &one = per.non_const_read_one();
        Read &two = per.non_const_read_one();
        if (!one.getDiscard() && !two.getDiscard()) {
            ++c["PE_Out"];
            R1_stats(one);
            R2_stats(two);
        } else if (!one.getDiscard() && !no_orphans) { //if stranded RC
            ++c["SE_Out"];
            ++c["R2_Discarded"];
            R1_stats(one);
        } else if (!two.getDiscard() && !no_orphans) { // Will never be RC
            ++c["SE_Out"];
            ++c["R1_Discarded"];
            R2_stats(two);
        } else {
            ++c["R1_Discarded"];
            ++c["R2_Discarded"];
        }
    }

    void output(SingleEndRead &ser) {
        Read &one = ser.non_const_read_one();
        if (!one.getDiscard()) {
            ++c["SE_Out"];
            SE_stats(one);
        } else {
            ++c["SE_Discarded"];
        }
    }
    
};

#endif
