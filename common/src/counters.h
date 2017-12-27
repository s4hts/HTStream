#ifndef COUNTERS_H
#define COUNTERS_H

#include <map>
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
        c["TotalFragmentsInput"] = 0;
        c["PE_In"] = 0;
        c["SE_In"] = 0;
        c["TotalFragmentsOutput"] = 0;
        c["PE_Out"] = 0;
        c["SE_Out"] = 0;
    }

    Counters() {
        Common();
    }
    
    virtual void input(SingleEndRead &read) {
        ++c["TotalFragmentsInput"];
        ++c["PE_In"]; 
    }

    virtual void input(PairedEndRead &read) {
        ++c["TotalFragmentsInput"];
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

    virtual void output(PairedEndRead &read, bool no_orphans = false) {
        ++c["TotalFragmentsOutput"];
        ++c["PE_Out"];
    }


    virtual void output(SingleEndRead &read) {
        ++c["TotalFragmentsOutput"];
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

    virtual void write_out(const std::string &statsFile, bool appendStats, std::string program_name, std::string notes) {
        
        std::ifstream testEnd(statsFile);
        int end = testEnd.peek();
        testEnd.close();
        
        std::fstream outStats;
        
        if (appendStats && end != -1) {
            //outStats.open(statsFile, std::ofstream::out | std::ofstream::app); //overwritte
            outStats.open(statsFile, std::ios::in | std::ios::out); //overwritte
            outStats.seekp(-6, std::ios::end );
            outStats << "  }, \"" << program_name << "_" << getpid()  << "\": {\n";
        } else {
            //outStats.open(statsFile, std::ofstream::out); //overwritte
            outStats.open(statsFile, std::ios::out | std::ios::trunc); //overwrite
            outStats << "{ \"" << program_name << "_" << getpid() <<  "\": {\n";
        }
        outStats << "    \"Notes\": \"" << notes << "\"";

        for (const auto &name : c) {
            outStats << ",\n"; //make sure json format is kept
            outStats << "    \"" << name.first << "\": " << name.second;
        }
        outStats << "\n  }\n}\n";
        outStats.flush();
    }
};


class OverlappingCounters : public Counters {

public:
    std::vector<unsigned long long int> insertLength;

    OverlappingCounters() {
        Common();
        c["sins"] = 0;
        c["mins"] = 0;
        c["lins"] = 0;
        c["SE_Discard"] = 0;
        c["PE_Discard"] = 0;
        c["Adapter_BpTrim"] = 0;
        insertLength.resize(1);
    }

    virtual void output(SingleEndRead &ser)  {
        if (ser.non_const_read_one().getDiscard()) {
            ++c["SE_Discard"];
        } else {
            ++c["TotalFragmentsOutput"];
            ++c["SE_Out"];
        }
    }

    virtual void output(PairedEndRead &per)  {
        if (per.non_const_read_one().getDiscard()) {
            ++c["PE_Discard"];
        } else {
            ++c["lins"];
            ++c["TotalFragmentsOutput"];
            ++c["PE_Out"];
        }
    }

    virtual void output(SingleEndRead &ser, unsigned int origLength)  {
        Read &one = ser.non_const_read_one();
        if (!one.getDiscard()) {
            if (one.getLength() < origLength) {
                ++c["sins"]; //adapters must be had (short insert)
                c["Adapter_BpTrim"] += (origLength - one.getLength());
            } else {
                ++c["mins"]; //must be a long insert
            }
            if ( one.getLength() + 1 > insertLength.size() ) {
                insertLength.resize(one.getLength() + 1);
            }
            ++insertLength[one.getLength()];
            ++c["SE_Out"];            
            ++c["TotalFragmentsOutput"];
        } else {
            ++c["PE_Discard"]; // originated as a PE read
        }
    }

    virtual void output(PairedEndRead &per, unsigned int overlapped) {
        //test lin or sin
        Read &one = per.non_const_read_one();
        Read &two = per.non_const_read_two();

        if (!one.getDiscard() && !two.getDiscard() ) {
            if (overlapped) {
                if (one.getLength() > overlapped || two.getLength() > overlapped ) {
                    ++c["sins"]; //adapters must be had (short insert)
                    c["Adapter_BpTrim"] += std::max((one.getLength() - one.getLengthTrue()),(two.getLength() - two.getLengthTrue()));

                } else {
                    ++c["mins"]; //must be a long insert
                }
                if ( overlapped + 1 > insertLength.size() ) {
                    insertLength.resize(overlapped + 1);
                }
                ++insertLength[overlapped];
                
            } else {
                ++c["lins"]; //lin
            }
            ++c["PE_Out"];
            ++c["TotalFragmentsOutput"];
        } else {
            ++c["PE_Discard"];
        }
    }

    virtual void write_out(const std::string &statsFile, bool appendStats, std::string program_name, std::string notes) {

        std::ifstream testEnd(statsFile);
        int end = testEnd.peek();
        testEnd.close();

        std::fstream outStats;

        if (appendStats && end != -1) {
            outStats.open(statsFile, std::ios::in | std::ios::out); //overwritte
            outStats.seekp(-6, std::ios::end );
            outStats << "  }, \"" << program_name << "_" << getpid()  << "\": {\n";
        } else {
            outStats.open(statsFile, std::ios::out | std::ios::trunc); //overwritt
            outStats << "{ \"" << program_name << "_" << getpid() <<  "\": {\n";
        }

        outStats << "    \"Notes\": \"" << notes << "\"";

        for (const auto &name : c) {
            outStats << ",\n"; //make sure json format is kept
            outStats << "    \"" << name.first << "\": " << name.second;
        }

        // embed instertLength (histogram) in sub json vector
        outStats << ",\n"; //make sure json format is kept
        outStats << "    \"histogram\": [";
        outStats << "[" << 1 << "," << insertLength[1] << "]";  // first, so as to keep the json comma convention

        for (size_t i = 2; i < insertLength.size(); ++i) {
            outStats << ", [" << i << "," << insertLength[i] << "]"; //make sure json format is kept
        }
        outStats << "]"; // finish off histogram
        outStats << "\n  }\n}\n";
        outStats.flush();
    }
};

class PhixCounters : public Counters {

public:
    PhixCounters() {
        Common();
        c["PE_hits"] = 0;
        c["SE_hits"] = 0;
        c["Inverse"] = 0;
    }
    void set_inverse() {
        c["Inverse"] = 1;
    }
    void inc_SE_hits() {
        ++c["SE_hits"];
    }
    void inc_PE_hits() {
        ++c["PE_hits"];
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
        c["PE_Discarded"] = 0;

    }

    void R1_stats(Read &one) {
        c["R1_Left_Trim"] += one.getLTrim();
        c["R1_Right_Trim"] += one.getRTrim();
    }

    void R2_stats(Read &two) {
        c["R2_Left_Trim"] += two.getLTrim();
        c["R2_Right_Trim"] += two.getRTrim();
    }

    void SE_stats(Read &se) {
        c["SE_Left_Trim"] += se.getLTrim();
        c["SE_Right_Trim"] += se.getRTrim();
    }

    void output(PairedEndRead &per, bool no_orphans = false) {
        Read &one = per.non_const_read_one();
        Read &two = per.non_const_read_two();
        if (!one.getDiscard() && !two.getDiscard()) {
            ++c["TotalFragmentsOutput"];
            ++c["PE_Out"];
            R1_stats(one);
            R2_stats(two);
        } else if (!one.getDiscard() && !no_orphans) { //if stranded RC
            ++c["TotalFragmentsOutput"];
            ++c["SE_Out"];
            ++c["R2_Discarded"];
            SE_stats(one);
        } else if (!two.getDiscard() && !no_orphans) { // Will never be RC
            ++c["TotalFragmentsOutput"];
            ++c["SE_Out"];
            ++c["R1_Discarded"];
            SE_stats(two);
        } else {
            ++c["PE_Discarded"];
        }
    }

    void output(SingleEndRead &ser) {
        Read &one = ser.non_const_read_one();
        if (!one.getDiscard()) {
            ++c["TotalFragmentsOutput"];
            ++c["SE_Out"];
            SE_stats(one);
        } else {
            ++c["SE_Discarded"];
        }
    }
    
};

class StatsCounters : public Counters {

public:
    Counter b;

    StatsCounters() {
        Common();
        b["A"] = 0;
        b["T"] = 0;
        b["C"] = 0;
        b["G"] = 0;
        b["N"] = 0;
        c["R1_BpLen"] = 0;
        c["R2_BpLen"] = 0;
        c["SE_BpLen"] = 0;
        c["R1_bQ30"] = 0;
        c["R2_bQ30"] = 0;
        c["SE_bQ30"] = 0;
    }
   
    void read_stats(Read &r) {
        for (auto bp : r.get_seq()) {
            switch (bp) {
                case 'A':
                    ++b["A"];
                    break;
                case 'T':
                    ++b["T"];
                    break;
                case 'C':
                    ++b["C"];
                    break;
                case 'G':
                    ++b["G"];
                    break;
                case 'N':
                    ++b["N"];
                    break;
                default:
                    throw std::runtime_error("Unknown bp in stats counter");
            }
        }
    }
 
    void q_stats(Read &r, unsigned long long int &val) {
    }
 
    void output(PairedEndRead &per) {
        Counters::output(per);
        Read &one = per.non_const_read_one();
        Read &two = per.non_const_read_two();
        read_stats(one);
        read_stats(two);
        int r1_q30bases=0;
        for (auto q : one.get_qual()) {
            r1_q30bases += (q - 33) >= 30;
        }
        int r2_q30bases=0;
        for (auto q : two.get_qual()) {
            r2_q30bases += (q - 33) >= 30;
        }
        c["R1_bQ30"] += r1_q30bases;
        c["R2_bQ30"] += r2_q30bases;
        c["R1_BpLen"] += one.getLength();
        c["R2_BpLen"] += two.getLength();
    }
    
    void output(SingleEndRead &ser) {
        Counters::output(ser);
        Read &one = ser.non_const_read_one();
        read_stats(one);
        int q30bases=0;
        for (auto q : one.get_qual()) {
            q30bases += (q - 33) >= 30;
        }
        c["SE_bQ30"] += q30bases;
        c["SE_BpLen"] += one.getLength();
    }

    virtual void write_out(const std::string &statsFile, bool appendStats, std::string program_name, std::string notes) {
        
        std::ifstream testEnd(statsFile);
        int end = testEnd.peek();
        testEnd.close();
        
        std::fstream outStats;
        
        if (appendStats && end != -1) {
            //outStats.open(statsFile, std::ofstream::out | std::ofstream::app); //overwritte
            outStats.open(statsFile, std::ios::in | std::ios::out); //overwritte
            outStats.seekp(-6, std::ios::end );
            outStats << "  }, \"" << program_name << "_" << getpid()  << "\": {\n";
        } else {
            //outStats.open(statsFile, std::ofstream::out); //overwritte
            outStats.open(statsFile, std::ios::out | std::ios::trunc); //overwrite
            outStats << "{ \"" << program_name << "_" << getpid() <<  "\": {\n";
        }
        outStats << "    \"Notes\": \"" << notes << "\"";

        for (const auto &name : c) {
            outStats << ",\n"; //make sure json format is kept
            outStats << "    \"" << name.first << "\": " << name.second;
        }
        // embed base composition as sub json vector
        outStats << ",\n"; //make sure json format is kept
        outStats << "    \"Base_Composition\": {\n";
        outStats << "        \"" << 'A' << "\": " << b["A"] << ",\n";
        outStats << "        \"" << 'T' << "\": " << b["T"] << ",\n";
        outStats << "        \"" << 'G' << "\": " << b["G"] << ",\n";
        outStats << "        \"" << 'C' << "\": " << b["C"] << ",\n";
        outStats << "        \"" << 'N' << "\": " << b["N"] << "\n";
        outStats << "    }"; // finish off histogram
        outStats << "\n  }\n}\n";
        outStats.flush();
    }
 
};

class SuperDeduperCounters : public Counters {
public:
    SuperDeduperCounters() {
        Common();
        c["Duplicate"] = 0;
        c["Ignored"] = 0;
    }

    void increment_replace() {
        ++c["Duplicate"];
    }

    void increment_ignored() {
        ++c["Ignored"];
    }
};


#endif
