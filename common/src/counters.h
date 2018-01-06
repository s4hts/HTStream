#ifndef COUNTERS_H
#define COUNTERS_H


#include <unistd.h>
#include <map>
#include <unordered_map>
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <tuple>
#include "read.h"
#include "typedefs.h"
#include <boost/filesystem/path.hpp>

class Counters {
public:
    std::fstream outStats;
    std::string fStats;
    bool aStats;
    std::string pName;
    std::string pNotes;

    std::vector<Label> generic;
    std::vector<Label> se;
    std::vector<Label> pe;

    uint64_t TotalFragmentsInput = 0;
    uint64_t TotalFragmentsOutput = 0;

    uint64_t SE_In = 0;
    uint64_t SE_Out = 0;

    uint64_t PE_In = 0;
    uint64_t PE_Out = 0;

    Counters(const std::string &statsFile, bool appendStats, std::string program_name, std::string notes) {
        fStats.assign(statsFile);
        check_write();
        aStats = appendStats;
        pName = program_name;
        pNotes = notes;

        generic.push_back(std::forward_as_tuple("totalFragmentsInput", TotalFragmentsInput));
        generic.push_back(std::forward_as_tuple("totalFragmentsOutput", TotalFragmentsOutput));

        se.push_back(std::forward_as_tuple("SE_In", SE_In));
        se.push_back(std::forward_as_tuple("SE_Out", SE_Out));

        pe.push_back(std::forward_as_tuple("PE_In", PE_In));
        pe.push_back(std::forward_as_tuple("PE_Out", PE_Out));
    }

    virtual ~Counters() {}

    virtual void input(const ReadBase &read) {
        const PairedEndRead *per = dynamic_cast<const PairedEndRead *>(&read);
        if (per) {
            ++PE_In;
        } else {
            const SingleEndRead *ser = dynamic_cast<const SingleEndRead *>(&read);
            if (ser) {
                ++SE_In; 
            } else {
                throw std::runtime_error("In utils.h output: read type not valid");
            }
        }
        ++TotalFragmentsInput;
    }

    virtual void output(ReadBase &read) {
        PairedEndRead *per = dynamic_cast<PairedEndRead *>(&read);
        if (per) {
            ++PE_Out;
        } else {
            SingleEndRead *ser = dynamic_cast<SingleEndRead *>(&read);
            if (ser) {
                ++SE_Out;
            } else {
                throw std::runtime_error("In utils.h output: read type not valid");
            }
        }
        ++TotalFragmentsOutput;
    }

    virtual void write_out() {
        
        initialize_json();

        write_labels(generic);
        write_sublabels("Single_end", se);
        write_sublabels("Paired_end", pe);

        finalize_json();        
    }

    virtual void initialize_json() {
        std::ifstream testEnd(fStats);
        int end = testEnd.peek();
        testEnd.close();
        
        if (aStats && end != -1) {
            outStats.open(fStats, std::ios::in | std::ios::out); //append
            outStats.seekp(-6, std::ios::end );
            outStats << "  }, \"" << pName << "_" << getpid()  << "\": {\n";
        } else {
            outStats.open(fStats, std::ios::out | std::ios::trunc); //overwrite
            outStats << "{ \"" << pName << "_" << getpid() <<  "\": {\n";
        }

        outStats << "    \"Notes\": \"" << pNotes << "\",\n";
        // initialize should always be followed by finalize_json()
    }

    virtual void write_labels(std::vector<Label> &labels, const unsigned int indent = 1) {
        std::string pad(4 * indent, ' ');
        size_t i;
        for (i=0; i < labels.size(); ++i) {
            outStats << pad << "\"" << std::get<0>(labels[i]) << "\": " << std::get<1>(labels[i]) << ",\n";
        }
    }

    virtual void write_sublabels(std::string labelStr, std::vector<Label> &labels, const unsigned int indent = 1) {
        std::string pad(4 * indent, ' ');
        outStats << pad << "\"" << labelStr << "\": {\n";
        write_labels(labels, indent+1);
        outStats.seekp(-2, std::ios::end );
        outStats << "\n" << pad << "},\n"; // finish off histogram
    }

    virtual void write_vector(std::string vector_name, std::vector<std::tuple<uint_fast64_t, uint_fast64_t>> vectortuple, const unsigned int indent = 1) {
        size_t i;
        std::string pad(4 * indent, ' ');
        outStats << pad << "\""<< vector_name << "\": [";
        for (i=0 ; i < vectortuple.size()-1; ++i) {
            outStats << " [" << std::get<0>(vectortuple[i]) << "," << std::get<1>(vectortuple[i]) << "],"; //make sure json format is kept
        }
        outStats << " [" << std::get<0>(vectortuple[i]) << "," << std::get<1>(vectortuple[i]) << "]";  // first, so as to keep the json comma convention
        outStats << " ],\n"; // finish off histogram
    }

    virtual void finalize_json() {
        outStats.seekp(-2, std::ios::end );
        outStats << "\n  }\n}\n";
        outStats.flush();
        outStats.close();
    }

private:
    virtual void check_write() {
        FILE* f = NULL;
        
        f = fopen(fStats.c_str(), "w");
        if (!f) {
            throw std::runtime_error("cannot write to " + fStats + ": " +  std::strerror( errno ));
        }
        fclose (f);
    }

};

class OverlappingCounters : public Counters {

public:
    std::vector<uint_fast64_t> insertLength;

    uint64_t sins = 0;
    uint64_t mins = 0;
    uint64_t lins = 0;
    uint64_t Adapter_BpTrim = 0;

    uint64_t SE_Discard = 0;

    uint64_t PE_Discard = 0;

    OverlappingCounters(const std::string &statsFile, bool appendStats, std::string program_name, std::string notes) : Counters::Counters(statsFile, appendStats, program_name, notes) {

        insertLength.resize(1);

        generic.push_back(std::forward_as_tuple("sins", sins));
        generic.push_back(std::forward_as_tuple("mins", mins));
        generic.push_back(std::forward_as_tuple("lins", lins));
        generic.push_back(std::forward_as_tuple("adapterBpTrim", Adapter_BpTrim));

        se.push_back(std::forward_as_tuple("SE_Discard", SE_Discard));

        pe.push_back(std::forward_as_tuple("PE_Discard", PE_Discard));
    }

    using Counters::output;
    virtual void output(SingleEndRead &ser)  {
        if (ser.non_const_read_one().getDiscard()) {
            ++SE_Discard;
        } else {
            ++TotalFragmentsOutput;
            ++SE_Out;
        }
    }

    virtual void output(PairedEndRead &per)  {
        if (per.non_const_read_one().getDiscard()) {
            ++PE_Discard;
        } else {
            ++lins;
            ++TotalFragmentsOutput;
            ++PE_Out;
        }
    }

    virtual void output(SingleEndRead &ser, uint_fast64_t origLength)  {
        Read &one = ser.non_const_read_one();
        if (!one.getDiscard()) {
            if (one.getLength() < origLength) {
                ++sins; //adapters must be had (short insert)
                Adapter_BpTrim += (origLength - one.getLength());
            } else {
                ++mins; //must be a long insert
            }
            if ( one.getLength() + 1 > insertLength.size() ) {
                insertLength.resize(one.getLength() + 1);
            }
            ++insertLength[one.getLength()];
            ++SE_Out;            
            ++TotalFragmentsOutput;
        } else {
            ++PE_Discard; // originated as a PE read
        }
    }

    virtual void output(PairedEndRead &per, uint_fast64_t overlapped) {
        //test lin or sin
        Read &one = per.non_const_read_one();
        Read &two = per.non_const_read_two();

        if (!one.getDiscard() && !two.getDiscard() ) {
            if (overlapped) {
                if (one.getLength() > overlapped || two.getLength() > overlapped ) {
                    ++sins; //adapters must be had (short insert)
                    Adapter_BpTrim += std::max((one.getLength() - one.getLengthTrue()),(two.getLength() - two.getLengthTrue()));

                } else {
                    ++mins; //must be a long insert
                }
                if ( overlapped + 1 > insertLength.size() ) {
                    insertLength.resize(overlapped + 1);
                }
                ++insertLength[overlapped];
                
            } else {
                ++lins; //lin
            }
            ++PE_Out;
            ++TotalFragmentsOutput;
        } else {
            ++PE_Discard;
        }
    }

    virtual void write_out() {

        std::vector<std::tuple<uint_fast64_t, uint_fast64_t>> iLength;
        for (size_t i = 1; i < insertLength.size(); ++i) {
            if (insertLength[i] > 0) {
                iLength.push_back(std::forward_as_tuple(i, insertLength[i]));
            }
        }

        initialize_json();

        write_labels(generic);
        write_vector("readlength_histogram",iLength);
        write_sublabels("Single_end", se);
        write_sublabels("Paired_end", pe);

        finalize_json();        
    }
};

class PhixCounters : public Counters {

public:
    uint64_t Inverse = 0;

    uint64_t SE_hits = 0;

    uint64_t PE_hits = 0;

    PhixCounters(const std::string &statsFile, bool appendStats, std::string program_name, std::string notes) : Counters::Counters(statsFile, appendStats, program_name, notes) {
        generic.push_back(std::forward_as_tuple("inverse", Inverse));

        se.push_back(std::forward_as_tuple("SE_hits", SE_hits));

        pe.push_back(std::forward_as_tuple("PE_hits", PE_hits));
    }

    void set_inverse() {
        Inverse = 1;
    }
    void inc_SE_hits() {
        ++SE_hits;
    }
    void inc_PE_hits() {
        ++PE_hits;
    }
};

class TrimmingCounters : public Counters {

public: 
    uint64_t SE_Right_Trim = 0;
    uint64_t SE_Left_Trim = 0;
    uint64_t SE_Discarded = 0;

    uint64_t R1_Left_Trim = 0;
    uint64_t R1_Right_Trim = 0;
    uint64_t R2_Left_Trim = 0;
    uint64_t R2_Right_Trim = 0;    
    uint64_t R1_Discarded = 0;
    uint64_t R2_Discarded = 0;
    uint64_t PE_Discarded = 0;

    TrimmingCounters(const std::string &statsFile, bool appendStats, std::string program_name, std::string notes) : Counters::Counters(statsFile, appendStats, program_name, notes) {
        se.push_back(std::forward_as_tuple("SE_rightTrimm", SE_Right_Trim));
        se.push_back(std::forward_as_tuple("SE_leftTrim", SE_Left_Trim));
        se.push_back(std::forward_as_tuple("SE_discarded", SE_Discarded));

        pe.push_back(std::forward_as_tuple("R1_leftTrim", R1_Left_Trim));
        pe.push_back(std::forward_as_tuple("R1_rightTrim", R1_Right_Trim));
        pe.push_back(std::forward_as_tuple("R2_leftTrim", R2_Left_Trim));
        pe.push_back(std::forward_as_tuple("R2_rightTrim", R2_Right_Trim));
        pe.push_back(std::forward_as_tuple("R1_discarded", R1_Discarded));
        pe.push_back(std::forward_as_tuple("R2_discarded", R2_Discarded));
        pe.push_back(std::forward_as_tuple("PE_discarded", PE_Discarded));
    }

    void R1_stats(Read &one) {
        R1_Left_Trim += one.getLTrim();
        R1_Right_Trim += one.getRTrim();
    }

    void R2_stats(Read &two) {
        R2_Left_Trim += two.getLTrim();
        R2_Right_Trim += two.getRTrim();
    }

    void SE_stats(Read &se) {
        SE_Left_Trim += se.getLTrim();
        SE_Right_Trim += se.getRTrim();
    }

    using Counters::output;
    void output(PairedEndRead &per, bool no_orphans = false) {
        Read &one = per.non_const_read_one();
        Read &two = per.non_const_read_two();
        if (!one.getDiscard() && !two.getDiscard()) {
            ++TotalFragmentsOutput;
            ++PE_Out;
            R1_stats(one);
            R2_stats(two);
        } else if (!one.getDiscard() && !no_orphans) { //if stranded RC
            ++TotalFragmentsOutput;
            ++SE_Out;
            ++R2_Discarded;
            SE_stats(one);
        } else if (!two.getDiscard() && !no_orphans) { // Will never be RC
            ++TotalFragmentsOutput;
            ++SE_Out;
            ++R1_Discarded;
            SE_stats(two);
        } else {
            ++PE_Discarded;
        }
    }

    void output(SingleEndRead &ser) {
        Read &one = ser.non_const_read_one();
        if (!one.getDiscard()) {
            ++TotalFragmentsOutput;
            ++SE_Out;
            SE_stats(one);
        } else {
            ++SE_Discarded;
        }
    }
    
};

#endif
