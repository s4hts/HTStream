#ifndef COUNTERS_H
#define COUNTERS_H

#include <map>
#include <unordered_map>
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <tuple>
#include "read.h"
#include "typedefs.h"
#include <unistd.h>

class Counters {
public:
    std::fstream outStats;

    uint64_t TotalFragmentsInput = 0;
    uint64_t PE_In = 0;
    uint64_t SE_In = 0;
    uint64_t TotalFragmentsOutput = 0;
    uint64_t PE_Out = 0;
    uint64_t SE_Out = 0;

    std::vector<Label> labels;
    Counters() {
        labels.push_back(std::forward_as_tuple("TotalFragmentsInput", TotalFragmentsInput));
        labels.push_back(std::forward_as_tuple("PE_In", PE_In));
        labels.push_back(std::forward_as_tuple("SE_In", SE_In));
        labels.push_back(std::forward_as_tuple("TotalFragmentsOutput", TotalFragmentsOutput));
        labels.push_back(std::forward_as_tuple("PE_Out", PE_Out));
        labels.push_back(std::forward_as_tuple("SE_Out", SE_Out));
    }

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

    virtual void output(PairedEndRead &read, bool no_orhpans = false) {
        ++TotalFragmentsOutput;
        ++PE_Out;
    }


    virtual void output(SingleEndRead &read) {
        ++TotalFragmentsOutput;
        ++SE_Out;
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
        
        initialize_json(statsFile, appendStats, program_name, notes);

        write_labels();

        finalize_json();        
    }

    virtual const void initialize_json(const std::string &statsFile, bool appendStats, std::string program_name, std::string notes) {
        std::ifstream testEnd(statsFile);
        int end = testEnd.peek();
        testEnd.close();
        
        if (appendStats && end != -1) {
            outStats.open(statsFile, std::ios::in | std::ios::out); //overwritte
            outStats.seekp(-6, std::ios::end );
            outStats << "  }, \"" << program_name << "_" << getpid()  << "\": {\n";
        } else {
            outStats.open(statsFile, std::ios::out | std::ios::trunc); //overwrite
            outStats << "{ \"" << program_name << "_" << getpid() <<  "\": {\n";
        }

        outStats << "    \"Notes\": \"" << notes << "\"";
    }

    virtual const void write_labels() {
        for (const auto &label : labels) {
            outStats << ",\n"; //make sure json format is kept
            outStats << "    \"" << std::get<0>(label) << "\": " << std::get<1>(label);
        }
    }

    virtual const void write_vector(std::string vector_name, std::vector<std::tuple<uint_fast64_t, uint_fast64_t>> vectortuple) {
        // embed duplicate saturation in sub json vector
        outStats << ",\n"; //make sure json format is kept
        outStats << "    \""<< vector_name << "\": [";
        outStats << "[" << std::get<0>(vectortuple[1]) << "," << std::get<1>(vectortuple[1]) << "]";  // first, so as to keep the json comma convention

        for (size_t i = 2; i < vectortuple.size(); ++i) {
            outStats << ", [" << std::get<0>(vectortuple[i]) << "," << std::get<1>(vectortuple[i]) << "]"; //make sure json format is kept
        }
        outStats << "]"; // finish off histogram
    }

    virtual const void finalize_json() {
        outStats << "\n  }\n}\n";
        outStats.flush();
    }

};


class OverlappingCounters : public Counters {

public:
    std::vector<unsigned long long int> insertLength;
    uint64_t sins = 0;
    uint64_t mins = 0;
    uint64_t lins = 0;
    uint64_t SE_Discard = 0;
    uint64_t PE_Discard = 0;
    uint64_t Adapter_BpTrim = 0;

    OverlappingCounters() {
        labels.push_back(std::forward_as_tuple("sins", sins));
        labels.push_back(std::forward_as_tuple("mins", mins));
        labels.push_back(std::forward_as_tuple("lins", lins));
        labels.push_back(std::forward_as_tuple("SE_Discard", SE_Discard));
        labels.push_back(std::forward_as_tuple("PE_Discard", PE_Discard));
        labels.push_back(std::forward_as_tuple("Adapter_BpTrim", Adapter_BpTrim));

        insertLength.resize(1);
    }

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

    virtual void output(SingleEndRead &ser, unsigned int origLength)  {
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

    virtual void output(PairedEndRead &per, unsigned int overlapped) {
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

    virtual void write_out(const std::string &statsFile, bool appendStats, std::string program_name, std::string notes) {

        std::vector<std::tuple<uint_fast64_t, uint_fast64_t>> iLength;

        initialize_json(statsFile, appendStats, program_name, notes);

        write_labels();

        for (size_t i = 1; i < insertLength.size(); ++i) {
            if (insertLength[i] > 0) {
                iLength.push_back(std::forward_as_tuple(i, insertLength[i]));
            }
        }
        write_vector("readlength_histogram",iLength);

        finalize_json();        
    }
};

class PhixCounters : public Counters {

public:
    uint64_t PE_hits = 0;
    uint64_t SE_hits = 0;
    uint64_t Inverse = 0;

    PhixCounters() {
        labels.push_back(std::forward_as_tuple("PE_hits", PE_hits));
        labels.push_back(std::forward_as_tuple("SE_hits", SE_hits));
        labels.push_back(std::forward_as_tuple("Inverse", Inverse));
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
    uint64_t R1_Left_Trim = 0;
    uint64_t R1_Right_Trim = 0;
    uint64_t R2_Left_Trim = 0;
    uint64_t R2_Right_Trim = 0;
    uint64_t SE_Right_Trim = 0;
    uint64_t SE_Left_Trim = 0;
    
    uint64_t R1_Discarded = 0;
    uint64_t R2_Discarded = 0;
    uint64_t SE_Discarded = 0;
    uint64_t PE_Discarded = 0;

    TrimmingCounters() {
        labels.push_back(std::forward_as_tuple("R1_Left_Trim", R1_Left_Trim));
        labels.push_back(std::forward_as_tuple("R1_Right_Trim", R1_Right_Trim));
        labels.push_back(std::forward_as_tuple("R2_Left_Trim", R2_Left_Trim));
        labels.push_back(std::forward_as_tuple("R2_Right_Trim", R2_Right_Trim));
        labels.push_back(std::forward_as_tuple("SE_Right_Trimm", SE_Right_Trim));
        labels.push_back(std::forward_as_tuple("SE_Left_Trimm", SE_Left_Trim));
        labels.push_back(std::forward_as_tuple("R1_Discarded", R1_Discarded));
        labels.push_back(std::forward_as_tuple("R2_Discarded", R2_Discarded));
        labels.push_back(std::forward_as_tuple("SE_Discarded", SE_Discarded));
        labels.push_back(std::forward_as_tuple("PE_Discarded", PE_Discarded));

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

class StatsCounters : public Counters {

public:
    uint64_t A = 0;
    uint64_t C = 0;
    uint64_t T = 0;
    uint64_t G = 0;
    uint64_t N = 0;
    uint64_t R1_BpLen = 0;
    uint64_t R2_BpLen = 0;
    uint64_t SE_BpLen = 0;
    uint64_t R1_bQ30 = 0;
    uint64_t R2_bQ30 = 0;
    uint64_t SE_bQ30 = 0;

    StatsCounters () {
        labels.push_back(std::forward_as_tuple("R1_BpLen", R1_BpLen));
        labels.push_back(std::forward_as_tuple("R2_BpLen", R2_BpLen));
        labels.push_back(std::forward_as_tuple("SE_BpLen", SE_BpLen));
        labels.push_back(std::forward_as_tuple("R1_bQ30", R1_bQ30));
        labels.push_back(std::forward_as_tuple("R2_bQ30", R2_bQ30));
        labels.push_back(std::forward_as_tuple("SE_bQ30", SE_bQ30));
    };


    void read_stats(Read &r) {
        for (auto bp : r.get_seq()) {
            switch (bp) {
                case 'A':
                    ++A;
                    break;
                case 'T':
                    ++T;
                    break;
                case 'C':
                    ++C;
                    break;
                case 'G':
                    ++G;
                    break;
                case 'N':
                    ++N;
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
        R1_bQ30 += r1_q30bases;
        R2_bQ30 += r2_q30bases;
        R1_BpLen += one.getLength();
        R2_BpLen += two.getLength();
    }
    
    void output(SingleEndRead &ser) {
        Counters::output(ser);
        Read &one = ser.non_const_read_one();
        read_stats(one);
        int q30bases=0;
        for (auto q : one.get_qual()) {
            q30bases += (q - 33) >= 30;
        }
        SE_bQ30 += q30bases;
        SE_BpLen += one.getLength();
    }

    virtual void write_out(const std::string &statsFile, bool appendStats, std::string program_name, std::string notes) {
        
        initialize_json(statsFile, appendStats, program_name, notes);

        write_labels();

        // embed base composition as sub json vector
        outStats << ",\n"; //make sure json format is kept
        outStats << "    \"Base_Composition\": {\n";
        outStats << "        \"" << 'A' << "\": " << A << ",\n";
        outStats << "        \"" << 'T' << "\": " << T << ",\n";
        outStats << "        \"" << 'G' << "\": " << G << ",\n";
        outStats << "        \"" << 'C' << "\": " << C << ",\n";
        outStats << "        \"" << 'N' << "\": " << N << "\n";
        outStats << "    }"; // finish off histogram

        finalize_json();        
    }
 
};

class SuperDeduperCounters : public Counters {

public:
    std::vector<std::tuple<uint_fast64_t, uint_fast64_t>> duplicateProportion;
    uint64_t Duplicate = 0;
    uint64_t Ignored = 0;

    SuperDeduperCounters() {
        labels.push_back(std::forward_as_tuple("Duplicate", Duplicate));
        labels.push_back(std::forward_as_tuple("Ignored", Ignored));
    }

    virtual void input(const ReadBase &read, size_t dup_freq) {

        if (dup_freq > 0 & (TotalFragmentsInput - Ignored) % dup_freq == 0){
            duplicateProportion.push_back(std::forward_as_tuple((TotalFragmentsInput - Ignored), Duplicate));
        }
        Counters::input(read);
    }

    void increment_replace() {
        ++Duplicate;
    }

    void increment_ignored() {
        ++Ignored;
    }

    virtual void write_out(const std::string &statsFile, bool appendStats, std::string program_name, std::string notes) {

        initialize_json(statsFile, appendStats, program_name, notes);

        write_labels();

        // record final input/dup
        duplicateProportion.push_back(std::forward_as_tuple((TotalFragmentsInput - Ignored), Duplicate));

        write_vector("duplicate_saturation",duplicateProportion);

        finalize_json();        
    }
};


#endif
