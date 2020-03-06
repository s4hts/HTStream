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

namespace bf = boost::filesystem;

class Counters {
public:
    std::fstream outStats;
    std::string fStats;
    bool force;
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

    Counters(const std::string &statsFile, bool force_, bool appendStats, const std::string &program_name, const std::string &notes):
            fStats(statsFile),
            force(force_),
            aStats(appendStats),
            pName(program_name),
            pNotes(notes) {

        check_write();

        generic.push_back(std::forward_as_tuple("totalFragmentsInput", TotalFragmentsInput));
        generic.push_back(std::forward_as_tuple("totalFragmentsOutput", TotalFragmentsOutput));

        se.push_back(std::forward_as_tuple("SE_in", SE_In));
        se.push_back(std::forward_as_tuple("SE_out", SE_Out));

        pe.push_back(std::forward_as_tuple("PE_in", PE_In));
        pe.push_back(std::forward_as_tuple("PE_out", PE_Out));
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

    virtual void output(PairedEndRead &read, bool no_orphans = false) {
        (void)read;  //ignore unused variable warning
        (void)no_orphans;  //ignore unused variable warning
        ++PE_Out;
        ++TotalFragmentsOutput;
     }

    virtual void output(SingleEndRead &read) {
        (void)read;  //ignore unused variable warning
        ++SE_Out;
        ++TotalFragmentsOutput;
     }
    
    virtual void output(ReadBase &read) {
        PairedEndRead *per = dynamic_cast<PairedEndRead *>(&read);
        if (per) {
            return output(*per);
        } else {
            SingleEndRead *ser = dynamic_cast<SingleEndRead *>(&read);
            if (ser) {
                return output(*ser);
            } else {
                throw std::runtime_error("In utils.h output: read type not valid");
            }
        }
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

    virtual void write_labels(const std::vector<Label> &labels, const unsigned int indent = 1) {
        std::string pad(4 * indent, ' ');
        for (auto& label : labels) {
            outStats << pad << "\"" << std::get<0>(label) << "\": " << std::get<1>(label) << ",\n";
        }
    }

    virtual void write_sublabels(const std::string &labelStr, const std::vector<Label> &labels, const unsigned int indent = 1) {
        std::string pad(4 * indent, ' ');
        outStats << pad << "\"" << labelStr << "\": {\n";
        write_labels(labels, indent+1);
        outStats.seekp(-2, std::ios::end );
        outStats << "\n" << pad << "},\n"; // finish off histogram
    }

    template <class T>
    void write_vector(const std::string &name, T &tuple, const unsigned int indent = 1) {
        std::string pad(4 * indent, ' ');
        if (tuple.size() == 0) return;
        size_t i;
        outStats << pad << "\"" << name << "\": [";
        for (i=0 ; i < tuple.size()-1; ++i) {
            outStats << " [" << std::get<0>(tuple[i]) << "," << std::get<1>(tuple[i]) << "],"; //make sure json format is kept
        }
        outStats << " [" << std::get<0>(tuple[i]) << "," << std::get<1>(tuple[i]) << "]";  // first, so as to keep the json comma convention
        outStats << " ],\n"; // finish off histogram
    }

    virtual void write_vector_slabel(const std::string &vector_name, const std::vector<sLabel> &labeltuple, const unsigned int indent = 1) {
        std::string pad(4 * indent, ' ');
        if (labeltuple.size() == 0) return;
        size_t i;
        outStats << pad << "\"" << vector_name << "\": [";
        for (i=0 ; i < labeltuple.size()-1; ++i) {
            outStats << " [\"" << std::get<0>(labeltuple[i]) << "\",\"" << std::get<1>(labeltuple[i]) << "\"],"; //make sure json format is kept
        }
        outStats << " [\"" << std::get<0>(labeltuple[i]) << "\",\"" << std::get<1>(labeltuple[i]) << "\"]";  // first, so as to keep the json comma convention
        outStats << " ],\n"; // finish off
    }

    virtual void finalize_json() {
        outStats.seekp(-2, std::ios::end );
        outStats << "\n  }\n}\n";
        outStats.flush();
        outStats.close();
    }

private:
    virtual void check_write() {
        bf::path p(fStats);
        outStats.open(fStats, std::ios::out | std::ios::app);

        if(outStats.is_open())
        {
            outStats.close();
        }
        else
        {
            throw std::runtime_error("Error: Cannot write to " + fStats + ": " +  std::strerror( errno ) + '\n');
        }
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

    TrimmingCounters(const std::string &statsFile, bool force, bool appendStats, const std::string &program_name, const std::string &notes) : Counters::Counters(statsFile, force, appendStats, program_name, notes) {
        se.push_back(std::forward_as_tuple("SE_rightTrim", SE_Right_Trim));
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
    virtual ~TrimmingCounters() {}

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
    virtual void output(PairedEndRead &per, bool no_orphans = false) {
        Read &one = per.non_const_read_one();
        Read &two = per.non_const_read_two();
        if (!one.getDiscard() && !two.getDiscard()) {
            ++TotalFragmentsOutput;
            ++PE_Out;
            R1_stats(one);
            R2_stats(two);
        } else if (!one.getDiscard() && !no_orphans) {
            ++TotalFragmentsOutput;
            ++SE_Out;
            ++R2_Discarded;
            SE_stats(one);
        } else if (!two.getDiscard() && !no_orphans) {
            ++TotalFragmentsOutput;
            ++SE_Out;
            ++R1_Discarded;
            SE_stats(two);
        } else {
            ++PE_Discarded;
        }
    }

    virtual void output(SingleEndRead &ser) {
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
