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
#include "version.h"
#include "read.h"
#include "typedefs.h"
#include <boost/filesystem/path.hpp>
#include <boost/program_options.hpp>

namespace bf = boost::filesystem;
namespace po = boost::program_options;

class Counters {
public:
    std::fstream outStats;
    std::string fStats = "/dev/null";
    bool force = false;
    bool aStats = false;

    std::string pName = "hts";
    std::string pNotes = "hts";
    po::variables_map vm;

    std::vector<Label> fragment;
    std::vector<Label> se;
    std::vector<Label> pe;
    std::vector<sLabel> pd;

    uint64_t TotalFragmentsInput = 0;
    uint64_t TotalFragmentsOutput = 0;

    uint64_t SE_In = 0;
    uint64_t SE_Out = 0;

    uint64_t PE_In = 0;
    uint64_t PE_Out = 0;

    Counters(const std::string &program_name_, po::variables_map vm_):
            pName(program_name_),
            vm(vm_) {

        if (vm.size()){
          fStats = vm["stats-file"].as<std::string>();
          force = vm["force"].as<bool>();
          aStats = vm["append-stats-file"].as<bool>();
          pNotes = vm["notes"].as<std::string>();
        }
        check_write();

        pd.push_back(std::forward_as_tuple("program", pName));
        pd.push_back(std::forward_as_tuple("version", VERSION));
        pd.push_back(std::forward_as_tuple("notes", pNotes));

        fragment.push_back(std::forward_as_tuple("totalFragmentsInput", TotalFragmentsInput));
        fragment.push_back(std::forward_as_tuple("totalFragmentsOutput", TotalFragmentsOutput));

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

        start_sublabel("Program_details");
        write_values(pd, 2);
        end_sublabel();

        start_sublabel("Fragment");
        write_values(fragment, 2);
        end_sublabel();

        start_sublabel("Single_end");
        write_values(se, 2);
        end_sublabel();

        start_sublabel("Paired_end");
        write_values(pe, 2);
        end_sublabel();

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
    }

    virtual void start_sublabel(const std::string &labelStr, const unsigned int indent = 1) {
        std::string pad(4 * indent, ' ');
        outStats << pad << "\"" << labelStr << "\": {\n";
    }

    virtual void end_sublabel(const unsigned int indent = 1) {
        std::string pad(4 * indent, ' ');
        outStats.seekp(-2, std::ios::end );
        outStats << "\n" << pad << "},\n"; // finish off histogram
    }

    template <class T>
    void write_values(T &tuples, const unsigned int indent = 1) {
        std::string pad(4 * indent, ' ');
        for (auto& tuple : tuples) {
            try
            {
                uint64_t value = boost::lexical_cast<uint64_t>(std::get<1>(tuple));
                outStats << pad << "\"" << std::get<0>(tuple) << "\": " << value << ",\n";
            }
            catch(boost::bad_lexical_cast &)
            {
                outStats << pad << "\"" << std::get<0>(tuple) << "\": \"" << std::get<1>(tuple) << "\",\n";
            }
        }
    }

    template <class T>
    void write_vector(const std::string &name, T &tuple, const unsigned int indent = 1) {
        std::string pad(4 * indent, ' ');
        if (tuple.size() == 0) return;
        size_t i;
        outStats << pad << "\"" << name << "\": [";
        for (i=0 ; i < tuple.size(); ++i) {
            try
            {
                uint64_t value = boost::lexical_cast<uint64_t>(std::get<1>(tuple[i]));
                outStats << " [" << std::get<0>(tuple[i]) << "," << value << "]"; //make sure json format is kept
            }
            catch(boost::bad_lexical_cast &)
            {
                outStats << " [\"" << std::get<0>(tuple[i]) << "\",\"" << std::get<1>(tuple[i]) << "\"]"; //make sure json format is kept
            }
            if (i < tuple.size()-1) outStats << ",";
        }
        outStats << " ],\n"; // finish off
    }

    virtual void write_matrix(const std::string &matrix_name, const Mat &data, const std::vector<std::string> &row_name, const std::vector<std::string> &col_name, const bool sparse = 0, const unsigned int indent = 1) {
        std::string pad(4 * indent, ' ');
        std::string pad2(4 * (indent + 1), ' ');
        if (data.size() == 0) return;

        if (data.size() < col_name.size()) throw std::runtime_error("In counters.h output: data matrix column size less than col_names size");
        if (data[0].size() < row_name.size()) throw std::runtime_error("In counters.h output: data matrix row size less than row_names size");

        outStats << pad << "\"" << matrix_name << "\": {\n";

        outStats << pad2 << "\"matrix_type\": \"" << ((sparse) ? "sparse" : "dense") << "\",\n";
        outStats << pad2 << "\"matrix_element_type\": \"int\",\n";
        outStats << pad2 << "\"shape\": [" << row_name.size() << ", " << col_name.size() << "],\n";
        outStats << pad2 << "\"row_names\": [";
        for (size_t i=0 ; i < row_name.size(); ++i) {
            outStats << " \"" << row_name[i] << "\""; //make sure json format is kept
            if (i < row_name.size()-1) outStats << ",";
        }
        outStats << " ],\n"; // finish off
        outStats << pad2 << "\"col_names\": [";
        for (size_t i=0 ; i < col_name.size(); ++i) {
            outStats << " \"" << col_name[i] << "\""; //make sure json format is kept
            if (i < col_name.size()-1) outStats << ",";
        }
        outStats << " ],\n"; // finish off
        outStats << pad2 << "\"" << "data" << "\": [";
        // sparse matrix
        if (sparse){
          bool first = 1;
          for (size_t i = 0 ; i < row_name.size(); ++i) {
              for (size_t j=0 ; j < col_name.size(); ++j ) {
                  if (data[j][i] != 0){
                      if (!first) outStats << ",";
                      outStats << " [" << i << "," << j << "," << data[j][i] << "]"; //make sure json format is kept
                      first = 0;
                  }
              }
          }
        } else {
          for (size_t i = 0 ; i < row_name.size(); ++i) {
              outStats << " [";
              for (size_t j=0 ; j < col_name.size(); ++j ) {
                  outStats << " " << data[j][i]; //make sure json format is kept
                  if (j < col_name.size()-1) outStats << ",";
              }
              outStats << " ]";
              if (i < row_name.size()-1) outStats << ",";
          }
        }
          outStats << " ]\n"; // finish off
        outStats << pad << "},\n";
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

    TrimmingCounters(const std::string &program_name, po::variables_map vm ) : Counters::Counters(program_name, vm) {
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
