#ifndef COUNTERS_H
#define COUNTERS_H


#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <tuple>
#include <unistd.h>
#include <unordered_map>
#include <vector>

#include "hts_exception.h"
#include "read.h"
#include "typedefs.h"
#include "version.h"

namespace bf = boost::filesystem;
namespace po = boost::program_options;

class Counters : public ReadVisitor {
public:
    std::fstream outStats;
    std::string fStats = "/dev/null";
    bool aStats = false;

    std::string pName = "hts";
    std::string pNotes = "hts";
    po::variables_map vm;

    std::vector<sLabel> pd;
    std::vector<Label> fragment;
    std::vector<Label> se;
    std::vector<Label> pe;
    std::vector<Label> r1;
    std::vector<Label> r2;

    uint64_t TotalFragmentsInput = 0;
    uint64_t TotalFragmentsOutput = 0;
    uint64_t TotalBasepairsInput = 0;
    uint64_t TotalBasepairsOutput = 0;

    uint64_t SE_In = 0;
    uint64_t SE_Out = 0;
    uint64_t SE_BpLen_In = 0;
    uint64_t SE_BpLen_Out = 0;

    uint64_t PE_In = 0;
    uint64_t PE_Out = 0;
    uint64_t R1_BpLen_In = 0;
    uint64_t R1_BpLen_Out = 0;
    uint64_t R2_BpLen_In = 0;
    uint64_t R2_BpLen_Out = 0;

    Counters(const std::string &program_name_, const po::variables_map& vm_):
            pName(program_name_),
            vm(vm_) {

        if (vm.size()){
          if (vm.count("append-stats-file")){
              fStats = vm["append-stats-file"].as<std::string>();
              aStats = true;
          } else {
              fStats = vm["stats-file"].as<std::string>();
          }
        }

        check_write();

        pd.push_back(std::forward_as_tuple("program", pName));
        pd.push_back(std::forward_as_tuple("version", VERSION));

        fragment.push_back(std::forward_as_tuple("in", TotalFragmentsInput));
        fragment.push_back(std::forward_as_tuple("out", TotalFragmentsOutput));
        fragment.push_back(std::forward_as_tuple("basepairs_in", TotalBasepairsInput));
        fragment.push_back(std::forward_as_tuple("basepairs_out", TotalBasepairsOutput));

        se.push_back(std::forward_as_tuple("in", SE_In));
        se.push_back(std::forward_as_tuple("out", SE_Out));
        se.push_back(std::forward_as_tuple("basepairs_in", SE_BpLen_In));
        se.push_back(std::forward_as_tuple("basepairs_out", SE_BpLen_Out));

        r1.push_back(std::forward_as_tuple("basepairs_in", R1_BpLen_In));
        r1.push_back(std::forward_as_tuple("basepairs_out", R1_BpLen_Out));
        r2.push_back(std::forward_as_tuple("basepairs_in", R2_BpLen_In));
        r2.push_back(std::forward_as_tuple("basepairs_out", R2_BpLen_Out));

        pe.push_back(std::forward_as_tuple("in", PE_In));
        pe.push_back(std::forward_as_tuple("out", PE_Out));
    }

    virtual ~Counters() {}

    virtual void input(PairedEndRead &per) {
        ++PE_In;
        Read &one = per.non_const_read_one();
        Read &two = per.non_const_read_two();
        R1_BpLen_In += one.getLength();
        R2_BpLen_In += two.getLength();
        TotalBasepairsInput += one.getLength();
        TotalBasepairsInput += two.getLength();
        ++TotalFragmentsInput;
    }

    virtual void input(SingleEndRead &ser) {
        ++SE_In;
        Read &one = ser.non_const_read_one();
        SE_BpLen_In += one.getLength();
        TotalBasepairsInput += one.getLength();
        ++TotalFragmentsInput;
    }

    virtual void input(ReadBase &read) {
        auto read_visit = make_read_visitor_func(
            [&](SingleEndRead *ser) {
                input(*ser);
            },
            [&](PairedEndRead *per) {
                input(*per);
            });
        read.accept(read_visit);
    }

    virtual void output(PairedEndRead &read) {
        Read &one = read.non_const_read_one();
        Read &two = read.non_const_read_two();
        R1_BpLen_Out += one.getLengthTrue();
        R2_BpLen_Out += two.getLengthTrue();
        TotalBasepairsOutput += one.getLengthTrue();
        TotalBasepairsOutput += two.getLengthTrue();
        ++PE_Out;
        ++TotalFragmentsOutput;
     }

    virtual void output(SingleEndRead &read) {
        Read &one = read.non_const_read_one();
        SE_BpLen_Out += one.getLengthTrue();
        TotalBasepairsOutput += one.getLengthTrue();
        ++SE_Out;
        ++TotalFragmentsOutput;
     }

    virtual void output(ReadBase &read) {
        read.accept(*this);
    }

    virtual void visit(PairedEndRead* per) {
        output(*per);
    }

    virtual void visit(SingleEndRead* ser) {
        output(*ser);
    }

    virtual void write_out() {

        initialize_json();

        start_sublabel("Program_details");
        write_values(pd, 2);
        start_sublabel("options", 2);
        write_options(3);
        end_sublabel(2);
        end_sublabel();

        start_sublabel("Fragment");
        write_values(fragment, 2);
        end_sublabel();

        start_sublabel("Single_end");
        write_values(se, 2);
        end_sublabel();

        start_sublabel("Paired_end");
        write_values(pe, 2);
        if (!r1.empty()){
            start_sublabel("Read1",2);
            write_values(r1, 3);
            end_sublabel(2);
        }
        if (!r2.empty()){
            start_sublabel("Read2",2);
            write_values(r2, 3);
            end_sublabel(2);
        }
        end_sublabel();

        finalize_json();
    }

    void initialize_json() {
        std::ifstream testEnd(fStats);
        int end = testEnd.peek();
        testEnd.close();

        if (aStats && end != -1) {
            outStats.open(fStats, std::ios::in | std::ios::out); //append
            outStats.seekp(-6, std::ios::end );
            outStats << "  }, {\n";
        } else {
            outStats.open(fStats, std::ios::out | std::ios::trunc); //overwrite
            outStats << "[ {\n";
        }
    }

    void start_sublabel(const std::string &labelStr, const unsigned int indent = 1) {
        std::string pad(4 * indent, ' ');
        outStats << pad << "\"" << labelStr << "\": {\n";
    }

    void end_sublabel(const unsigned int indent = 1) {
        std::string pad(4 * indent, ' ');
        outStats.seekp(-2, std::ios::end );
        outStats << "\n" << pad << "},\n"; // finish off histogram
    }

    void write_options(const unsigned int indent = 1){
        std::string pad(4 * indent, ' ');
        for (const auto& it : vm) {
            auto& value = it.second.value();
            //unfortunate hack
            if ((it.first == "stats-file") and vm.count("append-stats-file")) continue;
            outStats << pad << "\"" << it.first.c_str() << "\": ";
            if (auto v = boost::any_cast<std::string>(&value))
                outStats << "\"" << *v << "\"";
            else if (auto v = boost::any_cast<char>(&value)) // options case for chars (hts_ExtractUMI)
                outStats << "\"" << *v << "\"";
            else if (auto v = boost::any_cast<bool>(&value))    
                outStats << ((*v) ? "true" : "false");
            else if (auto v = boost::any_cast<size_t>(&value))
                outStats << *v;
            else if (auto v = boost::any_cast<double>(&value))
                outStats << *v ;
            else if (auto v = boost::any_cast<std::vector<std::string>>(&value)){
                outStats << "[ ";
                for (std::vector<std::string>::const_iterator x = v->begin(); x != v->end(); x++){
                    outStats << "\"" << *x << "\", ";
                }
                outStats.seekp(-2, std::ios::end );
                outStats << "]";
            } else
                throw HtsRuntimeException("In counters.h write_options: options type not seen in the switch, need to add new options type.");

            outStats << ",\n";
        }
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

    void write_matrix(const std::string &matrix_name, const Mat &data, const std::vector<std::string> &row_name, const std::vector<std::string> &col_name, const bool sparse = 0, const unsigned int indent = 1) {
        std::string pad(4 * indent, ' ');
        std::string pad2(4 * (indent + 1), ' ');
        if (data.size() == 0) return;

        if (data.size() < col_name.size()) throw HtsRuntimeException("In counters.h output: data matrix column size less than col_names size");
        if (data[0].size() < row_name.size()) throw HtsRuntimeException("In counters.h output: data matrix row size less than row_names size");

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

    void finalize_json() {
        outStats.seekp(-2, std::ios::end );
        outStats << "\n  }\n]\n";
        outStats.flush();
        outStats.close();
    }

private:
    void check_write() {
        std::fstream out;
        out.open(fStats, std::ios::out | std::ios::app);

        if(out.is_open())
        {
            out.close();
        }
        else
        {
            throw HtsIOException("Error: Cannot write to " + fStats + ": " +  std::strerror( errno ) + '\n');
        }
    }

};

class TrimmingCounters : public Counters {

public:
    uint64_t SE_Right_Trim = 0;
    uint64_t SE_Left_Trim = 0;

    uint64_t R1_Left_Trim = 0;
    uint64_t R1_Right_Trim = 0;
    uint64_t R2_Left_Trim = 0;
    uint64_t R2_Right_Trim = 0;

    TrimmingCounters(const std::string &program_name, po::variables_map vm ) : Counters::Counters(program_name, vm) {
        se.push_back(std::forward_as_tuple("rightTrim", SE_Right_Trim));
        se.push_back(std::forward_as_tuple("leftTrim", SE_Left_Trim));

        r1.push_back(std::forward_as_tuple("leftTrim", R1_Left_Trim));
        r1.push_back(std::forward_as_tuple("rightTrim", R1_Right_Trim));
        r2.push_back(std::forward_as_tuple("leftTrim", R2_Left_Trim));
        r2.push_back(std::forward_as_tuple("rightTrim", R2_Right_Trim));
    }

    virtual ~TrimmingCounters() {}

    virtual void R1_stats(Read &one) {
        R1_Left_Trim += one.getLTrim();
        R1_Right_Trim += one.getRTrim();
    }

    virtual void R2_stats(Read &two) {
        R2_Left_Trim += two.getLTrim();
        R2_Right_Trim += two.getRTrim();
    }

    virtual void SE_stats(Read &se) {
        SE_Left_Trim += se.getLTrim();
        SE_Right_Trim += se.getRTrim();
    }

    using Counters::output;
    virtual void output(PairedEndRead &per) {
        Counters::output(per);
        Read &one = per.non_const_read_one();
        Read &two = per.non_const_read_two();
        R1_stats(one);
        R2_stats(two);
    }

    virtual void output(SingleEndRead &ser) {
        Counters::output(ser);
        Read &one = ser.non_const_read_one();
        SE_stats(one);
    }

};

#endif
