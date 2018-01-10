#ifndef UTILS_H
#define UTILS_H

#include "ioHandler.h"
#include <unordered_map>
#include <string>
#include <iostream>
#include <fstream>
#include <memory>
#include <vector>
#include <boost/program_options.hpp>
#include <boost/bind.hpp>
#include "version.h"

namespace po = boost::program_options;

void outputWriters(std::shared_ptr<OutputWriter> &pe, std::shared_ptr<OutputWriter> &se, bool fastq_out, bool tab_out, bool interleaved_out, bool unmapped_out,  bool force, bool gzip_out, bool std_out, std::string &prefix);

void version_or_help(std::string program_name, std::string app_description, po::options_description &desc, po::variables_map vm, bool error = false);

po::options_description setStandardOptions();
po::options_description setInputOptions();
po::options_description setOutputOptions(std::string program_name);

void setDefaultParamsCutting(po::options_description &desc);
void setDefaultParamsTrim(po::options_description &desc);
void setDefaultParamsOverlapping(po::options_description &desc);

template<typename T>
void check_range(const std::string& name, const T& value, const T& min, const T& max)
{
   if (value < min || value > max)
   {
        throw po::validation_error(po::validation_error::invalid_option_value, name);
      // throw exception
   }
}


char rc(const char &bp);

bool threshold_mismatches(std::string::const_iterator r1, std::string::const_iterator r2, size_t length, size_t max); 

typedef std::unordered_multimap<std::string, std::size_t> seqLookup;

seqLookup readOneMap(std::string seq1, const size_t kmer, const size_t kmerOffset);

#endif
