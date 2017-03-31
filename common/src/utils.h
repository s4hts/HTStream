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
#include "version.h"
#include "typedefs.h"

namespace po = boost::program_options;


void setupCounter(Counter &c);
void write_stats(const std::string &statsFile, const bool &appendStats, const Counter &c, const std::string &program_name);

void outputWriters(std::shared_ptr<OutputWriter> &pe, std::shared_ptr<OutputWriter> &se, bool fastq_out, bool tab_out, bool interleaved_out, bool unmapped_out,  bool force, bool gzip_out, bool std_out, std::string &prefix);

void version_or_help(std::string program_name, po::options_description &desc, po::variables_map vm);

void setDefaultParams(po::options_description &desc, std::string program_name);
void setDefaultParamsCutting(po::options_description &desc);
void setDefaultParamsTrim(po::options_description &desc);

#endif
