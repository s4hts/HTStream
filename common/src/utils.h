#ifndef UTILS_H
#define UTILS_H

#include <unordered_map>
#include <string>
#include <iostream>
#include <fstream>

typedef std::unordered_map <std::string, size_t> Counter;

void setupCounter(Counter &c);
void write_stats(const std::string &statsFile, const bool &appendStats, const Counter &c, const std::string &program_name);

#endif
