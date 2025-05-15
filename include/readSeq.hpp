#ifndef READSEQ_HPP
#define READSEQ_HPP
#include <string>
#include <vector>
#include <fstream>
#include "Seq.hpp"


std::vector<Seq> get_sequence(const std::string& filename , int K);

#endif