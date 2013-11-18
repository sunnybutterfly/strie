#ifndef MENUFREQ_HPP_
#define MENUFREQ_HPP_

// STL libraries.
#include <string>
#include <vector>

#include "boost_libraries.hpp"

using namespace std;

void MenuFreqSrp(int argc, char *argv[]);
void SetParametersFreqSrp(po::variables_map &vm);
void ShowMenuFreqSrp(const po::options_description &desc);

#endif /* MENUFREQ_HPP_ */
