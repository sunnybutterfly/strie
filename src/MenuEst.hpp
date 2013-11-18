#ifndef STR_MENU_EST_H
#define STR_MENU_EST_H

// STL libraries.
#include <map>
#include <string>
#include <vector>

#include "boost_libraries.hpp"

using namespace std;

void MenuEstSrp(int argc, char *argv[]);
void SetParametersEstSrp(po::variables_map &vm);
void ShowMenuEstSrp(const po::options_description &desc);

#endif // STR_MENU_EST_H
