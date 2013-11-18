#ifndef STR_MENU_AUX_H
#define STR_MENU_AUX_H

// STL.
#include <string>
#include <vector>

// External lib.
#include "boost/lexical_cast.hpp"

// Own files.
#include "structs.hpp"

using namespace std;
using namespace boost;

int GetCiEstimator(string &ciEstimatorStr, int ciEstimator);
int GetInterval(int &valMin, int &valMax, const string &intervalStr);
int GetFloat(double &val, const string &valStr);
int GetInt(int &val, const string &valStr);
int GetPath(string &out, const std::vector<string> &in, int &begin, int end);
int GetRegion(Region &region, const string &val);
int GetSkipLib(std::vector<string> &skipLibVec, string val);
int MissingArgument(const string &flag, int idxCurr, int argReq, int idxMax);

#endif // STR_MENU_AUX_H
