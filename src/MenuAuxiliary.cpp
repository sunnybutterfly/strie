#include "MenuAuxiliary.hpp"

int GetCiEstimator(string &ciEstimatorStr, int ciEstimator)
{
  if (ciEstimator < 0 || ciEstimator > 2)
  {
    cout << "Confidence interval estimator must be denoted by an integer:"
        << endl;
    cout << "  0: none" << endl;
    cout << "  1: credible interval" << endl;
    cout << "  2: PDF percentile" << endl;
  }

  switch (ciEstimator)
  {
  case 0:
    ciEstimatorStr = "none";
    break;
  case 1:
    ciEstimatorStr = "credible interval";
    break;
  case 2:
    ciEstimatorStr = "PDF percentile";
    break;
  }

  return 0;
}

int GetFloat(double &val, const string &valStr)
{
  try
  {
    val = lexical_cast<double>(valStr.c_str());
  } catch (bad_lexical_cast &)
  {
    cout << "Not a float: \"" << valStr << "\"" << endl;
    return 1;
  }

  return 0;
}

int GetInterval(int &valMin, int &valMax, const string &intervalStr)
{
  size_t pos = intervalStr.find("-");
  string refBeg;
  string refEnd;
  size_t dist;
  int returnVal1;
  int returnVal2;

  if (pos == string::npos)
  {
    cout << "No interval entered" << endl;
    return 1;
  }

  dist = intervalStr.size() - min(pos + 1, intervalStr.size()) + 1;
  refBeg = intervalStr.substr(0, pos);
  refEnd = intervalStr.substr(pos + 1, dist);
  returnVal1 = GetInt(valMin, refBeg);
  returnVal2 = GetInt(valMax, refEnd);

  if (returnVal1 != 0 || returnVal2 != 0 || valMin >= valMax)
  {
    cout << "Improper interval: \"" << intervalStr << "\"" << endl;
    cout << "Correct interval: \"val1-val2\", where val1 < val2" << endl;
    return 2;
  }

  return 0;
}

int GetInt(int &val, const string &valStr)
{
  try
  {
    val = lexical_cast<int>(valStr.c_str());
  } catch (bad_lexical_cast &)
  {
    cout << "Not an integer: \"" << valStr << "\"" << endl;
    return 1;
  }

  return 0;
}

int GetPath(string &out, const std::vector<string> &in, int &begin, int end)
{
  string path;
  string val = in[begin + 1];
  char first = val[0];
  char last = val[val.size() - 1];

  if (first == '\"' && last != '\"')
  {
    path = val.substr(1);

    for (int i1 = begin + 1; i1 < end; i1++)
    {
      size_t pos = in[i1].find("\"");

      if (pos != string::npos)
      {
        path += in[i1];
        begin = i1;
        break;
      }

      path += in[i1];
    }
  }
  else if (first == '\"' && last == '\"')
  {
    path.assign(in[begin].begin() + 1, in[begin].end() - 1);
  }
  else
    path = val;

  out = path;

  return 0;
}

int GetRegion(Region &region, const string &val)
{
  size_t pos1 = val.find(":");
  size_t pos2 = val.find("-");
  bool error;

  if (pos1 > 0 && pos1 + 1 < pos2 && pos2 != val.size() - 1
      && pos1 != string::npos && pos2 != string::npos)
  {
    string refBegStr = val.substr(pos1 + 1, pos2 - pos1 - 1);
    string refEndStr = val.substr(pos2 + 1);

    region.refIdStr = val.substr(0, pos1);
    error = GetInt(region.refBeg, refBegStr);
    if (error)
    {
      cout << "Incorrect region: \"" << val << "\"" << endl;
      return 1;
    }
    error = GetInt(region.refEnd, refEndStr);
    if (error)
    {
      cout << "Incorrect region: \"" << val << "\"" << endl;
      return 1;
    }
  }
  else if (pos1 < val.size() - 1 && pos1 != string::npos
      && pos2 == string::npos)
  {
    string refBegStr = val.substr(pos1 + 1);

    region.refIdStr = val.substr(0, pos1);
    error = GetInt(region.refBeg, refBegStr);
    if (error)
    {
      cout << "Incorrect region: \"" << val << "\"" << endl;
      return 1;
    }
  }
  else if (val.size() > 0 && pos1 != string::npos && pos2 == string::npos)
  {
    region.refIdStr = val;
  }

  return 0;
}

int GetSkipLib(std::vector<string> &skipLibVec, string val)
{

  return 0;
}

int MissingArgument(const string &flag, int idxCurr, int argReq, int idxMax)
{
  if (idxCurr + argReq >= idxMax)
  {
    cout << "The flag \"" << flag << "\" requires " << argReq
        << " no. of arguments." << endl;
    return 1;
  }

  return 0;
}

