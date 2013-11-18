#ifndef STRINPUT_HPP_
#define STRINPUT_HPP_

// STL libraries
#include <cmath>        // log(), exp()
#include <cstdlib>      // abs(), rand
#include <fstream>      // open(), close()
#include <iomanip>      // Insertion and extraction operators.
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

// External libraries
#include <boost/multi_array.hpp>
#include <boost/lexical_cast.hpp>
#include <sys/time.h>   // gettimeofday

// Own files
#include "type.hpp"

// Namespace
using namespace std;
using namespace boost;

// Global variables.

// Class
class StrInput
{
public:
  StrInput();
  ~StrInput();
  StrInput(const StrInput & rhs);

// Protected variables.
protected:

  int ReadInsertSize(
      std::vector<StrLibrary> & libVec,
      int & insertSizeMin,
      int & insertSizeMax,
      int & mapqMin,
      const filesystem::path &f_path);
  int ReadLocus(
      deque<Region> & strLocDeque,
      const filesystem::path &f_path);
  int ReadLocusRefDon(
      deque<pair<Region, Region> > & strLocDeque,
      const filesystem::path &f_path);
  int ReadPriorPdf(
      list<T> & priorList,
      const filesystem::path &f_path,
      int nLineSkip,
      int indelSizeMax,
      int indelStepSize);
};

StrInput::StrInput()
{

}

StrInput::~StrInput()
{

}

StrInput::StrInput(const StrInput & rhs)
{

}

int StrInput::ReadInsertSize(
    std::vector<StrLibrary> & libVec,
    int & insertSizeMin,
    int & insertSizeMax,
    int & mapqMin,
    const filesystem::path &f_path)
{
  string fLoc = f_path.string();
  ifstream fInsertSize(fLoc.c_str());
  timeval t1, t2;
  const int bufferSize = 16384;
  char *line = new char[bufferSize];
  //int len;

  int nLib = 0;
  int count = 0;
  vector < string > libTagVec;

  if (!fInsertSize.is_open())
  {
    cout << "Cannot open the file containing insert size frequencies: ";
    cout << "\"" << fLoc << "\"." << endl;
    return 1;
  }

  // Get length of file.
  fInsertSize.seekg(0, ios::end);
  //len = fInsertSize.tellg();
  fInsertSize.seekg(0, ios::beg);

  // Start timer.
  gettimeofday(&t1, 0);

  while (!fInsertSize.eof())
  {
    string str;
    string substr;
    istringstream iss;

    for (int i0 = 0; i0 < bufferSize; i0++)
      line[i0] = '\0';

    fInsertSize.getline(line, bufferSize);
    str.assign(line);

    if (str.size() == 0)
      continue;

    if (str[0] == '>')
    {
      int pos = str.find(':');
      iss.str(str.substr(pos + 2));

      if (str.find("Libraries:") != string::npos)
      {
        string libTag;

        while (true)
        {
          iss >> libTag;

          if (libTag.size() > 0)
            libTagVec.push_back(libTag);
          else
            break;

          libTag.clear();
          nLib++;
        }
      }
      else if (str.find("MAPQ minimum:") != string::npos)
      {
        string mapqMinStr;

        iss >> mapqMinStr;

        try
        {
          mapqMin = lexical_cast<int>(mapqMinStr);
        } catch (bad_lexical_cast &)
        {
          cout << "Error at line " << count + 1 << " when reading file \""
              << fLoc << "\"." << endl;
          return 1;
        }
      }
      else if (str.find("Insert size interval:") != string::npos)
      {
        substr = str.substr(pos + 2);
        pos = substr.find('[');
        substr = substr.substr(pos + 1);
        pos = substr.find(',');
        try
        {
          insertSizeMin = lexical_cast<int>(substr.substr(0, pos));
        } catch (bad_lexical_cast &)
        {
          cout << "Error at line " << count + 1 << " when reading file \""
              << fLoc << "\"." << endl;
          return 1;
        }

        substr = substr.substr(pos + 1);
        pos = substr.find(']');
        try
        {
          insertSizeMax = lexical_cast<int>(substr.substr(0, pos));
        } catch (bad_lexical_cast &)
        {
          cout << "Error at line " << count + 1 << " when reading file \""
              << fLoc << "\"." << endl;
          return 1;
        }
      }
      else
      {
        //cout << "\">\" sth: "  << str << endl;
      }
    }
    else
    {
      iss.str(str);
      int libNo = 0;
      int insertSize;
      int freq;

      substr.clear();
      iss >> substr;

      try
      {
        insertSize = lexical_cast<int>(substr);
      } catch (bad_lexical_cast &)
      {
        cout << "Error at line " << count + 1 << " when reading file \"" << fLoc
            << "\"." << endl;
        return 1;
      }
      while (true)
      {
        substr.clear();
        iss >> substr;

        if (substr.size() == 0)
          break;

        try
        {
          freq = lexical_cast<int>(substr);
        } catch (bad_lexical_cast &)
        {
          cout << "Error at line " << count + 1 << " when reading file \""
              << fLoc << "\"." << endl;
          return 1;
        }

        libVec[libNo].insertSizeFreqMap.insert(make_pair(insertSize, freq));

        libNo++;
      }

    }

    count++;
  }

  // Stop timer.
  gettimeofday(&t2, 0);

  // Show time.
  cout << endl << "Duration of loading insert size frequencies: ";
  cout << (t2.tv_sec - t1.tv_sec) + (t2.tv_usec - t1.tv_usec) / 1e6;
  cout << " s." << endl;

  delete line;

  return 0;
}

int StrInput::ReadLocus(
    deque<Region> & strLocDeque,
    const filesystem::path &f_path)
{
  string fLoc = f_path.string();
  ifstream locInput;
  timeval t1, t2;
  const int bufferSize = 256;
  char *line;
  int rowNo = 1;
  int len;
  size_t pos;

  locInput.open(fLoc.c_str(), ifstream::in);

  // If failing to open file.
  if ((void*) locInput == 0)
  {
    cout << "Failed to open file with locus coordinates: \"";
    cout << fLoc << "\"." << endl;
    return 1;
  }

  // Get length of file.
  locInput.seekg(0, ios::end);
  len = locInput.tellg();
  locInput.seekg(0, ios::beg);

  // Skip comment lines.
  while (true)
  {
    char c;
    locInput.get(c);

    if (c != '>')
    {
      pos = locInput.tellg();

      if (pos > 0)
      {
        pos--;
        locInput.seekg(pos);
      }

      break;
    }

    while (true)
    {
      char c;
      locInput.get(c);
      if (c == '\n' || c == '\0')
        break;
    }

    rowNo++;
  }

  // Create and initialize char vector.
  line = new char[bufferSize];

  for (int i0 = 0; i0 < bufferSize; i0++)
  {
    line[i0] = '\0';
  }

  // Start timer.
  gettimeofday(&t1, 0);

  // Read STR locus coordinates from file.
  while (!locInput.eof()) // && locInput.tellg() != len)
  {
    if (int(locInput.tellg()) == len)
      break;
    Region strLoc = { -1, "", -1, -1 };
    string str;
    istringstream is;

    locInput.getline(line, bufferSize);
    str.append(line);
    is.str(str);

    is >> strLoc.refIdStr;

    try
    {
      is >> str;
      strLoc.refBeg = lexical_cast<int>(str);
      str.clear();
      is >> str;
      strLoc.refEnd = lexical_cast<int>(str);
    } catch (bad_lexical_cast &)
    {
      cout << "Row " << rowNo << " in the file \"" << fLoc;
      cout << "\" does not contain valid coordinates of an STR locus." << endl;
      return 1;
    }

    if (strLoc.refIdStr.size() == 0 || strLoc.refBeg >= strLoc.refEnd
        || strLoc.refBeg < 0)
    {
      cout << "Row " << rowNo << " in the file \"" << fLoc;
      cout << "\" does not contain valid coordinates of an STR locus." << endl;
      return 1;
    }

    strLocDeque.push_back(strLoc);
    rowNo++;
  }

  // Stop timer.
  gettimeofday(&t2, 0);

  // Show time.
  cout << endl << "Duration of loading STR locus coordinates: ";
  cout << (t2.tv_sec - t1.tv_sec) + (t2.tv_usec - t1.tv_usec) / 1e6;
  cout << " s." << endl;

  // Delete heap memory.
  delete line;

  // Close files.
  locInput.close();

  return 0;
}

int StrInput::ReadLocusRefDon(
    deque<pair<Region, Region> > & strLocDeque,
    const filesystem::path &f_path)
{
  string fLoc = f_path.string();
  ifstream locInput;
  timeval t1, t2;
  const int bufferSize = 256;
  char *line;
  int rowNo = 1;
  int len;
  size_t pos;

  locInput.open(fLoc.c_str(), ifstream::in);

  // If failing to open file.
  if ((void*) locInput == 0)
  {
    cout << "Failed to open file with locus coordinates: \"";
    cout << fLoc << "\"." << endl;
    return 1;
  }

  // Get length of file.
  locInput.seekg(0, ios::end);
  len = locInput.tellg();
  locInput.seekg(0, ios::beg);

  // Skip comment lines.
  while (true)
  {
    char c;
    locInput.get(c);

    if (c != '>')
    {
      pos = locInput.tellg();

      if (pos > 0)
      {
        pos--;
        locInput.seekg(pos);
      }

      break;
    }

    while (true)
    {
      char c;
      locInput.get(c);
      if (c == '\n' || c == '\0')
        break;
    }

    rowNo++;
  }

  // Create and initialize char vector.
  line = new char[bufferSize];

  for (int i0 = 0; i0 < bufferSize; i0++)
  {
    line[i0] = '\0';
  }

  // Start timer.
  gettimeofday(&t1, 0);

  // Read STR locus coordinates from file.
  while (!locInput.eof()) // && locInput.tellg() != len)
  {
    if (int(locInput.tellg()) == len)
      break;
    Region strLocRef = { -1, "", -1, -1 };
    Region strLocDon = strLocRef;
    string str;
    istringstream is;

    locInput.getline(line, bufferSize);
    str.append(line);
    is.str(str);

    //is >> strLoc.refIdStr;

    try
    {
      is >> strLocRef.refIdStr;
      is >> str;
      strLocRef.refBeg = lexical_cast<int>(str);
      str.clear();
      is >> str;
      strLocRef.refEnd = lexical_cast<int>(str);
      strLocDon.refIdStr = strLocRef.refIdStr;
      is >> str;
      strLocDon.refBeg = lexical_cast<int>(str);
      str.clear();
      is >> str;
      strLocDon.refEnd = lexical_cast<int>(str);
    } catch (bad_lexical_cast &)
    {
      cout << "Row " << rowNo << " in the file \"" << fLoc;
      cout << "\" does not contain valid coordinates of an STR locus." << endl;
      return 1;
    }

    if (strLocRef.refIdStr.size() == 0 || strLocRef.refBeg >= strLocRef.refEnd
        || strLocRef.refBeg < 0)
    {
      cout << "Row " << rowNo << " in the file \"" << fLoc;
      cout << "\" does not contain valid coordinates of an STR locus." << endl;
      return 1;
    }

    strLocDeque.push_back(make_pair(strLocRef, strLocDon));
    rowNo++;
  }

  // Stop timer.
  gettimeofday(&t2, 0);

  // Show time.
  cout << endl << "Duration of loading STR locus coordinates: ";
  cout << (t2.tv_sec - t1.tv_sec) + (t2.tv_usec - t1.tv_usec) / 1e6;
  cout << " s." << endl;

  // Delete heap memory.
  delete line;

  // Close files.
  locInput.close();

  return 0;
}

int StrInput::ReadPriorPdf(
    list<T> & priorList,
    //const string & fPath,
    //const string & fName,
    const filesystem::path &f_path,
    int nLineSkip,
    int indelSizeMax,
    int indelStepSize)
{
  string fLoc = f_path.string(); //fPath + fName;
  ifstream priorInput;
  timeval t1, t2;
  char *line;
  int rowNo = 0;
  int len;
  T sum = 0;
  const int bufferSize = 256;
  int nElem = 2 * (indelSizeMax / indelStepSize) + 1;

  priorInput.open(fLoc.c_str(), fstream::in);

  if (!priorInput.is_open())
  {
    cout << "Failed to open file with prior probabilities: \"";
    cout << fLoc << "\"." << endl;
    return 1;
  }

  // Skip some lines.
  for (int i0 = 0; i0 < nLineSkip; i0++)
  {
    while (true)
    {
      char c;
      priorInput.get(c);
      if (c == '\n' || c == '\0')
        break;
    }

    rowNo++;
  }

  // Create and initialize char vector.
  line = new char[bufferSize];

  for (int i0 = 0; i0 < bufferSize; i0++)
  {
    line[i0] = '\0';
  }

  // Get length of file.
  priorInput.seekg(0, ios::end);
  len = size_t(priorInput.tellg());
  priorInput.seekg(0, ios::beg);

  // Start timer.
  gettimeofday(&t1, 0);

  // Read STR locus coordinates from file.
  while (!priorInput.eof()) // && priorInput.tellg() != len)
  {
    if (int(priorInput.tellg()) == len)
      break;

    string str;
    istringstream is;
    T val;

    priorInput.getline(line, bufferSize);
    str.append(line);
    is.str(str);

    try
    {
      is >> str;
      is >> str;
      val = lexical_cast<T>(str);
      str.clear();
    } catch (bad_lexical_cast &)
    {
      cout << "Row " << rowNo << " in the file \"" << fLoc;
      cout << "\" does not contain a prior probability." << endl;
      return 1;
    }

    if (val < 0)
    {
      cout << "All values in the file \"" << fLoc;
      cout << "\" must contain non-negative prior probabilities." << endl;
      return 1;
    }

    priorList.push_back(val);
    sum += val;
    rowNo++;
  }

  // Stop timer.
  gettimeofday(&t2, 0);

  // Show time.
  cout << endl << "Duration of loading prior probabilities: ";
  cout << (t2.tv_sec - t1.tv_sec) + (t2.tv_usec - t1.tv_usec) / 1e6;
  cout << " s." << endl;

  if (sum <= 0)
  {
    cout << "The sum of prior probabilities in the file \"" << fLoc;
    cout << "\" must be positive." << endl;
    return 1;
  }

  if (nElem != int(priorList.size()))
  {
    cout << "The number of elements in the file of prior probabilities  \""
        << fLoc;
    cout << "\" must be 2*(indel size maximum / indel step size) + 1." << endl;
    return 1;
  }

  delete line;

  return 0;
}

#endif /* STRINPUT_HPP_ */
