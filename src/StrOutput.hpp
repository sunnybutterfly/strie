#ifndef STROUTPUT_HPP_
#define STROUTPUT_HPP_

// STL libraries
#include <cmath>        // log(), exp()
#include <cstdlib>      // abs()
#include <fstream>      // open(), close()
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <queue>
#include <set>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

// External libraries
#include "boost_libraries.hpp"

// Own files
#include "type.hpp"
#include "StrLibrary.hpp"
#include "StrUtil.hpp"

// Namespace
using namespace boost;
using namespace std;

// Structure

// Class
class StrOutput
{
public:
  StrOutput();
  ~StrOutput();
  StrOutput(const StrOutput &rhs);

  int CloseLocOutput();

  string GetLocOutputFile();
  int GetPrintLibCount();
  int GetPrintNMax();
  int GetPrintProb();

  int IsOpenLocOutput();

  int SetPrecision(const string &location,
  //const string    &floatField,
      int precision,
      bool point);
  int SetPrecision(
      const char * location,
      int locationLen,
      int precision,
      bool point);
  int SetPrintLibCount(int printLibCountNew);
  int SetPrintNMax(int printNMaxNew);
  int SetPrintProb(int printProbNew);

protected:
  int CloseInsertSizeOutput();
  int CloseMpersPdfOutput();
  int CloseOverlapOutput();
  int ClosePriorPdfOutput();
  int CloseReadLenOutput();

  string GetInsertSizeOutputFile();
  string GetMpersPdfOutputFile();

  int IsOpenInsertSizeOutput();
  int IsOpenOverlapOutput();
  int IsOpenMpersPdfOutput();
  int IsOpenSimEvalOrOutput();
  int IsOpenSimEvalSrpOutput();

  int OpenInsertSizeOutput(const filesystem::path &f_path);
  int OpenLocOutput(
      const filesystem::path &f_path,
      const std::vector<StrLibrary> &libVec,
      int mapQMin);
  
  int OpenMpersPdfOutput(const filesystem::path &bam_fp);
  int OpenOverlapOutput(
      const filesystem::path &f_path,
      const std::vector<StrLibrary> &libVec);
  int OpenPriorPdfOutput(const filesystem::path &bam_fp);
  int OpenReadLenOutput(const filesystem::path &bam_fp);

  int WriteInsertSize(
      const std::vector<StrLibrary> &libVec,
      int mapQMpersMin,
      int insertSizeMin,
      int insertSizeMax);
  int WriteLocus(
      const Region &strLoc,
      const int strLen,
      const multi_array<T, 2> &indelLogl,
      const int indelStepSize,
      const bool isNan,
      const bool isWrongSize,
      const int nSrp,
      const T logOddsMin,
      const map<string*, int> &libCount,
      const map<string*, multiset<int> > &libFreq);
  int WriteLocusError(const int errNo, const Region &strLoc, const int strLen);
  int WriteMpersPdf(const multi_array<T, 2> &mpersPdf, const string &libTag);

  int WriteOverlap(
      const Region &strLoc,
      const int strLen,
      std::multimap<int, Read> &orMap,
      const int nOr,
      const map<string*, int> &libCount,
      const map<string*, multiset<int> > &libFreq);
  int WriteOverlapError(
      const int errNo,
      const Region &strLoc,
      const int strLen);
  int WritePriorPdf(const multi_array<T, 2> &priorPdf);
  int WriteReadLen(
      const multi_array<int, 2> &readLenFreqArr,
      const std::vector<StrLibrary> &libVec,
      int readLenMax,
      int readLenMapQMin);

  filesystem::path insertSizeOutput_path;
  ofstream insertSizeOutputF;

  filesystem::path locOutput_path;
  ofstream locOutputF;

  filesystem::path mpersPdfOutput_path;
  ofstream mpersPdfOutputF;

  filesystem::path overlapOutput_path;
  ofstream overlapOutputF;

  filesystem::path priorPdfOutput_path;
  ofstream priorPdfOutputF;

  filesystem::path readLenOutput_path;
  ofstream readLenOutputF;

  int printLibCount;
  int printNMax;
  int printProb;

private:
  StrUtil strUtil;    // Extra utilities.
};

StrOutput::StrOutput()
{

  printLibCount = 1;
  printNMax = 0;
  printProb = 0;

  //cout << "locOutputF.precision " << locOutputF.precision() << endl;
  //locOutputF.precision(4);
  //locOutputF.setf(0, ios::precision);
  SetPrecision("locus", 4, true);
  SetPrecision("readLen", 4, true);
  SetPrecision("overlap", 4, true);
}

StrOutput::~StrOutput()
{
  if (insertSizeOutputF.is_open() == true)
    insertSizeOutputF.close();

  if (mpersPdfOutputF.is_open() == true)
    mpersPdfOutputF.close();

  if (locOutputF.is_open() == true)
    locOutputF.close();
}

StrOutput::StrOutput(const StrOutput &rhs)
{
  this->printLibCount = rhs.printLibCount;
  this->printNMax = rhs.printNMax;
  this->printProb = rhs.printProb;
}

int StrOutput::CloseInsertSizeOutput()
{
  if (insertSizeOutputF.is_open() == false)
  {
    cout << "Insert size output file is already closed." << endl;
    return 1;
  }

  insertSizeOutputF.close();
  cout << "";

  return 0;
}

int StrOutput::CloseLocOutput()
{
  if (locOutputF.is_open() == false)
  {
    cout << "Locus output file is already closed." << endl;
    return 1;
  }

  locOutputF.close();

  return 0;
}

int StrOutput::CloseMpersPdfOutput()
{
  if (mpersPdfOutputF.is_open() == false)
  {
    cout << "Library output file is already closed." << endl;
    return 1;
  }

  mpersPdfOutputF.close();

  return 0;
}

int StrOutput::CloseOverlapOutput()
{
  if (overlapOutputF.is_open() == false)
  {
    cout << "Library output file is already closed." << endl;
    return 1;
  }

  overlapOutputF.close();

  return 0;
}

int StrOutput::ClosePriorPdfOutput()
{
  if (priorPdfOutputF.is_open() == false)
  {
    cout << "Prior output file is already closed." << endl;
    return 1;
  }

  priorPdfOutputF.close();

  return 0;
}

int StrOutput::CloseReadLenOutput()
{
  if (readLenOutputF.is_open() == false)
  {
    cout << "Read length output file is already closed." << endl;
    return 1;
  }

  readLenOutputF.close();

  return 0;
}

string StrOutput::GetInsertSizeOutputFile()
{
  return insertSizeOutput_path.string();
}

string StrOutput::GetMpersPdfOutputFile()
{
  return mpersPdfOutput_path.string();
}

string StrOutput::GetLocOutputFile()
{
  return locOutput_path.string();
}

int StrOutput::GetPrintLibCount()
{
  return printLibCount;
}

int StrOutput::GetPrintNMax()
{
  return printNMax;
}

int StrOutput::GetPrintProb()
{
  return printProb;
}

int StrOutput::IsOpenInsertSizeOutput()
{
  if (insertSizeOutputF.is_open() == true)
    return 1;

  return 0;
}

int StrOutput::IsOpenLocOutput()
{
  if (locOutputF.is_open() == true)
    return 1;

  return 0;
}

int StrOutput::IsOpenMpersPdfOutput()
{
  if (mpersPdfOutputF.is_open() == true)
    return 1;

  return 0;
}

int StrOutput::IsOpenOverlapOutput()
{
  if (overlapOutputF.is_open())
    return 1;

  return 0;
}

int StrOutput::OpenInsertSizeOutput(const filesystem::path &f_path)
{
  const string fName= f_path.string();

  if (insertSizeOutputF.is_open() == true)
  {
    cout << "Insert size output file is already open." << endl;
    return 1;
  }

  insertSizeOutputF.open(fName.c_str());

  if (insertSizeOutputF.is_open() == false)
  {
    cout << "Error opening insert size output file \"" << fName << "\"" << endl;
    return 2;
  }

  return 0;
}

int StrOutput::OpenLocOutput(
    const filesystem::path &f_path,
    const std::vector<StrLibrary> &libVec,
    int mapQMin)
{
  set<string>::const_iterator it1;
  set < string > libTagSet;
  string fName;
  char libChar = 'A';
  int pos;

  if (locOutputF.is_open() == true)
  {
    cout << "Locus output file is already open." << endl;
    return 1;
  }

  fName = f_path.string();
  locOutputF.open(fName.c_str());

  if (locOutputF.is_open() == false)
  {
    cout << "Error opening locus output file." << endl;
    return 2;
  }

  locOutputF << "> " << "Columns:" << " ";
  locOutputF << "CHROM" << " " << "BEGIN" << " " << "END" << " " << "STRLEN"
      << " " << "A1" << " " << "A2" << " " << "POSTPROB" << " " << "LOGODDS"
      << " " << "RPFREQ" << " " << "ISIZE" << endl;

  for (int i0 = 1; i0 < int(libVec.size()); i0++)
    libTagSet.insert(libVec[i0].libTag);

  locOutputF << "> Libraries: ";

  for (it1 = libTagSet.begin(); it1 != libTagSet.end(); it1++)
  {
    locOutputF << libChar << ":" << *it1 << ", ";
    libChar++;
  }

  pos = locOutputF.tellp();
  locOutputF.seekp(pos - 2);

  locOutputF << endl;
  locOutputF << "> " << "MAPQ minimum: " << mapQMin << endl;

  return 0;
}

int StrOutput::OpenMpersPdfOutput(
    const filesystem::path &f_path)
{
  string fName;

  if (mpersPdfOutputF.is_open() == true)
  {
    cout << "Library output file is already open." << endl;
    return 1;
  }

  fName = f_path.string(); //fPathNew + fNameNew;
  mpersPdfOutputF.open(fName.c_str());

  if (mpersPdfOutputF.is_open() == false)
  {
    cout << "Error opening library output file." << endl;
    return 2;
  }

  return 0;
}

int StrOutput::OpenOverlapOutput(
    const filesystem::path &f_path,
    const std::vector<StrLibrary> &libVec)
{
  set<string>::const_iterator it1;
  set < string > libTagSet;
  string fName;
  char libChar = 'A';
  int pos;

  if (overlapOutputF.is_open() == true)
  {
    cout << "Overlap output file is already open." << endl;
    return 1;
  }

  fName = f_path.string();
  overlapOutputF.open(fName.c_str());

  if (overlapOutputF.is_open() == false)
  {
    cout << "Error opening overlap output file." << endl;
    return 2;
  }

  overlapOutputF << "> " << "Columns:" << " ";
  overlapOutputF << "CHROM" << " " << "BEGIN" << " " << "END" << " " << "STRLEN"
      << " " << "LIB" << " " << "QNAME" << " " << "FLAG" << " " << "POS" << " "
      << "MPOS" << " " << "MAPQ" << " " << "RLEN" << " " << "ISIZE" << endl;

  for (int i0 = 1; i0 < int(libVec.size()); i0++)
    libTagSet.insert(libVec[i0].libTag);

  overlapOutputF << "> Libraries: ";

  for (it1 = libTagSet.begin(); it1 != libTagSet.end(); it1++)
  {
    overlapOutputF << libChar << ":" << *it1 << ", ";
    libChar++;
  }

  pos = overlapOutputF.tellp();
  overlapOutputF.seekp(pos - 2);
  overlapOutputF << endl;

  return 0;
}

int StrOutput::OpenPriorPdfOutput(const filesystem::path &f_path)
{
  string fLoc;

  if (priorPdfOutputF.is_open() == true)
  {
    cout << "Prior output file is already open." << endl;
    return 1;
  }

  fLoc = f_path.string();
  priorPdfOutputF.open(fLoc.c_str());

  if (priorPdfOutputF.is_open() == false)
  {
    cout << "Error opening prior output file." << endl;
    return 2;
  }

  return 0;
}

int StrOutput::OpenReadLenOutput(const filesystem::path &f_path)
{
  string fName;

  if (readLenOutputF.is_open() == true)
  {
    cout << "Read length output file is already open." << endl;
    return 1;
  }

  fName = f_path.string();
  readLenOutputF.open(fName.c_str());

  if (readLenOutputF.is_open() == false)
  {
    cout << "Error opening read length output file \"" << fName << "\"" << endl;
    return 2;
  }

  return 0;
}

int StrOutput::SetPrecision(const string &location,
    int precision,
    bool point)
{
  if (precision < 1)
  {
    cout << "Precision must be a positive integer." << endl;
    return 1;
  }

  if (location != "locus" && location != "library" && location != "readLen"
      && location != "overlap")
  {
    cout
        << "\"location\" must be \"locus\", \"library\", \"readLen\" or \"overlap\"."
        << endl;
    return 2;
  }

  string floatField = "fixed";

  if ((floatField == "scientific" || floatField == "fixed"
      || floatField == "variable") == false)
  {
    cout << "\"floatField\" must be \"scientific\", \"fixed\" or \"variable\"."
        << endl;
    return 3;
  }

  if (location == "locus")
  {
    locOutputF.precision(precision);

    if (floatField == "scientific")
      locOutputF.setf(ios::scientific, ios::floatfield);
    else if (floatField == "fixed")
      locOutputF.setf(ios::fixed, ios::floatfield);
    else if (floatField == "variable")
      locOutputF.unsetf(ios::floatfield);

    if (point)
      locOutputF.setf(ios::showpoint);
    else
      locOutputF.unsetf(ios::showpoint);
  }
  else if (location == "library")
  {
    mpersPdfOutputF.precision(precision);

    if (floatField == "scientific")
      mpersPdfOutputF.setf(ios::scientific, ios::floatfield);
    else if (floatField == "fixed")
      mpersPdfOutputF.setf(ios::fixed, ios::floatfield);
    else if (floatField == "variable")
      mpersPdfOutputF.unsetf(ios::floatfield);

    if (point)
      mpersPdfOutputF.setf(ios::showpoint);
    else
      mpersPdfOutputF.unsetf(ios::showpoint);
  }
  else if (location == "readLen")
  {
    readLenOutputF.precision(precision);

    if (floatField == "scientific")
      readLenOutputF.setf(ios::scientific, ios::floatfield);
    else if (floatField == "fixed")
      readLenOutputF.setf(ios::fixed, ios::floatfield);
    else if (floatField == "variable")
      readLenOutputF.unsetf(ios::floatfield);

    if (point)
      readLenOutputF.setf(ios::showpoint);
    else
      readLenOutputF.unsetf(ios::showpoint);
  }
  else if (location == "overlap")
  {
    overlapOutputF.precision(precision);

    if (floatField == "scientific")
      overlapOutputF.setf(ios::scientific, ios::floatfield);
    else if (floatField == "fixed")
      overlapOutputF.setf(ios::fixed, ios::floatfield);
    else if (floatField == "variable")
      overlapOutputF.unsetf(ios::floatfield);

    if (point)
      overlapOutputF.setf(ios::showpoint);
    else
      overlapOutputF.unsetf(ios::showpoint);
  }

  return 0;
}

int StrOutput::SetPrecision(
    const char * location,
    int locationLen,
    int precision,
    bool point)
{
  int result;
  string location2(location, locationLen);

  result = StrOutput::SetPrecision(location2, precision, point);

  return result;
}

int StrOutput::SetPrintLibCount(int printLibCountNew)
{
  if (printLibCountNew != 0 && printLibCountNew != 1)
  {
    cout << "Print library count must be 0 (false) or 1 (true)." << endl;
    return 1;
  }

  this->printLibCount = printLibCountNew;

  return 0;
}

int StrOutput::SetPrintNMax(int printNMaxNew)
{
  if (printNMaxNew < 0)
  {
    cout << "Print max values must be at least 0." << endl;
    return 1;
  }

  this->printNMax = printNMaxNew;

  return 0;
}

int StrOutput::SetPrintProb(int printProbNew)
{
  if (printProbNew != 0 && printProbNew != 1)
  {
    cout << "Print probability must be 0 (false) or 1 (true)." << endl;
    return 1;
  }

  this->printProb = printProbNew;

  return 0;
}

int StrOutput::WriteInsertSize(
    const std::vector<StrLibrary> &libVec,
    int mapQMpersMin,
    int insertSizeMin,
    int insertSizeMax)
{
  if (insertSizeOutputF.is_open() == false)
  {
    cout << "Insert size output is not open." << endl;
    return 1;
  }

  std::vector<StrLibrary>::const_iterator it0;
  map<int, int>::const_iterator it1;
  int pos;

  // Write library tags.
  insertSizeOutputF << "> Libraries: ";
  for (it0 = libVec.begin(); it0 != libVec.end(); it0++)
    insertSizeOutputF << it0->libTag << " ";

  pos = insertSizeOutputF.tellp();
  insertSizeOutputF.seekp(pos - 1);
  insertSizeOutputF << endl;

  // Write library statistics.
  for (int i0 = 0; i0 < 10; i0++)
  {
    insertSizeOutputF << "> ";
    switch (i0)
    {
    case 0:
      insertSizeOutputF << "No. of read-pairs: ";
      for (it0 = libVec.begin(); it0 != libVec.end(); it0++)
        insertSizeOutputF << it0->nRead << " ";
      break;
    case 1:
      insertSizeOutputF << "p0: ";
      for (it0 = libVec.begin(); it0 != libVec.end(); it0++)
        insertSizeOutputF << it0->pct0 << " ";
      break;
    case 2:
      insertSizeOutputF << "p2.5: ";
      for (it0 = libVec.begin(); it0 != libVec.end(); it0++)
        insertSizeOutputF << it0->pct2p5 << " ";
      break;
    case 3:
      insertSizeOutputF << "p25: ";
      for (it0 = libVec.begin(); it0 != libVec.end(); it0++)
        insertSizeOutputF << it0->pct25 << " ";
      break;
    case 4:
      insertSizeOutputF << "p50: ";
      for (it0 = libVec.begin(); it0 != libVec.end(); it0++)
        insertSizeOutputF << it0->pct50 << " ";
      break;
    case 5:
      insertSizeOutputF << "p75: ";
      for (it0 = libVec.begin(); it0 != libVec.end(); it0++)
        insertSizeOutputF << it0->pct75 << " ";
      break;
    case 6:
      insertSizeOutputF << "p97.5: ";
      for (it0 = libVec.begin(); it0 != libVec.end(); it0++)
        insertSizeOutputF << it0->pct97p5 << " ";
      break;
    case 7:
      insertSizeOutputF << "p100: ";
      for (it0 = libVec.begin(); it0 != libVec.end(); it0++)
        insertSizeOutputF << it0->pct100 << " ";
      break;
    case 8:
      insertSizeOutputF << "Mean: ";
      for (it0 = libVec.begin(); it0 != libVec.end(); it0++)
        insertSizeOutputF << it0->mean << " ";
      break;
    case 9:
      insertSizeOutputF << "Std: ";
      for (it0 = libVec.begin(); it0 != libVec.end(); it0++)
        insertSizeOutputF << it0->stDev << " ";
      break;
    }

    pos = insertSizeOutputF.tellp();
    insertSizeOutputF.seekp(pos - 1);
    insertSizeOutputF << endl;
  }

  // Write other information.
  insertSizeOutputF << "> " << "MAPQ minimum: " << mapQMpersMin << endl;
  insertSizeOutputF << "> " << "Insert size interval: [" << insertSizeMin << ","
      << insertSizeMax << "]" << endl;

  //pos = insertSizeOutputF.tellp();
  //insertSizeOutputF.seekp(pos - 1);

  // Write insert sizes and frequencies of libraries.
  //insertSizeOutputF << endl;

  for (int i0 = insertSizeMin; i0 < insertSizeMax + 1; i0++)
  {
    insertSizeOutputF << i0 << ' ';

    for (it0 = libVec.begin(); it0 != libVec.end(); it0++)
    {
      it1 = it0->insertSizeFreqMap.find(i0);

      if (it1 == it0->insertSizeFreqMap.end())
        insertSizeOutputF << 0;
      else
        insertSizeOutputF << it1->second;

      insertSizeOutputF << ' ';
    }

    pos = insertSizeOutputF.tellp();
    insertSizeOutputF.seekp(pos - 1);
    insertSizeOutputF << endl;
  }

  return 0;
}

// Write information about locus to file. It is assumed that the file object locOutputF has
// a file open already and that this has been checked by the calling function.
int StrOutput::WriteLocus(
    const Region &strLoc,
    const int strLen,
    const multi_array<T, 2> &indelLogl,
    const int indelStepSize,
    const bool isNan,
    const bool isWrongSize,
    const int nSrp,
    const T logOddsMin,
    const map<string*, int> &libCount,
    const map<string*, multiset<int> > &libFreq)
{

  map<T, pair<int, int> >::iterator it1;
  map<T, pair<int, int> >::reverse_iterator rit1;
  map<string, int>::const_iterator it2;
  map<string*, int>::const_iterator it3;
  map<string*, multiset<int> >::const_iterator it4;
  multiset<int>::const_iterator it5;

  map<T, pair<int, int> > probMax;

  //string mainTag = "all";
  int nRow = indelLogl.size();
  int center = (indelLogl.size() - 1) / 2;
  int valueMid = (nRow - 1) * indelStepSize / 2;
  T logOdds;

  // Find the largest probabilities.
  for (int i0 = 0; i0 < printNMax + 1; i0++)
    probMax.insert(make_pair(T(-i0 - 1), make_pair(-1, -1)));

  it1 = probMax.begin();

  for (int i0 = 0; i0 < int(indelLogl.shape()[0]); i0++)
  {
    for (int i1 = 0; i1 < i0 + 1; i1++)
    {
      T val = indelLogl[i0][i1];

      if (it1->first < val)
      {
        probMax.erase(it1);
        probMax.insert(
            make_pair(val,
                make_pair(i0 * indelStepSize - valueMid,
                    i1 * indelStepSize - valueMid)));
        it1 = probMax.begin();
      }
    }
  }

  it1 = probMax.end();
  it1--;

  // Log-odds is undefined if there are no spanning reads.
  if (isNan == true)
    logOdds = numeric_limits < T > ::signaling_NaN();
  else
    logOdds = log(it1->first / indelLogl[center][center]);

  // If the log-odds ratio is smaller than the minimum, report an error.
  // (If log-odds is NaN, the if statement evaluates to false.)
  if (logOdds < logOddsMin)
  {
    WriteLocusError(1, strLoc, strLen);
    return 1;
  }

  // Print output of things that must be printed.
  locOutputF << strLoc.refIdStr << "\t";
  locOutputF << strLoc.refBeg << "\t";
  locOutputF << strLoc.refEnd << "\t";
  locOutputF << strLen << "\t";
  locOutputF << it1->second.second << "\t";
  locOutputF << it1->second.first << "\t";
  locOutputF << it1->first << "\t";
  locOutputF << logOdds << "\t";
  locOutputF << nSrp;
  //locOutputF << endl;

  // Print optional information.
  if (isNan == false)
  {
    // Print library count of spanning read-pairs.
    if (this->GetPrintLibCount() == true)
    {
      string libInsertSize;
      int pos;
      char libChar = 'A';

      libInsertSize.reserve(8 * nSrp);
      locOutputF << "\t(";
      libInsertSize += "(";

      // If MPERS distribution computed using individual libraries was used.

      for (it3 = libCount.begin(), it4 = libFreq.begin(); it3 != libCount.end();
          it3++, it4++)
      {
        if (it3->second)
        {
          locOutputF << it3->second << libChar << ',';

          libInsertSize.push_back(libChar);
          libInsertSize.push_back(':');

          for (it5 = it4->second.begin(); it5 != it4->second.end(); it5++)
          {
            libInsertSize += lexical_cast<string>(*it5) + ",";
          }

          pos = locOutputF.tellp();
          locOutputF.seekp(pos - 1);
          locOutputF.write(",", 1);

          pos = libInsertSize.size();
          libInsertSize.replace(pos - 1, 1, ";");
        }

        libChar++;
      }

      pos = locOutputF.tellp();
      locOutputF.seekp(pos - 1);
      locOutputF.write(")", 1);

      //pos = ss.tellg();
      //ss.seekg(pos-1);
      //ss << ')';
      //locOutputF << ss.str();
      //locOutputF.write (")", 1);

      pos = libInsertSize.size();
      libInsertSize.replace(pos - 1, 1, ")");
      locOutputF << "\t" << libInsertSize;
    }

    // Print posterior probabilities.
    if (this->GetPrintProb() == true)
    {
      locOutputF << endl;

      for (int i0 = int(indelLogl.shape()[0]) - 1; i0 > -1; i0--)
      {
        for (int i1 = 0; i1 < i0; i1++)
          locOutputF << indelLogl[i0][i1] << " ";

        locOutputF << indelLogl[i0][i0];
      }
    }
  }

  locOutputF << endl;

  return 0;
}

int StrOutput::WriteLocusError(
    const int errNo,
    const Region &strLoc,
    const int strLen)
{
  string errStr = "";

  switch (errNo)
  {
  case 1:
    errStr = "Err 1: Log-odds ratio too small.";
    break;
  case 2:
    errStr = "Err 2: No reads.";
    break;
  case 3:
    errStr = "Err 3: No spanning read-pairs.";
    break;
  case 4:
    errStr = "Err 4: Too few spanning read-pairs.";
    break;
  case 5:
    errStr = "Err 5: STR too large or too small.";
    break;
  }

  if (errStr.empty() == true)
  {
    cout << "Wrong error no." << endl;
    return 1;
  }

  // Print output of things that must be printed.
  locOutputF << strLoc.refIdStr << "\t";
  locOutputF << strLoc.refBeg << "\t";
  locOutputF << strLoc.refEnd << "\t";
  locOutputF << strLen << "\t";
  locOutputF << "# " << errStr << endl;

  return 0;
}

int StrOutput::WriteMpersPdf(
    const multi_array<T, 2> &mpersPdf,
    const string &libTag)
{
  if (mpersPdfOutputF.is_open() == false)
  {
    cout << "Library output is not open." << endl;
    return 1;
  }

  int m = int(mpersPdf.shape()[0]);
  int n = int(mpersPdf.shape()[1]);

  // Write library tag.
  mpersPdfOutputF << libTag << endl;

  for (int i0 = 0; i0 < m; i0++)
  {
    for (int i1 = 0; i1 < n - 1; i1++)
      mpersPdfOutputF << mpersPdf[i0][i1] << " ";

    mpersPdfOutputF << mpersPdf[i0][n - 1] << endl;
  }

  return 0;
}

int StrOutput::WriteOverlap(
    const Region &strLoc,
    const int strLen,
    std::multimap<int, Read> &orMap,
    const int nOr,
    const map<string*, int> &libCount,
    const map<string*, multiset<int> > &libFreq)
{
  std::multimap<int, Read>::iterator cit0;
  std::multimap<int, int>::iterator it1;
  std::multimap<int, int>::reverse_iterator rit2;
  multimap<int, int> indelSizeFreq;
  int nReadValid = 0;
  int pos;

  for (cit0 = orMap.begin(); cit0 != orMap.end(); cit0++)
  {
    if (cit0->second.overlapType == 0)
    {
      int overlapIndelSize = cit0->second.overlapIndelSize;
      bool found = false;

      // Look for indel sizes in STR locus.
      for (it1 = indelSizeFreq.begin(); it1 != indelSizeFreq.end(); it1++)
      {
        if (it1->second == overlapIndelSize)
        {
          pair<int, int> tmp = make_pair(it1->first + 1, it1->second);
          indelSizeFreq.erase(it1);
          indelSizeFreq.insert(tmp);
          nReadValid++;
          found = true;
          break;
        }
      }

      if (found == false)
      {
        indelSizeFreq.insert(make_pair(1, overlapIndelSize));
        nReadValid++;
      }
    }
  }

  // Find error.
  if (indelSizeFreq.empty())
  {
    WriteOverlapError(4, strLoc, strLen);
    return 0;
  }

  // Print information about overlapping reads.
  // Locus information.
  overlapOutputF << strLoc.refIdStr << " ";
  overlapOutputF << strLoc.refBeg << " ";
  overlapOutputF << strLoc.refEnd << " ";
  overlapOutputF << strLen << " ";
  overlapOutputF << orMap.size() << " ";
  //overlapOutputF << nReadValid;

  // Indel information in STR locus.
  // (frequencies, probabilities, standard deviations)

  for (rit2 = indelSizeFreq.rbegin(); rit2 != indelSizeFreq.rend(); rit2++)
  {
    T prob = rit2->first / T(nReadValid);
    overlapOutputF << "{" << rit2->second << "}:"; // indel size
    overlapOutputF << rit2->first << ","; // frequency
    overlapOutputF << fixed;
    overlapOutputF << setprecision(4) << prob << " "; // probability
    //overlapOutputF << setprecision(4) << pow(prob*(1-prob),0.5) << ";"; // standard deviation
  }

  pos = overlapOutputF.tellp();
  overlapOutputF.seekp(pos - 1);
  //overlapOutputF.write(")", 1);
  overlapOutputF << endl;

  return 0;
}

int StrOutput::WriteOverlapError(
    const int errNo,
    const Region &strLoc,
    const int strLen)
{
  string errStr;

  switch (errNo)
  {
  case 1:
    errStr = "Err 1: STR too large or too small.";
    break;
  case 2:
    errStr = "Err 2: No reads.";
    break;
  case 3:
    errStr = "Err 3: No overlapping reads.";
    break;
  case 4:
    errStr = "Err 4: No valid overlapping reads.";
    break;

  }

  if (errStr.empty() == true)
  {
    cout << "Wrong error no." << endl;
    return 1;
  }

  // Print output of things that must be printed.
  overlapOutputF << strLoc.refIdStr << " ";
  overlapOutputF << strLoc.refBeg << " ";
  overlapOutputF << strLoc.refEnd << " ";
  overlapOutputF << strLen << " ";
  overlapOutputF << "# " << errStr << endl;

  return 0;
}

int StrOutput::WritePriorPdf(const multi_array<T, 2> &priorPdf)
{
  if (priorPdfOutputF.is_open() == false)
  {
    cout << "Library output is not open." << endl;
    return 1;
  }

  int m = int(priorPdf.shape()[0]);

  // Write library tag.

  for (int i0 = 0; i0 < m; i0++)
  {
    for (int i1 = 0; i1 < i0; i1++)
      priorPdfOutputF << priorPdf[i0][i1] << " ";

    priorPdfOutputF << priorPdf[i0][i0] << endl;
  }

  return 0;
}

int StrOutput::WriteReadLen(const multi_array<int, 2> &readLenFreqArr,
    const std::vector<StrLibrary> &libVec,
    int readLenMax,
    int readLenMapQMin)
{
  string func = "SO::WriteReadLen ";

  if (readLenOutputF.is_open() == false)
  {
    cout << "Insert size output is not open." << endl;
    return 1;
  }

  std::vector<StrLibrary>::const_iterator cit0;
  map<pair<int, int>, int>::const_iterator cit1;
  map<int, int>::const_iterator cit2;
  int nLib = int(libVec.size());
  int nCol = 0;
  multi_array<int, 2> readLenPercArr;
  vector<T> readLenNRead;
  vector<T> readLenMeanArr;
  vector<T> readLenStdArr;
  //int colNoPrev = 0;
  int colNo = -1;

  T percArr[] = { 0, 2.5, 25, 50, 75, 97.5, 100 };
  int percArrLen = 7;

  for (int i0 = 0; i0 < nLib; i0++)
  {
    nCol += int(libVec[i0].rLenNRead.size());

    for (cit2 = libVec[i0].rLenNRead.begin();
        cit2 != libVec[i0].rLenNRead.end(); cit2++)
      readLenNRead.push_back(cit2->second);
  }

  readLenPercArr.resize(extents[readLenMax + 1][nCol]);
  readLenMeanArr.resize(nCol);
  readLenStdArr.resize(nCol);

  for (int i0 = 0; i0 < nLib; i0++)
  {
    int rLenPrev = -1;
    //colNo = colNoPrev;

    for (cit1 = libVec[i0].overlapFfMap.begin();
        cit1 != libVec[i0].overlapFfMap.end(); cit1++)
    {
      int rLen = cit1->first.second;
      int freq = cit1->second;

      if (rLenPrev != cit1->first.first)
      {
        rLenPrev = cit1->first.first;
        colNo++;
      }

      readLenPercArr[rLen][colNo] = freq;
    }
  }

  // Write library tags.
  readLenOutputF << "> Libraries:";
  for (cit0 = libVec.begin(); cit0 != libVec.end(); cit0++)
  {

    for (cit2 = cit0->rLenNRead.begin(); cit2 != cit0->rLenNRead.end(); cit2++)
    {
      int rLen = cit2->first;
      readLenOutputF << " " << cit0->libTag << ":" << rLen;
    }
    //readLenOutputF << " " << it0->libTag;
  }

  readLenOutputF << endl;

  // Write sum of reads per library.
  readLenOutputF << "> " << "No. of reads:";

  for (int i0 = 0; i0 < nLib; i0++)
  {
    int sum = 0;

    for (cit2 = libVec[i0].rLenNRead.begin();
        cit2 != libVec[i0].rLenNRead.end(); cit2++)
    {
      sum = cit2->second;
      readLenOutputF << " " << sum;
    }
  }

  readLenOutputF << endl;

  for (int i0 = 0; i0 < percArrLen; i0++)
  {
    T perc = percArr[i0] / 100;
    readLenOutputF << "> p" << perc << ":";

    for (int i1 = 0; i1 < nCol; i1++)
    {
      int nRead = 0;
      T nReadSum = readLenNRead[i1];

      for (int i2 = 0; i2 < readLenMax + 1; i2++)
      {
        nRead += readLenPercArr[i2][i1];
        T perc_now = (T(nRead) + 0.5) / nReadSum;

        if (perc_now >= perc && nRead > 0)
        {
          readLenOutputF << " " << i2;
          break;
        }
      }
    }

    readLenOutputF << endl;
  }

  // Write mean and stdev.

  readLenOutputF << "> Mean:";

  for (int i0 = 0; i0 < nCol; i0++)
  {
    uint64_t sum = 0;

    for (int i1 = 0; i1 < readLenMax + 1; i1++)
    {
      sum += i1 * readLenPercArr[i1][i0];
    }

    readLenMeanArr[i0] = sum / T(readLenNRead[i0]);
    readLenOutputF << " " << readLenMeanArr[i0];
  }

  readLenOutputF << "\n> Std:";

  for (int i0 = 0; i0 < nCol; i0++)
  {
    T sum = 0;

    for (int i1 = 0; i1 < readLenMax; i1++)
    {
      T diff = (i1 - readLenMeanArr[i0]);
      sum += diff * diff * readLenPercArr[i1][i0];
    }

    readLenStdArr[i0] = pow(sum / (readLenNRead[i0] - 1), 0.5);
    readLenOutputF << " " << readLenStdArr[i0];
  }

  readLenOutputF << endl;

  // Write other information.
  readLenOutputF << "> " << "MAPQ minimum: " << readLenMapQMin << endl;
  readLenOutputF << "> " << "Read length interval: [" << 1 << "," << readLenMax
      << "]" << endl;

  // Write read length sizes and frequencies of libraries.
  for (int i0 = 1; i0 < readLenMax + 1; i0++)
  {
    readLenOutputF << i0;

    for (int i1 = 0; i1 < nCol; i1++)
    {
      readLenOutputF << " " << readLenPercArr[i0][i1];
    }

    readLenOutputF << endl;
  }

  return 0;
}

#endif /* STROUTPUT_HPP_ */
