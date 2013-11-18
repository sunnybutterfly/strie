#ifndef STRINDELSIZE_HPP_
#define STRINDELSIZE_HPP_

// STL libraries
#include <algorithm>    // swap()
#include <cmath>        // log(), exp()
#include <cassert>      // assert()
#include <cstdlib>      // abs(), rand()
#include <cstring>      // abs(), rand()
#include <deque>
#include <fstream>      // open(), close()
#include <iostream>
#include <limits>
#include <list>
#include <map>
#include <queue>
#include <set>
#include <memory>     // unique_ptr
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

// External libraries
#include "boost_libraries.hpp"
#include "api/BamAlignment.h"
#include "api/BamAux.h"
#include "api/BamIndex.h"
#include "api/BamReader.h"
#include "api/SamHeader.h"
#include <sys/time.h>   // gettimeofday

// Own files
#include "type.hpp"
#include "StrCalc.hpp"
#include "StrInput.hpp"
#include "StrLibrary.hpp"
#include "StrOutput.hpp"
#include "StrPd.hpp"
#include "StrUtil.hpp"

// Namespace
using namespace std;
using namespace BamTools;
using namespace boost;
using namespace boost::numeric::ublas;

// Global files.

// Class
class StrIndelSize: public StrInput,
    public StrOutput,
    public StrPd,
    public StrCalc
{
public:
  StrIndelSize();
  ~StrIndelSize();
  StrIndelSize(const StrIndelSize &rhs);

  int CalcMpersPdf();

  int ClearInsertSizeMap();
  int ClearLib();
  int ClearStrLoc();

  int CloseBam();

  int DelLib(string libOld);
  int DelLib(const char * libOld, int libOldLen);

  int GetAnchorLen();
  int GetBamFileLoc();
  int GetBufferSize();
  string GetCiEstimator();
  bool GetCheckNeighbor();
  bool GetClipping();
  int GetInsertSize(int &insertSizeMin, int &insertSizeMax);
  int GetInsertSizeMaxGrand();
  int GetLibInsertSize(int &libInsertSizeMax, int libNo);
  int GetLibInsertSizeFreq(int * insertSizeArr, int libNo);
  std::vector<string> GetLibPicked();
  int GetLibPicked(char *&libPicked, int &libPickedLen);
  int GetLibSizeMin();
  int GetLibTag(char *&libTag, int &libTagLen, int libNo);
  T GetLogOddsMin();
  int GetMapQMin();
  int GetMapQMpersMin();
  int GetMpersPdf(
      T * mpersPdf,
      const int nRowMpersPdf,
      const int nColMpersPdf,
      const int libNo);
  int GetNLib();
  int GetNSrpMin();
  int GetNStrLoc();
  int GetNVarIter();
  int GetReadCountMax();
  int GetRefSize(int &refSizeMin, int &refSizeMax);
  int GetRefCount();
  int GetRefData(string &refName, int refNo);
  int GetRefDataFull(char *& refName, int &refNameLen, int refNo);
  int GetStrLoc(string &refIdStr, int &refBeg, int &refEnd, int no);

  int InsertLib(string libNew);
  int InsertLib(const char * libNew, int libNewLen);

  int IsOpenBam();

  int LaunchStrEstimation();
  int LoadLibInsertSizeFreq(
      const string &libTagNew,
      std::map<int, int> &insertSizeMapNew);
  int LoadLibInsertSizeFreq(
      const char * libTagNew,
      int libTagNewLen,
      const int * insertSizeArr,
      const int * freqArr,
      int freqArrLen);
  int LoadStrLoc(string refIdStr, int refBeg, int refEnd);
  int LoadStrLoc(
      const char * refIdStr,
      int refIdStrLen,
      int refBeg,
      int refEnd);

  int OpenBam(const filesystem::path &bam_fp);
  int OpenLocOutput(const filesystem::path &f_path);

  int PickLib(string pickedLibNew);
  int PickLib(const char * pickedLibNew, int pickedLibNewLen);
  int PrintLibInfo(int type);

  int ReadBamHeader();
  int ReadInsertSize(const filesystem::path &f_path);
  int ReadInsertSize(const string &refIdStr, int refBeg, int refEnd);
  int ReadInsertSize(
      const char * refIdStr,
      int refIdStrLen,
      int refBeg,
      int refEnd);
  int ReadLocus(const filesystem::path &f_path);
  int ReadPrior(
      const filesystem::path &f_path,
      int nLineSkip,
      int indelSizeMaxNew,
      int indelStepSizeNew);
  int ReadReadLen(
      const filesystem::path &f_path,
      const Region &region,
      int readLenMax,
      int readLenMapQMin,
      int nReadMax);

  int SetAnchorLen(int anchorLenNew);
  int SetBufferSize(int bufferSizeNew);
  int SetCiEstimator(const string &ciEstimator);
  int SetClipping(bool clippingNew);
  int SetCheckNeighbor(bool checkNeighborNew);
  int SetInsertSize(int insertSizeMinNew, int insertSizeMaxNew);
  int SetLibInsertSize(int libInsertSizeMaxNew, const string &libTag);
  int SetLibInsertSize(
      int libInsertSizeMaxNew,
      const char * libTag,
      int libTagLen);
  int SetLibInsertSizeFreq(
      const int * insertSizeArrNew,
      const int arrLen,
      const int libNo);
  int SetLibSizeMin(int libSizeMinNew);
  int SetLogOddsMin(T logOddsMinNew);
  int SetMapQMin(int mapQMinNew);
  int SetMapQMpersMin(int mapQMpersMinNew);
  int SetMpersPdf(
      const char * libTag,
      int libTagLen,
      const T * mpersPdfNew,
      int m,
      int n);
  int SetNSrpMin(int nSrpNew);
  int SetNVarIter(int nVarIter);
  int SetReadCountMax(int readCountMaxNew);
  int SetRefSize(int refSizeMinNew, int refSizeMaxNew);

  int UnpickLib(const string &unpickedLibNew);
  int UnpickLib(const char * unpickedLibNew, int unpickedLibNewLen);

  int WriteInsertSize(const filesystem::path &f_path);
  int WriteMpersPdf(const filesystem::path &f_path);
  int WritePriorPdf(const filesystem::path &f_path);

private:
  T CalcDuration(timeval t1, timeval t2);
  int CalcLibInsertSize();
  int CalcVar(
      multi_array<T, 2> &indelLoglVar,
      const std::vector<multi_array<T, 2> > &indelLogl,
      int nIter);

  int ClearSrp(std::vector<Read> &readTo, const std::list<Read> &readFrom);

  int CompMean(
      multi_array<T, 2> &libMat,
      const Read &read,
      const int beg,
      const int end);

  int EstStrSize(
      const std::vector<Region> &strLocVec,
      const map<string, pair<std::vector<int>, std::vector<int> > > &neighborMap);

  int FindRead(
      std::list<Read> &readTo,
      const std::vector<Read> &readFrom,
      const Region strLoc);
  int FindSrp(std::list<Read> &readTo,
      std::list<Read> &readFrom,
      const Region strLoc);

  int GenIdLibPosMap();
  int GenLibTagLibPosMap();
  int GetAlignmentNext(Read &read);

  int Percentile(
      std::vector<T> &percentile,
      std::map<int, int> &samp,
      std::vector<T> &percentage,
      string libTag);

  int ScanRegion(std::list<Read> &readVec, const Region strLoc);
  int ScanStrLoc(std::vector<Read> &readVec, const Region strLoc);
  int SetRegion(Region region);

  int UpdateInsertSizeMaxGrand();

  multimap<string, string> libIdMultiMap; // Library tags (key) and read-group identifiers (values).
  map<string, string> idLibMap; // Read-group identifiers (keys) and library tags (values).
  map<string, int> idLibPosMap; // Library ID (keys) and library positions (values).
  map<string, int> libTagLibPosMap; // Library tag (keys) and library positions (values).
  set<string> libSet;             // Container of library tags.

  std::vector<Read> readList;   // Container of reads.
  multi_array<T, 2> indelLogl;  // Array of indel log-likelihood values.
  std::vector<StrLibrary> libVec;         // Container of libraries.
  list<Read> srp;        // Container of spanning read-pairs.
  deque<Region> strLocDeque; // Container of coordinates to visit.

  StrUtil strUtil;    // Extra utilities.

  filesystem::path bam_path;

  string refIdStrLatest;  // Latest BAM reference of searching.
  int begLatest;          // Latest beginning of BAM file region searched.
  int endLatest;          // Latest end of BAM file region searched.

  int anchorLen;      // Anchor length.
  int bufferSize;     // Buffer size.
  bool checkNeighbor; // Check whether read-pairs span more than one locus.
  bool clipping;      // Remove clipping when in case of STR overlap.
  int insertSizeMin;  // Minimum insert size.
  int insertSizeMax;  // Maximum insert size.
  int insertSizeMaxGrand;  // Maximum insert size of all libraries.
  int libSizeMin;     // Minimum library size.
  T logOddsMin;     // Minimum log-odds for STR locus to be reported.
  int mapQMin;  // Minimum MAPQ score for read-pairs in estimating indel length.
  int mapQMpersMin; // Minimum MAPQ score for read-pairs in computing MPERS distributions from
                    // insert size frequencies.
  int nSrpMin; // Minimum no. of spanning read-pairs for STR locus to be reported.
  int nStep;          // No. of indel steps.
  int nVarIter;       // Calculate variance.
  int readCountMax; // Maximum no. of reads to sample for computing MPERS distributions.
  int refSizeMin; // Minimum locus reference size allowed for computing inferred indel sizes.
  int refSizeMax; // Maximum locus reference size allowed for computing inferred indel sizes.
  int refBeg_now;

  unique_ptr<BamReader> bamReader;
};

StrIndelSize::StrIndelSize()
{
  bamReader.reset(new BamReader());

  strUtil = StrUtil();

  refIdStrLatest = "";
  begLatest = 0;
  endLatest = 0;

  anchorLen = 1;
  bufferSize = 100;
  checkNeighbor = true;
  clipping = true;
  insertSizeMin = 36;
  insertSizeMax = 450;
  insertSizeMaxGrand = insertSizeMax;
  libSizeMin = 5000;
  logOddsMin = 0;
  mapQMin = 30;
  mapQMpersMin = 30;
  nSrpMin = 1;
  nVarIter = 0;
  readCountMax = 1000000;
  refSizeMin = 15;
  refSizeMax = 250;

  indelLogl.resize(extents[0][0]);

  libVec.resize(1);
  libVec[0].libTag = "all";

  CalcMpersPdf();

}

StrIndelSize::~StrIndelSize()
{
  bamReader->Close();
}

StrIndelSize::StrIndelSize(const StrIndelSize &rhs)
{

}

T StrIndelSize::CalcDuration(timeval t1, timeval t2)
{
  return (t2.tv_sec + t2.tv_usec / 1e6) - (t1.tv_sec + t1.tv_usec / 1e6);
}

int StrIndelSize::CalcLibInsertSize()
{
  map<int, int>::iterator it0;
  map<int, int>::iterator it1;
  std::vector<StrLibrary>::iterator it2;
  int nRead = 0;

  // Count the aggregate no. of reads for each insert size.
  libVec[0].insertSizeFreqMap.clear();

  for (int i0 = 1; i0 < int(libVec.size()); i0++)
  {
    int nRead = 0;
    it2 = libVec.begin() + i0;

    for (it0 = it2->insertSizeFreqMap.begin();
        it0 != it2->insertSizeFreqMap.end(); it0++)
    {
      it1 = libVec[0].insertSizeFreqMap.find(it0->first);

      if (it1 == libVec[0].insertSizeFreqMap.end())
        libVec[0].insertSizeFreqMap.insert(make_pair(it0->first, it0->second));
      else
        it1->second += it0->second;

      nRead += it0->second;
    }

    libVec[i0].nRead = nRead;
  }

  // Count the aggregate no. of reads for each library.
  for (it0 = libVec[0].insertSizeFreqMap.begin();
      it0 != libVec[0].insertSizeFreqMap.end(); it0++)
    nRead += it0->second;

  libVec[0].nRead = nRead;

  return 0;
}

int StrIndelSize::CalcMpersPdf()
{
  std::vector<int>::iterator it1;
  std::vector<T> pcg;

  int VPosStart = insertSizeMin;
  int libInsertSizeMin = insertSizeMin;

  // Calculate MPERS probability distribution.
  for (int i0 = 0; i0 < int(libVec.size()); i0++)
  {
    map<int, int>::iterator it2;
    map<int, int> * insertSizeMap_tmp = &libVec[i0].insertSizeFreqMap;
    map<int, int> insertSizeMapRobust;
    multi_array<T, 2> * mpersPdf_tmp = &libVec[i0].mpersPdfArr;

    int libInsertSizeMax = libVec[i0].libInsertSizeMax;

    // Copy counts inside of interval [insertSizeMin, insertSizeMax].
    for (int i1 = libInsertSizeMin; i1 < libInsertSizeMax + 1; i1++)
    {
      int count;

      it2 = insertSizeMap_tmp->find(i1);
      count = it2 == insertSizeMap_tmp->end() ? 0 : it2->second;
      insertSizeMapRobust[i1] = count + 1;
    }

    // Resize MPERS probability distribution array.
    mpersPdf_tmp->resize(
        extents[refSizeMax - refSizeMin + 1][libInsertSizeMax - libInsertSizeMin
            + 1]);

    // Calculate MPERS probability distribution.
    for (int i1 = refSizeMin; i1 < refSizeMax + 1; i1++)
    {
      map<int, T>::iterator it3;

      map<int, T> probTemp1;
      map<int, T> probTemp2;

      int refSize = i1;
      int cutoffDist = refSize + 2 * anchorLen;
      int HPos = refSize - refSizeMin;
      int count = 0;
      long double normalSum = 0;

      for (it2 = insertSizeMapRobust.begin(); it2 != insertSizeMapRobust.end();
          it2++)
        if (it2->first >= cutoffDist)
          count += it2->second;

      for (it2 = insertSizeMapRobust.begin(); it2 != insertSizeMapRobust.end();
          it2++)
      {
        if (it2->first >= cutoffDist)
        {
          probTemp1[it2->first] = T(it2->second) / count;
        }
      }

      for (it3 = probTemp1.begin(); it3 != probTemp1.end(); it3++)
      {
        int p_m = it3->first - cutoffDist + 1;
        T lenProb = T(p_m * it3->second);

        probTemp2[it3->first] = lenProb;
        normalSum += lenProb;
      }

      for (it3 = probTemp2.begin(); it3 != probTemp2.end(); it3++)
      {
        (*mpersPdf_tmp)[HPos][it3->first - VPosStart] = it3->second / normalSum;
      }
    }
  }

  // Calculate mean and variance.
  //cout << "Mean and variance of library insert sizes." << endl;

  for (int i0 = 0; i0 < int(libVec.size()); i0++)
  {
    map<int, int>::iterator it2;
    double momZ1 = 0;
    double momC2 = 0;
    double stDev;
    int nRead = 0;

    // Calculate first zero moment.
    for (it2 = libVec[i0].insertSizeFreqMap.begin();
        it2 != libVec[i0].insertSizeFreqMap.end(); it2++)
    {
      nRead += it2->second;
      momZ1 += it2->first * it2->second;
    }

    if (nRead == 0)
    {
      libVec[i0].mean = numeric_limits < T > ::quiet_NaN();
      libVec[i0].stDev = numeric_limits < T > ::quiet_NaN();
      continue;
    }

    momZ1 /= nRead;

    // Calculate second central moment.
    for (it2 = libVec[i0].insertSizeFreqMap.begin();
        it2 != libVec[i0].insertSizeFreqMap.end(); it2++)
    {
      double diff = (it2->first - momZ1);
      double diffSq = diff * diff;
      momC2 += diffSq * it2->second;
    }

    momC2 /= (nRead - 1);
    stDev = pow(momC2, 0.5);

    libVec[i0].mean = momZ1;
    libVec[i0].stDev = stDev;
  }

  // Calculate percentiles.
  pcg.push_back(0.0);
  pcg.push_back(2.5);
  pcg.push_back(25);
  pcg.push_back(50);
  pcg.push_back(75);
  pcg.push_back(97.5);
  pcg.push_back(100.0);

  for (int i0 = 0; i0 < int(libVec.size()); i0++)
  {
    std::vector<T> pct;
    pct.resize(pcg.size());

    int result = StrIndelSize::Percentile(pct, libVec[i0].insertSizeFreqMap,
        pcg, libVec[i0].libTag);

    if (result != 0 && result != 3)
      continue;

    libVec[i0].pct0 = pct[0];
    libVec[i0].pct2p5 = pct[1];
    libVec[i0].pct25 = pct[2];
    libVec[i0].pct50 = pct[3];
    libVec[i0].pct75 = pct[4];
    libVec[i0].pct97p5 = pct[5];
    libVec[i0].pct100 = pct[6];
  }

  // Print statistics.
  /*
   cout << "Print statistics of libraries" << endl;

   for (int i0 = 0; i0 < int(libVec.size()); i0++)
   {
   cout <<
   libVec[i0].libTag << " " <<
   libVec[i0].nRead << " " <<
   libVec[i0].mean << " " <<
   libVec[i0].stDev << " " <<
   libVec[i0].pct2p5 << " " <<
   libVec[i0].pct25 << " " <<
   libVec[i0].pct50 << " " <<
   libVec[i0].pct75 << " " <<
   libVec[i0].pct97p5 << " " <<
   endl;
   }
   */

  return 0;
}

int StrIndelSize::CalcVar(
    multi_array<T, 2> &indelLoglVar,
    const std::vector<multi_array<T, 2> > &indelLogl,
    int nIter)
{
  int nRead = int(indelLogl.size());
  int nRow = indelLogl[0].shape()[0];

  for (int i0 = 0; i0 < nIter; i0++)
  {
    int readRand = rand() % nRead;

    for (int i1 = 0; i1 < nRow; i1++)
      for (int i2 = 0; i2 < i1 + 1; i2++)
        indelLoglVar[i1][i2] += indelLogl[readRand][i1][i2];
  }

  return 0;
}

int StrIndelSize::ClearInsertSizeMap()
{
  for (int i0 = 0; i0 < int(libVec.size()); i0++)
  {
    libVec[i0].insertSizeFreqMap.clear();
    libVec[i0].nRead = 0;
  }

  return 0;
}

int StrIndelSize::ClearLib()
{
  libVec.resize(1);
  ClearInsertSizeMap();

  libSet.clear();

  idLibMap.clear();
  libIdMultiMap.clear();
  idLibPosMap.clear();
  libTagLibPosMap.clear();

  GenIdLibPosMap();
  GenLibTagLibPosMap();

  CalcMpersPdf();

  return 0;
}

// Function to clear reads that do not have a library or whose library is not in use due to user
// settings (e.g., the user not wanting to use particular library).
int StrIndelSize::ClearSrp(
    std::vector<Read> &readTo,
    const std::list<Read> &readFrom)
{
  std::list<Read>::const_iterator cit0;
  map<string, int>::const_iterator cit1;
  std::vector<Read>::iterator it0;

  for (cit0 = readFrom.begin(); cit0 != readFrom.end(); cit0++)
  {
    // Find ID string.
    int libPosNo;
    int insertSize = cit0->insertSize; //xisize;
    int libInsertSizeMax;

    // Check whether ID of read belongs to an existing library.
    cit1 = idLibPosMap.find(cit0->id);

    if (cit1 == idLibPosMap.end())
    {
      cout << "No library found for read group" << cit0->id << "." << endl;
      continue;
    }

    // Check whether insert size is within range of library.
    libPosNo = cit1->second;
    libInsertSizeMax = libVec[libPosNo].libInsertSizeMax;

    if (insertSize < insertSizeMin || insertSize > libInsertSizeMax)
    {
      continue;
    }

    // Check whether library has been picked.
    if (libVec[libPosNo].isPicked == false)
    {
      continue;
    }

    // Check whether library is big enough.
    if (libVec[libPosNo].nRead < libSizeMin)
    {
      continue;
    }

    // Approve read-pair for use in calculating probabilities of indel size.
    readTo.push_back(*cit0);
    it0 = readTo.end();
    it0--;
    it0->libPosNo = libPosNo;
  }

  return 0;
}

int StrIndelSize::ClearStrLoc()
{
  int strLocSize = int(strLocDeque.size());

  if (strLocSize > 0)
  {
    strLocDeque.clear();
    cout << strLocSize << " STR locus coordinates have been deleted." << endl;
  }
  else
  {
    cout << "There are no STR locus coordinates loaded." << endl;
  }

  return 0;
}

int StrIndelSize::CloseBam()
{
  // Close BAM file.
  if (bamReader->IsOpen() == 1)
  {
    bamReader->Close();
    cout << "Bam file \"" << bam_path.string() << "\" closed successfully." << endl;
  }
  else if (bamReader->IsOpen() == 0)
  {
    cout << "No Bam file is currently open." << endl;
    return 1;
  }

  // There is no SAM file.
  //bamFileName = "";
  //bamFilePath = "";
  //bamFileLoc = "";

  // Clear all libraries.
  //ClearLib();

  return 0;
}

// Delete a library given a library tag.
int StrIndelSize::DelLib(string libOld)
{
  std::vector<StrLibrary>::iterator it0;
  bool found = false;

  for (it0 = libVec.begin(); it0 != libVec.end(); it0++)
  {
    if (it0->libTag == libOld)
    {
      libVec.erase(it0);
      found = true;

      break;
    }
  }

  if (found == false)
  {
    cout << "DelLib: Library \"" << libOld << "\" could not be found." << endl;
    return 1;
  }

  cout << "Library \"" << libOld << "\" is deleted." << endl;

  // This should be run, but I cannot find the function (deleted???).
  //StrIndelSize::GenIdLibPtMap();
  //StrIndelSize::GenLibTagLibPtMap();

  return 0;
}

int StrIndelSize::DelLib(const char *libOld, int libOldLen)
{
  int result;
  string libOld2(libOld, libOldLen);

  result = this->StrIndelSize::DelLib(libOld2);

  return result;
}

// Estimate STR size of all STR loci.
int StrIndelSize::EstStrSize(
    const std::vector<Region> &strLocVec,
    const map<string, pair<std::vector<int>, std::vector<int> > > &neighborMap)
{

  // Declare variables.
  int nStrLoc = int(strLocVec.size());

  double valMax;

  // Declare variables.
  map<string, int>::iterator it0;
  map<string*, int>::iterator it1;
  map<string*, multiset<int> >::iterator it2;
  list<Read>::iterator it3;
  list<Read>::iterator it4;

  std::map<string*, int> libCount;
  std::map<string*, multiset<int> > libFreq;
  std::vector < std::vector<Read> > strLocReadVec(nStrLoc);
  std::vector < multi_array<T, 2> > bigIndelLoglEstArr(nStrLoc);
  std::vector < multi_array<T, 2> > indelLoglArr;
  multi_array<char, 2> isVisitedArr;

  int nIndelElem = StrPd::nIndelElem;

  // Fetch addresses of library tags.
  for (int i0 = 1; i0 < int(libVec.size()); i0++)
  {
    libCount.insert(make_pair(&libVec[i0].libTag, 0));
    libFreq.insert(make_pair(&libVec[i0].libTag, multiset<int>()));
  }

  // Load reads from all STR loci.
  for (int i0 = 0; i0 < nStrLoc; i0++)
  {
    ScanStrLoc(strLocReadVec[i0], strLocVec[i0]);
  }

  // Calculate indel size estimate and variance of each locus.
  for (int i0 = 0; i0 < nStrLoc; i0++)
  {
    refBeg_now = strLocVec[i0].refBeg;

    std::list<Read> srp1;
    std::list<Read> srp1b;
    std::vector<Read> srp2;

    int nSrp;
    int strLen = strLocVec[i0].refEnd - strLocVec[i0].refBeg + 1;
    bool isNan = false;
    bool isWrongSize =
        (strLen < refSizeMin || strLen > refSizeMax) ? true : false;
    IndelPoint indelPoint = { 0, 0 };

    // Print result if STR locus is of wrong size.
    if (isWrongSize == true)
    {
      StrOutput::WriteLocusError(5, strLocVec[i0], strLen);

      continue;
    }

    // Print result if there are no reads.
    if (strLocReadVec[i0].empty() == true)
    {
      StrOutput::WriteLocusError(2, strLocVec[i0], strLen);

      continue;
    }

    // Initialize library count.
    for (it1 = libCount.begin(); it1 != libCount.end(); it1++)
      it1->second = 0;

    for (it2 = libFreq.begin(); it2 != libFreq.end(); it2++)
    {
      it2->second.clear();
    }

    // Get rid of read-pairs that are inside neighbors or that span a neighboring STR locus
    // given the current STR locus under investigation.
    nSrp = int(strLocReadVec[i0].size());

    if (checkNeighbor)
    {
      int strNearestOnLeft;
      int strNearestOnRight;

      strUtil.FindNeighborNearest(strNearestOnLeft, strNearestOnRight,
          strLocVec[i0], neighborMap);

      for (int i1 = 0; i1 < nSrp; i1++)
      {
        if (strLocReadVec[i0][i1].pos + strLocReadVec[i0][i1].length
            > strNearestOnLeft + anchorLen
            && strLocReadVec[i0][i1].pos < strNearestOnRight - anchorLen)
        {
          srp1.push_back(strLocReadVec[i0][i1]);
        }
      }
    }
    else
    {
      for (int i1 = 0; i1 < nSrp; i1++)
        srp1.push_back(strLocReadVec[i0][i1]);
    }

    // Discard of reads that have no mate.
    FindSrp(srp1b, srp1, strLocVec[i0]);

    // Clear reads that do not satisfy user preferences.

    srp2.reserve(int(srp1b.size()));
    ClearSrp(srp2, srp1b);
    nSrp = int(srp2.size());

    // Print result if no spanning read-pairs have been found.
    if (nSrp == 0)
    {
      StrOutput::WriteLocusError(3, strLocVec[i0], strLen);

      continue;
    }

    // Print result if too few spanning read-pairs have been found.
    if (nSrp < nSrpMin && nSrp > 0)
    {
      StrOutput::WriteLocusError(4, strLocVec[i0], strLen);

      continue;
    }

    // Resize and initialize the array of probabilities for each read-pair.
    indelLoglArr.resize(nSrp);

    for (int i1 = 0; i1 < nSrp; i1++)
    {
      if (int(indelLoglArr[i1].shape()[0]) != nIndelElem
          || int(indelLoglArr[i1].shape()[1]) != nIndelElem)
        indelLoglArr[i1].resize(extents[nIndelElem][nIndelElem]);
      else
      {
        for (int i2 = 0; i2 < nIndelElem; i2++)
          for (int i3 = 0; i3 < i2 + 1; i3++)
            indelLoglArr[i1][i2][i3] = 0;
      }
    }

    // Resize and initialize the array of probabilities for the posterior probabilities.
    if (int(bigIndelLoglEstArr[i0].shape()[0]) != nIndelElem)
      bigIndelLoglEstArr[i0].resize(extents[nIndelElem][nIndelElem]);

    for (int i1 = 0; i1 < nIndelElem; i1++)
      for (int i2 = 0; i2 < nIndelElem; i2++)
        bigIndelLoglEstArr[i0][i1][i2] = 0;

    // Resize and initialize the array of visit status.
    if (int(isVisitedArr.shape()[0]) != nIndelElem)
      isVisitedArr.resize(extents[nIndelElem][nIndelElem]);

    for (int i1 = 0; i1 < nIndelElem; i1++)
      for (int i2 = 0; i2 < nIndelElem; i2++)
        isVisitedArr[i1][i2] = false;

    // For each read-pair, compute conditional probabilities.
    StrCalc::EstIndelSize(indelLoglArr, isVisitedArr, libCount, libFreq, srp2,
        strLocVec[i0], libVec, idLibPosMap, StrPd::indelShift,
        StrPd::nIndelElem, StrPd::indelStepSize, insertSizeMin, refSizeMin,
        refSizeMax, anchorLen);

    // Compute the estimated posterior probability of log-likelihood estimate.
    StrCalc::EstStrSize(bigIndelLoglEstArr[i0], indelLoglArr, nSrp, nIndelElem);

    // Find greatest nonzero posterior probability.
    valMax = -1e10;
    StrCalc::FindValMax(valMax, indelPoint, bigIndelLoglEstArr[i0],
        isVisitedArr, nIndelElem);

    // Exponentiate and normalize posterior probabilities, while avoiding numerical underflow
    // and overflow.
    StrCalc::Normalize(bigIndelLoglEstArr[i0], isNan, StrPd::priorPdf,
        isVisitedArr, valMax, nIndelElem);

    // Find the greatest probability after adding prior probabilities and store the coordinates
    // in indelPoint.
    valMax = -1;
    StrCalc::FindProbMax(valMax, indelPoint, bigIndelLoglEstArr[i0],
        nIndelElem);

    // Compute the variance of the posterior probability of log-likelihood estimate.
    //vector<IndelPoint> indelPointVec;
    /*
     StrCalc::EstStrSizeCredibleInterval(
     indelInterval,
     ciProb,
     indelPoint,
     bigIndelLoglEstArr[i0],
     nIndelElem);


     cout <<
     strLocDeque[i0].refBeg << " " <<
     -StrPd::indelSizeMax + StrPd::indelStepSize * indelInterval.i1Low << "," <<
     -StrPd::indelSizeMax + StrPd::indelStepSize * indelInterval.i1 << "," <<
     -StrPd::indelSizeMax + StrPd::indelStepSize * indelInterval.i1High << " " <<
     -StrPd::indelSizeMax + StrPd::indelStepSize * indelInterval.i2Low << "," <<
     -StrPd::indelSizeMax + StrPd::indelStepSize * indelInterval.i2 << "," <<
     -StrPd::indelSizeMax + StrPd::indelStepSize * indelInterval.i2High << endl;
     */

    //indelPointVec.clear();
    /*
     StrCalc::EstStrSizePercentilePdf(
     indelPointVec,
     ciProb,
     indelPoint,
     bigIndelLoglEstArr[i0],
     indelLoglArr,
     isVisitedArr,
     StrPd::priorPd,
     nIndelElem);

     int nIndelPoint = int(indelPointVec.size());


     for (int i1 = 0; i1 < nIndelPoint; i1++)
     {
     cout <<
     (-StrPd::indelSizeMax + StrPd::indelStepSize * indelPointVec[i1].i1) << " " <<
     (-StrPd::indelSizeMax + StrPd::indelStepSize * indelPointVec[i1].i2) <<
     endl;
     }

     indelPointVec.clear();
     */

    //cout << func << 1.7b << endl;
    // Write results to file.
    StrOutput::WriteLocus(strLocVec[i0], strLen, bigIndelLoglEstArr[i0],
        StrPd::indelStepSize, isNan, isWrongSize, nSrp, logOddsMin, libCount,
        libFreq);

    // Free memory of spanning read-pairs of locus.
    strLocReadVec[i0].clear();

    // Resize array that is no longer needed.
    bigIndelLoglEstArr[i0].resize(extents[0][0]);
  }

  return 0;
}

// Function to find spanning read-pairs given a set of reads.
int StrIndelSize::FindSrp(
    std::list<Read> &readTo,
    std::list<Read> &readFrom,
    const Region strLoc)
{
  //string func = "SIS::FindSrp ";

  //set<string>        readProcessed;
  map<string, Read *> qRead;
  map<string, Read *> mRead;

  //std::vector<Read>::const_iterator  cit0;
  std::list<Read>::iterator it0;
  //std::list<Read>::iterator                 it0;
  map<string, Read *>::iterator it1;
  map<string, Read *>::iterator it2;

  //cout << func << "Starting iteration over reads." << endl;

  // Check whether reads satisfy technical criteria.
  for (it0 = readFrom.begin(); it0 != readFrom.end(); it0++)
  {
    int mapQ;
    //int onSameStrand;
    int bitFlag;
    Read *r = &*it0;

    // Check whether mapping quality is sufficient.
    mapQ = it0->mapq;

    if (mapQ < mapQMin)
    {
      continue;
    }

    // Check whether query read belongs to a 'proper pair', that it and its mate are on the
    // same strand.
    //onSameStrand = (cit0->b.core.flag & 0x2) >> 1;
    //
    //if (onSameStrand != 1)
    //    continue;

    // Check bit flags.
    bitFlag = it0->flag;
    if ((bitFlag == 83 || bitFlag == 99 || bitFlag == 147 || bitFlag == 163)
        == false)
    {
      continue;
    }

    // Check whether read has mate with same QNAME. Save it.

    if (qRead.find(it0->qname) == qRead.end())
      qRead.insert(make_pair(it0->qname, r));
    else
      mRead.insert(make_pair(it0->qname, r));
  }

  // Find read-pairs.
  for (it1 = mRead.begin(); it1 != mRead.end(); it1++)
  {
    Read *r = it1->second;

    it2 = qRead.find(it1->first);

    // Check whether read does not have a mate.
    //if (it2 == qRead.end())
    //    continue;

    // Find positions of reads of read pair.
    int x1 = it1->second->pos + 1; //int(it1->second.b.core.pos) + 1;           // Leftmost coordinate of read 1.
                                   // Adjust by one, as SAMtools uses
                                   // zero-coordinates.
    int x2 = x1 + it1->second->length - 1; //int(it1->second.b.core.l_qseq) - 1;   // Rightmost coordinate of read 1.
    int x3 = it2->second->pos + 1; // int(it2->second.b.core.pos) + 1;
    int x4 = x3 + it2->second->length - 1; //int(it2->second.b.core.l_qseq) - 1;

    //int insertSize = it1->second.b.core.isize;
    //int libInsertSizeMax = libVec[it1->second.libPosNo].libInsertSizeMax;
    int insertSize = abs(it1->second->insertSize);

    if (x3 - x1 < 0)
    {
      swap(x1, x3);
      swap(x2, x4);
    }

    //insertSize = abs(it1->second.insertSize); //(x4 - x1);

    // Check that reads are on separate sides of STR and that they have one nucleotide outside
    // of the STR.
    if (x1 > strLoc.refBeg - anchorLen || x4 < strLoc.refEnd + anchorLen)
    {
      continue;
    }

    //if ((x1 > strLoc.refBeg - anchorLen || x2 > strLoc.refEnd) ||
    //    (x3 < strLoc.refBeg || x4 < strLoc.refEnd + anchorLen))
    //{
    //    continue;
    //}

    // Check whether insert size (x4-x1, not that given by SAMtools) is within range of library.
    //if (insertSize < insertSizeMin || insertSize > libInsertSizeMax)
    //{
    //    continue;
    //}

    // Now, a read-pair is assumed to be spanning an STR locus.

    readTo.push_back(*r);
    readTo.rbegin()->insertSize = insertSize;
  }

  return 0;
}

int StrIndelSize::GenIdLibPosMap()
{
  map<string, string>::iterator it0;

  idLibPosMap.clear();

  for (it0 = idLibMap.begin(); it0 != idLibMap.end(); it0++)
  {
    string id = it0->first;
    string libTag = it0->second;

    for (int i1 = 1; i1 < int(libVec.size()); i1++)
      if (libTag == libVec[i1].libTag)
        idLibPosMap.insert(make_pair(id, i1));
  }

  return 0;
}

int StrIndelSize::GenLibTagLibPosMap()
{
  std::vector<StrLibrary>::iterator it0;

  libTagLibPosMap.clear();

  for (int i1 = 1; i1 < int(libVec.size()); i1++)
    libTagLibPosMap.insert(make_pair(libVec[i1].libTag, i1));

  return 0;
}

int StrIndelSize::GetAlignmentNext(Read &read)
{
  BamAlignment bamAl;
  static const string idTag = "RG";

  if (bamReader->GetNextAlignment(bamAl))
  {
    read.qname = bamAl.Name;
    bamAl.GetTag(idTag, read.id);
    read.refId = bamAl.RefID;
    read.length = bamAl.Length;
    read.insertSize = bamAl.InsertSize;
    read.pos = bamAl.Position;
    read.mpos = bamAl.MatePosition;
    read.mapq = bamAl.MapQuality;
    read.flag = bamAl.AlignmentFlag;
    read.cigarData = bamAl.CigarData;

    return 1;
  }

  return 0;
}

int StrIndelSize::GetAnchorLen()
{
  return anchorLen;
}

int StrIndelSize::GetBamFileLoc()
{
  cout << bam_path.string() << endl;
  return 0;
}

int StrIndelSize::GetBufferSize()
{
  return bufferSize;
}

bool StrIndelSize::GetCheckNeighbor()
{
  return checkNeighbor;
}

bool StrIndelSize::GetClipping()
{
  return clipping;
}

int StrIndelSize::GetInsertSize(int &insertSizeMin, int &insertSizeMax)
{
  insertSizeMin = this->insertSizeMin;
  insertSizeMax = this->insertSizeMax;

  return 0;
}

int StrIndelSize::GetInsertSizeMaxGrand()
{
  return insertSizeMaxGrand;
}

int StrIndelSize::GetLibInsertSize(int &libInsertSizeMax, int libNo)
{
  //insertSizeMax = this->libInsertSizeMax;

  if (libNo < 0 || libNo > int(libVec.size()) - 1)
  {
    cout << "\"libNo\" must be a number between 0 and the number of libraries ("
        << libVec.size() << ")." << endl;

    return 1;
  }

  libInsertSizeMax = libVec[libNo].libInsertSizeMax;

  return 0;
}

int StrIndelSize::GetLibInsertSizeFreq(int * insertSizeArr, int libNo)
{
  if (libNo < 0 || libNo > int(libVec.size()) - 1)
  {
    cout << "\"libNo\" must be a number between 0 and the number of libraries ("
        << libVec.size() << ")." << endl;

    return 1;
  }

  map<int, int>::iterator it0;
  int libInsertSizeMin = insertSizeMin;
  int libInsertSizeMax = insertSizeMax;
  int libInsertSizeRange = libInsertSizeMax - libInsertSizeMin + 1;

  for (int i0 = 0; i0 < libInsertSizeRange; i0++)
  {
    int rowNo = i0 * 2;
    insertSizeArr[rowNo] = i0 + libInsertSizeMin;
  }

  for (it0 = libVec[libNo].insertSizeFreqMap.begin();
      it0 != libVec[libNo].insertSizeFreqMap.end(); it0++)
  {
    int rowNo = (it0->first - libInsertSizeMin) * 2;
    insertSizeArr[rowNo + 1] = it0->second;
  }

  return 0;
}

std::vector<string> StrIndelSize::GetLibPicked()
{
  std::vector < string > libPicked;

  for (int i0 = 0; i0 < int(libVec.size()); i0++)
    if (libVec[i0].isPicked == true)
      libPicked.push_back(libVec[i0].libTag);

  return libPicked;
}

int StrIndelSize::GetLibSizeMin()
{
  return libSizeMin;
}

int StrIndelSize::GetLibTag(char *&libTag, int &libTagLen, int libNo)
{
  if (libNo < 0 || libNo > int(libVec.size()))
  {
    cout
        << "\"libTagNo\" must be a number between 0 and the number of libraries."
        << endl;
    libTag = 0;

    return 1;
  }

  //if (libTagLen != 0)
  //    delete libTag;

  libTagLen = libVec[libNo].libTag.size();
  libTag = new char[libTagLen];

  for (int i0 = 0; i0 < libTagLen; i0++)
    libTag[i0] = libVec[libNo].libTag[i0];

  return 0;
}

T StrIndelSize::GetLogOddsMin()
{
  return logOddsMin;
}

int StrIndelSize::GetMapQMin()
{
  return mapQMin;
}

int StrIndelSize::GetMapQMpersMin()
{
  return mapQMpersMin;
}

int StrIndelSize::GetMpersPdf(
    T * mpersPdf,
    const int nRowMpersPdf,
    const int nColMpersPdf,
    const int libNo)
{
  int nRow = libVec[libNo].mpersPdfArr.shape()[0];
  int nCol = libVec[libNo].mpersPdfArr.shape()[1];

  if (libNo < 0 || libNo > int(libVec.size()))
  {
    cout << "This library does not exist." << endl;
    return 1;
  }

  if (nRowMpersPdf != nRow || nColMpersPdf != nCol)
  {
    cout << "The dimensions of the MPERS PDF must be (" << nRow << ", " << nCol
        << ") and not (" << nRowMpersPdf << ", " << nColMpersPdf << ")" << endl;
    return 1;
  }

  for (int i0 = 0; i0 < nRow; i0++)
  {
    int elem_now = nCol * i0;

    for (int i1 = 0; i1 < nCol; i1++)
      mpersPdf[elem_now + i1] = libVec[libNo].mpersPdfArr[i0][i1];
  }

  return 0;
}

int StrIndelSize::GetNLib()
{
  return libVec.size();
}

int StrIndelSize::GetNSrpMin()
{
  return nSrpMin;
}

int StrIndelSize::GetNStrLoc()
{
  return int(strLocDeque.size());
}

int StrIndelSize::GetNVarIter()
{
  return nVarIter;
}

int StrIndelSize::GetReadCountMax()
{
  return readCountMax;
}

int StrIndelSize::GetRefCount()
{
  if (!this->bamReader->IsOpen())
  {
    cout << "Reference count not available." << endl;
    return -1;
  }

  int refCount = bamReader->GetReferenceCount();

  return refCount;
}

int StrIndelSize::GetRefData(string &refName, int refNo)
{
  RefVector refData;
  int refCount = this->GetRefCount();
  //bool region;
  BamAlignment bamAl;

  if (refCount < 1)
  {
    return 1;
  }

  if (refNo < 0 || refNo >= refCount)
  {
    cout << "Reference must lie in the interval [" << 0 << "," << refCount - 1
        << "]." << endl;
    return 1;
  }

  //region = bamReader->Jump(refNo);
  bamReader->GetNextAlignment(bamAl);

  if (bamAl.RefID != refNo)
  {
    refName = "";
    return 2;
  }

  refData = bamReader->GetReferenceData();
  refName = refData[refNo].RefName;

  return 0;
}

int StrIndelSize::GetRefDataFull(char *& refName, int &refNameLen, int refNo)
{
  return 0;
}

int StrIndelSize::GetRefSize(int &refSizeMin, int &refSizeMax)
{
  refSizeMin = this->refSizeMin;
  refSizeMax = this->refSizeMax;

  return 0;
}

int StrIndelSize::GetStrLoc(string &refIdStr, int &refBeg, int &refEnd, int no)
{
  if (no < 0 || no >= int(strLocDeque.size()))
  {
    cout << "Number in StrLoc must be greater than 0." << endl;
    return 1;
  }

  refIdStr = strLocDeque[no].refIdStr;
  refBeg = strLocDeque[no].refBeg;
  refEnd = strLocDeque[no].refEnd;

  return 0;
}

int StrIndelSize::InsertLib(string libNew)
{
  std::map<std::string, StrLibrary> libMap_temp;
  std::map<std::string, StrLibrary>::iterator it0;
  int count = 0;

  // Check whether library name already exists.
  if (libSet.count(libNew) != 0)
  {
    cout << "Library \"" << libNew << "\" already exists." << endl;
    return 1;
  }

  // Check whether library name is empty.
  if (libNew.size() == 0)
  {
    cout << "A new library tag must contain at least one character." << endl;
    return 2;
  }

  // Copy old libraries.
  for (int i0 = 0; i0 < int(libVec.size()); i0++)
    libMap_temp.insert(make_pair(libVec[i0].libTag, libVec[i0]));

  // Insert new library.
  libMap_temp.insert(make_pair(libNew, StrLibrary()));
  it0 = libMap_temp.find(libNew);
  it0->second.libTag = libNew;

  // Copy libraries in alphabetical order.
  libVec.clear();
  libVec.resize(libMap_temp.size());

  for (it0 = libMap_temp.begin(); it0 != libMap_temp.end(); it0++)
    libVec[count] = it0->second;

  libSet.insert(libNew);

  // Generate map values.
  StrIndelSize::GenIdLibPosMap();
  StrIndelSize::GenLibTagLibPosMap();

  // Tell what has been done.
  cout << "Library \"" << libNew << "\" has been created." << endl;

  return 0;
}

int StrIndelSize::InsertLib(const char * libNew, int libNewLen)
{
  int result;
  string libNew2(libNew, libNewLen);

  result = this->StrIndelSize::InsertLib(libNew2);

  return result;
}

int StrIndelSize::IsOpenBam()
{
  int returnVal = bamReader->IsOpen();

  if (returnVal == 1)
    cout << "Bam file \"" << bam_path.string() << "\" is open." << endl;
  else
    cout << "No Bam file is open." << endl;

  return returnVal;
}

int StrIndelSize::LaunchStrEstimation()
{
  set<string>::iterator it0;
  map<string, pair<std::vector<int>, std::vector<int> > >::iterator it1;
  std::vector<Region> strLocVec;
  map < string, pair<std::vector<int>, std::vector<int> > > neighborMap;
  set < string > refIdSet;
  timeval t1, t2;
  int nRefId = GetRefCount();
  int nStrLoc = int(strLocDeque.size());
  int isMpersPdf = false;

  // Check that at least one library has an MPERS distribution.
  for (int i0 = 0; i0 < int(libVec.size()); i0++)
  {
    if (libVec[i0].mpersPdfArr.shape()[0] > 0)
      isMpersPdf = true;
  }

  if (isMpersPdf == false)
  {
    cout << "There is no MPERS probability distribution." << endl;
    return 1;
  }

  // Check that the locus output file is open.
  if (StrOutput::IsOpenLocOutput() == false)
  {
    cout << "Locus output file is not open." << endl;
    return 2;
  }

  // Check that there are STR loci.
  if (strLocDeque.empty())
  {
    cout << "There are no STR loci stored." << endl;
    return 3;
  }

  // Check that BAM file is open.
  if (this->bamReader->IsOpen() == 0)
  {
    cout << "BAM file must be open." << endl;
    return 4;
  }

  // Check that suggested STR loci are in the BAM file.
  for (int i0 = 0; i0 < nRefId; i0++)
  {
    string refStr;
    if (GetRefData(refStr, i0) == 0)
      refIdSet.insert(refStr);
  }

  // If checking for STRs spanning neighbors of the locus under investigation. Put left and right
  // coordinates of all STR loci in separate vectors and sort them. These coordinates are used to
  // find the nearest neighboring STR locus given the locus under investigation.
  if (checkNeighbor)
  {
    strUtil.GenNeighborMap(neighborMap, strLocDeque, refIdSet);
  }

  // Start timer.
  gettimeofday(&t1, 0);

  // Compute indel sizes of STR loci.
  for (int i0 = 0; i0 < nStrLoc; i0++)
  {
    //set<string>::iterator refIdPos = refIdSet.find(strLocDeque[i0].refIdStr);
    it0 = refIdSet.find(strLocDeque[i0].refIdStr);

    if (it0 != refIdSet.end())
    {
      strLocVec.push_back(strLocDeque[i0]);

      if (bufferSize == int(strLocVec.size()))
      {
        EstStrSize(strLocVec, neighborMap);
        strLocVec.clear();
        strLocVec.reserve(bufferSize);
      }
    }
  }

  if (strLocVec.size() > 0)
  {
    EstStrSize(strLocVec, neighborMap);
  }

  // End timer.
  gettimeofday(&t2, 0);

  // Show duration.
  cout << endl << "Duration of computing indel sizes of STR loci: "
      << (t2.tv_sec - t1.tv_sec) - (t2.tv_usec - t1.tv_usec) / 1e6 << "s."
      << endl;

  return 0;
}

int StrIndelSize::LoadLibInsertSizeFreq(
    const string &libTagNew,
    std::map<int, int> &insertSizeFreqMapNew)
{
  std::map<string, StrLibrary> libOld;

  std::map<string, StrLibrary>::iterator cit0;
  std::map<int, int>::iterator cit1;

  for (int i0 = 0; i0 < int(libVec.size()); i0++)
    libOld.insert(make_pair(libVec[i0].libTag, libVec[i0]));

  cit0 = libOld.find(libTagNew);

  if (cit0 == libOld.end())
  {
    StrLibrary strLibNew;
    strLibNew.libTag = libTagNew;
    strLibNew.insertSizeFreqMap = insertSizeFreqMapNew;
    int count = 0;

    libOld.insert(make_pair(libTagNew, strLibNew));

    libVec.resize(libOld.size());
    cit0 = libOld.find("all");
    libVec[count] = cit0->second;
    //count = 1;

    for (cit0 = libOld.begin(); cit0 != libOld.end(); cit0++)
    {
      count++;

      if (cit0->second.libTag != libVec[0].libTag)
        libVec[count] = cit0->second;
    }
  }
  else
  {
    //int libPos;
    int nRead = 0;
    int nReadAll = 0;

    // Update current library.
    for (int i0 = 1; i0 < int(libVec.size()); i0++)
    {
      if (libVec[i0].libTag == libTagNew)
      {
        for (cit1 = insertSizeFreqMapNew.begin();
            cit1 != insertSizeFreqMapNew.end(); cit1++)
          nRead += cit1->second;
        libVec[i0].insertSizeFreqMap = insertSizeFreqMapNew;
        libVec[i0].nRead = nRead;
        break;
      }
    }

    // Update "all" library.
    libVec[0].insertSizeFreqMap.clear();

    for (int i0 = insertSizeMin; i0 < insertSizeMax; i0++)
    {
      nRead = 0;

      for (int i1 = 1; i1 < int(libVec.size()); i1++)
      {
        cit1 = libVec[i1].insertSizeFreqMap.find(i0);

        if (cit1 != libVec[i1].insertSizeFreqMap.end())
          nRead += cit1->second;
      }

      if (nRead > 0)
        libVec[0].insertSizeFreqMap.insert(make_pair(i0, nRead));

      nReadAll += nRead;
    }

    libVec[0].nRead = nReadAll;
  }

  return 0;
}

int StrIndelSize::LoadLibInsertSizeFreq(
    const char * libTagNew,
    int libTagNewLen,
    const int * insertSizeArr,
    const int * freqArr,
    int freqArrLen)
{
  if (libTagNewLen < 1)
  {
    cout
        << "LoadLibInsertSizeFreq: The name of the new library must contain at least one character."
        << endl;
    return 1;
  }

  if (freqArrLen < 1)
  {
    cout
        << "LoadLibInsertSizeFreq: The length of the array of insert size frequencies ("
        << freqArrLen << ") must be at least one." << endl;
    return 2;
  }

  string libTagNew2(libTagNew, libTagNewLen);
  std::map<int, int> insertSizeMapNew;
  int returnVal;
  int libInsertSizeMin = insertSizeMin;

  for (int i0 = 0; i0 < freqArrLen; i0++)
  {
    int insertSize = insertSizeArr[i0];
    int freq = freqArr[i0];

    if (insertSize >= libInsertSizeMin && freq > 0)
      insertSizeMapNew.insert(make_pair(insertSize, freq));
  }

  returnVal = StrIndelSize::LoadLibInsertSizeFreq(libTagNew2, insertSizeMapNew);

  StrIndelSize::GenIdLibPosMap();
  StrIndelSize::GenLibTagLibPosMap();

  return returnVal;
}

// Store data on STR locus (chromosome, interval).
int StrIndelSize::LoadStrLoc(string refIdStr, int refBeg, int refEnd)
{
  //int strLen = refEnd - refBeg + 1;

  if (refBeg <= 0 || refEnd <= 0)
  {
    cout << "For an STR region, refBeg > 0 and refEnd > 0." << endl;
    return 1;
  }

  /*
   if (strLen < refSizeMin || strLen > refSizeMax)
   {
   cout << "For an STR locus, the size must be greater than the minimum reference size " <<
   "(" << refSizeMin << " nt)" <<
   "and smaller than the maximum reference size " <<
   "(" << refSizeMax << " nt)" <<
   "." << endl;
   return 2;
   }
   */

  Region strLoc = { 0, refIdStr, refBeg, refEnd };

  strLocDeque.push_back(strLoc);

  return 0;
}

int StrIndelSize::LoadStrLoc(
    const char * refIdStr,
    int refIdStrLen,
    int refBeg,
    int refEnd)
{
  string refIdStr2(refIdStr, refIdStrLen);
  int returnVal = LoadStrLoc(refIdStr2, refBeg, refEnd);

  return returnVal;
}

int StrIndelSize::OpenBam(const filesystem::path &bam_fp)
{
  string fLoc = bam_fp.string(); // filePath + fileName;
  int locateIndex;

  // Close any previous BAM file that was open.

  if (bamReader->IsOpen())
  {
    bamReader->Close();
    cout << "Closing old BAM file " << bam_path.string() << endl;
  }

  // Opening new BAM file.
  bamReader.reset(new BamReader());

  if (bamReader->Open(fLoc))
  {
    cout << "Has opened BAM file \"" << fLoc << "\"." << endl;
  }
  else
  {
    cout << "Failed to open BAM file \"" << fLoc << "\"" << endl;
    return 1;
  }

  locateIndex = bamReader->LocateIndex();

  if (locateIndex == false)
  {
    cout << "Index file not located." << endl;
    cout << "Closing Bam file." << endl;
    this->CloseBam();
    return 2;
  }
  else
  {
    cout << "Index file located." << endl;
  }

  // Save filename and its path.
  //this->bamFilePath = filePath;
  //this->bamFileName = fileName;
  //this->bamFileLoc = fLoc;
  bam_path = bam_fp;

  return 0;
}

int StrIndelSize::OpenLocOutput(const filesystem::path &loc_fp)
{
  int returnVal = StrOutput::OpenLocOutput(loc_fp, libVec, mapQMin);

  return returnVal;
}

int StrIndelSize::Percentile(
    std::vector<T> &percentile,
    std::map<int, int> &samp,
    std::vector<T> &percentage,
    string libTag)
{
  std::map<int, int>::iterator it1;
  map<T, int> pFreqMap;

  int nPerc = percentage.size();
  uint64_t nObs = 0;
  T tmp_1;

  // Check that vector of percentages is not empty.
  if (nPerc == 0)
  {
    cout
        << "There are no percentiles in the object \"percentage\" of library \""
        << libTag << "\"." << endl;
    return 1;
  }

  // Check that values of vector of percentages are not empty.
  for (int i0 = 0; i0 < nPerc; i0++)
  {
    if (percentage[i0] < 0 || percentage[i0] > 100)
    {
      cout << "Any percentile must be between 0 and 100." << endl;
      return 2;
    }
  }

  // Check that there are observations.
  if (samp.size() == 0)
  {
    //cout << "Percentile: There are no observations in the library \"" << libTag << "\"." << endl;

    for (int i0 = 0; i0 < int(percentage.size()); i0++)
      percentile[i0] = numeric_limits < T > ::quiet_NaN();

    return 3;
  }

  // Count no. of observations.
  for (it1 = samp.begin(); it1 != samp.end(); it1++)
  {
    nObs += it1->second;
  }

  // Precalculate some values.
  tmp_1 = 100 / T(nObs);
  nObs = 0;

  for (it1 = samp.begin(); it1 != samp.end(); it1++)
  {
    int freq = it1->second;
    int no = it1->first;
    T p_n;

    nObs += freq;
    p_n = tmp_1 * (nObs - 0.5);

    pFreqMap.insert(make_pair(p_n, no));
  }

  // Calculate percentiles.
  //percentile.resize(nPerc);

  for (int i0 = 0; i0 < nPerc; i0++)
  {
    std::map<T, int>::iterator itLow; //, itUp;
    T p = percentage[i0];

    itLow = pFreqMap.lower_bound(p);
    //itUp = pFreqMap.upper_bound (p);

    /*
     if (itLow == pFreqMap.begin())
     {
     percentile[i0] = itLow->second;
     }
     else if (itLow == pFreqMap.end())
     {
     itLow--;
     percentile[i0] = itLow->second;
     }
     else
     {
     percentile[i0] = itLow->second;
     }
     */

    if (itLow == pFreqMap.end())
      itLow--;

    percentile[i0] = itLow->second;
  }

  return 0;
}

int StrIndelSize::PickLib(string pickedLibNew)
{
  bool isLibPicked = false;

  for (int i0 = 0; i0 < int(libVec.size()); i0++)
  {
    if (pickedLibNew == libVec[i0].libTag)
    {
      libVec[i0].isPicked = true;
      isLibPicked = true;
      break;
    }
  }

  if (isLibPicked == false)
  {
    cout << "PickLib: Library " << pickedLibNew << " could not be found."
        << endl;
    return 1;
  }

  cout << "Library \"" << pickedLibNew << "\" is now picked." << endl;

  return 0;
}

int StrIndelSize::PickLib(const char * pickedLibNew, int pickedLibNewLen)
{
  string pickedLibNew2(pickedLibNew, pickedLibNewLen);
  int returnVal = PickLib(pickedLibNew2);

  return returnVal;
}

int StrIndelSize::PrintLibInfo(int type)
{
  int nLibTagLenMax = 0;

  for (int i0 = 0; i0 < int(libVec.size()); i0++)
  {
    nLibTagLenMax = max(nLibTagLenMax, int(libVec[i0].libTag.size()));
  }

  nLibTagLenMax = max(nLibTagLenMax, 7) + 2;

  cout << "Library" << setw(nLibTagLenMax - 7) << "";
  cout << "no. of reads" << endl;

  for (int i0 = 0; i0 < int(libVec.size()); i0++)
  {
    int libTagLen = libVec[i0].libTag.size();
    cout << libVec[i0].libTag << setw(nLibTagLenMax - libTagLen) << "";
    cout << libVec[i0].nRead << endl;
  }

  return 0;
}

int StrIndelSize::ReadBamHeader()
{
  if (bamReader->IsOpen() == false)
  {
    cout << "BAM file must be open for reading header." << endl;

    return 1;
  }

  // Define variables
  set<string>::iterator it1;
  std::map<int, int> allLib;
  std::vector<StrLibrary> libOldVec(libVec);
  int libCount = 0;

  SamReadGroupDictionary rgDict;
  SamReadGroupIterator rgIt0;

  SamHeader samHeader = bamReader->GetHeader();
  //string header = bamReader->GetHeaderText();
  rgDict = samHeader.ReadGroups;

  // Find unique read group identifier and their associated libraries.
  for (rgIt0 = rgDict.Begin(); rgIt0 != rgDict.End(); rgIt0++)
  {
    string id = rgIt0->ID;
    string lb = rgIt0->Library;

    // Insert elements.
    if (id.length() > 0 && lb.length() > 0)
    {
      idLibMap.insert(make_pair(id, lb));
      libIdMultiMap.insert(make_pair(lb, id));
      libSet.insert(lb);
      libCount++;
    }
  }

  // Return if no new library could be found.
  if (libCount == 0)
  {
    cout << "No new library could be found." << endl;
    return 2;
  }

  // Create objects for libraries.
  it1 = libSet.begin();
  libVec.resize(libSet.size() + 1);

  libVec[0].libTag = "all";

  for (int i0 = 1; i0 < int(libVec.size()); i0++)
  {
    bool libFound = false;

    // Existing library has been found.
    for (int i1 = 1; i1 < int(libOldVec.size()); i1++)
    {
      if (libOldVec[i1].libTag == *it1)
      {
        libVec[i0] = libOldVec[i1];
        libFound = true;
        break;
      }
    }

    // Add new library.
    if (libFound == false)
      libVec[i0].libTag = *it1;

    // Increase iterator.
    it1++;
  }

  // Calculate the number of reads in all libraries.
  //StrIndelSize::CalcLibInsertSize();

  // Generate maps IDs and library tags to library memory addresses.
  StrIndelSize::GenIdLibPosMap();
  StrIndelSize::GenLibTagLibPosMap();

  return 0;
}

//int StrIndelSize::ReadInsertSize(const string &fPath, const string &fName)
int StrIndelSize::ReadInsertSize(const filesystem::path &f_path)
{
  int insertSizeMinFile;
  int insertSizeMaxFile;
  int mapQMinFile;

  int returnVal = StrInput::ReadInsertSize(libVec, insertSizeMinFile,
      insertSizeMaxFile, mapQMinFile, f_path);
  if (returnVal)
    return returnVal;

  returnVal = SetInsertSize(insertSizeMinFile, insertSizeMaxFile);
  if (returnVal)
    return returnVal;

  returnVal = SetLibInsertSize(insertSizeMaxFile, "general");
  if (returnVal)
    return returnVal;

  returnVal = SetMapQMin(mapQMin);
  if (returnVal)
    return returnVal;

  // Calculate the number of reads.
  StrIndelSize::CalcLibInsertSize();

  // Update ?
  UpdateInsertSizeMaxGrand();

  return returnVal;
}

int StrIndelSize::ReadInsertSize(
    const string &refIdStr,
    const int refBeg = 0,
    const int refEnd = 0)
{
  if (bamReader->IsOpen() == 0)
  {
    cout << "BAM file must be open." << endl;
    return 1;
  }

  if (libVec.size() == 1)
  {
    cout << "There are no libraries." << endl;
    return 2;
  }

  if (readCountMax <= libVec[0].nRead)
  {
    cout << "The maximum allowed no. of read-pairs to sample is "
        << readCountMax << ". " << libVec[0].nRead
        << " reads have already been sampled." << endl;
    return 3;
  }

  int count = libVec[0].nRead;
  timeval t1, t2;

  map<int, int>::iterator it0;
  map<string, int>::iterator it1;

  BamAlignment bamAl;
  Region reg = { 0, refIdStr, refBeg, refEnd };
  bool region = SetRegion(reg);
  string tag = "RG";

  if (region == false)
  {
    return 4;
  }

  // Start timer.
  gettimeofday(&t1, 0);

  // Obtain reads.
  while (bamReader->GetNextAlignment(bamAl) && count < readCountMax)
  {
    int insertSize_tmp = bamAl.InsertSize;
    int mapQ_tmp = bamAl.MapQuality;
    int insertSizeAbs_tmp;
    string id_tmp;
    string tmp_1;
    int libPos;

    // If insert size is too small or too big.
    if (insertSize_tmp < insertSizeMin || insertSize_tmp > insertSizeMax)
      continue;

    // If MAPQ score for MPERS is too small.
    if (mapQ_tmp < mapQMpersMin)
      continue;

    // Ensure that insert size is positive
    insertSizeAbs_tmp = abs(insertSize_tmp);

    if (bamAl.GetTag(tag, id_tmp) == false)
      continue;

    // Put insert size of read into correct library.
    it1 = idLibPosMap.find(id_tmp);

    if (it1 == idLibPosMap.end())
      continue;

    libPos = it1->second;

    it0 = libVec[libPos].insertSizeFreqMap.find(insertSizeAbs_tmp);

    if (it0 == libVec[libPos].insertSizeFreqMap.end())
      libVec[libPos].insertSizeFreqMap[insertSizeAbs_tmp] = 1;  // New object
    else
    {
      it0->second += 1;                       // Existing object
    }

    count++;
  }

  // End timer.
  gettimeofday(&t2, 0);

  // Show duration.
  cout << "\nDuration of counting insert size frequencies: "
      << CalcDuration(t1, t2) << " sec." << endl;

  // Calculate the number of reads.
  StrIndelSize::CalcLibInsertSize();

  return 0;
}

int StrIndelSize::ReadInsertSize(
    const char * refIdStr,
    int refIdStrLen,
    int refBeg,
    int refEnd)
{
  int returnVal;
  string refIdStr2(refIdStr, refIdStrLen);

  returnVal = StrIndelSize::ReadInsertSize(refIdStr2, refBeg, refEnd);

  return returnVal;
}

//int StrIndelSize::ReadLocus(const string &fPath, const string &fName)
int StrIndelSize::ReadLocus(const filesystem::path &f_path)
{
  int returnVal = StrInput::ReadLocus(strLocDeque, f_path);

  return returnVal;
}

int StrIndelSize::ReadPrior(
    const filesystem::path &f_path,
    int nLineSkip,
    int indelSizeMaxNew,
    int indelStepSizeNew)
{
  list<T> priorList;

  int returnVal = StrInput::ReadPriorPdf(priorList, f_path, nLineSkip,
      indelSizeMaxNew, indelStepSizeNew);

  if (returnVal)
    return returnVal;

  returnVal = StrPd::SetPrior(priorList, indelSizeMaxNew, indelStepSizeNew);

  return returnVal;
}

// Find reads in all STR loci.
int StrIndelSize::ReadReadLen(
    const filesystem::path &bam_fp,
    const Region &region,
    int readLenMax,
    int readLenMapQMin,
    int nReadMax)
{
  if (readLenMax < 1)
  {
    cout << "Read length must lie in the interval [1,1000]." << endl;
    return 1;
  }

  if (readLenMapQMin < 0 || readLenMapQMin > 100)
  {
    cout << "MAPQ minimum must lie in the interval [0,100]." << endl;
    return 2;
  }

  if (nReadMax < 1)
  {
    cout << "The maximum number of reads must be equal to or greater than one."
        << endl;
    return 3;
  }

  bool regionValid = SetRegion(region);

  if (!regionValid)
  {
    return 4;
  }

  // Open output file.
  if (StrOutput::OpenReadLenOutput(bam_fp))
  {
    return 5;
  }

  int nLib = int(libVec.size());
  multi_array<int, 2> readLenArr;
  map<string, int>::iterator it0;
  int (StrUtil::*lenFunc)(const std::vector<CigarOp> &cigarData, int readLen);

  BamAlignment bamAl;
  const string tag = "RG";
  int count = 0;
  timeval t1, t2;
  lenFunc =
      clipping ?
          &StrUtil::ReadLenWithoutClipping : &StrUtil::ReadLenWithClipping;
  //readLenArr.resize(extents[readLenMax+1][nLib]);

  // Start timer.
  gettimeofday(&t1, 0);

  // Obtain reads.
  while (bamReader->GetNextAlignment(bamAl) && count < nReadMax)
  {
    int readLenUnclipped = bamAl.Length;    // Read length.
    int readLen_tmp = (strUtil.*lenFunc)(bamAl.CigarData, bamAl.Length); // Remove clipping if applicable.

    string id_tmp;
    int libPos;
    pair<int, int> pair_tmp;

    // If insert size is too small or too big.
    if (readLen_tmp > readLenMax)
      continue;

    // Get tag of read for finding correct library.
    if (bamAl.GetTag(tag, id_tmp) == false)
      continue;

    // Find correct library.
    it0 = idLibPosMap.find(id_tmp);

    if (it0 == idLibPosMap.end())
      continue;

    libPos = it0->second;

    // Increase counter.
    pair_tmp.first = readLenUnclipped;
    pair_tmp.second = readLen_tmp;

    //readLenArr[readLen_tmp][libPos]++;
    libVec[libPos].overlapFfMap[pair_tmp]++;

    count++;
  }

  // End timer.
  gettimeofday(&t2, 0);

  // Show duration.
  cout << "\nDuration of counting read length frequencies: "
      << CalcDuration(t1, t2) << " sec." << endl;

  // Compute sum of all libraries for each read length.
  /*
   for(int i0 = 0; i0 < readLenMax+1; i0++)
   {
   int sum = 0;

   for(int i1 = 1; i1 < nLib; i1++)
   {
   sum += readLenArr[i0][i1];
   }

   readLenArr[i0][0] = sum;
   }
   */

  for (int i0 = 0; i0 < nLib; i0++)
  {
    map<pair<int, int>, int>::iterator it1;
    //int rLenPrev = -1;

    for (it1 = libVec[i0].overlapFfMap.begin();
        it1 != libVec[i0].overlapFfMap.end(); it1++)
    {
      int rLen = it1->first.first;
      int freq = it1->second;
      libVec[i0].rLenNRead[rLen] += freq;
    }
  }

  // Print result to file.
  StrOutput::WriteReadLen(readLenArr, libVec, readLenMax, readLenMapQMin);

  // Close output file.
  StrOutput::CloseReadLenOutput();

  return 0;
}

// Find all reads in a region (e.g., flanking regions of an STR locus) and a stores them in a list
// ("readVec").
int StrIndelSize::ScanRegion(std::list<Read> &readVec, const Region strLoc)
{
  Read read;
  string refIdStr = strLoc.refIdStr;
  int refId;
  int refBeg = strLoc.refBeg;
  int refEnd = strLoc.refEnd;
  int result;

  string bamRegion = strLoc.refIdStr + ":" + lexical_cast<string>(refBeg) + "-"
      + lexical_cast<string>(refEnd);

  // Parse region.
  refId = bamReader->GetReferenceID(refIdStr);
  result = bamReader->SetRegion(refId, refBeg, refId, refEnd);

  if (result == false)
  {
    cout << "Invalid region in BAM file: " << bamRegion << "." << endl;
    return 1;
  }

  // Get reads.
  while (GetAlignmentNext(read))
  {
    readVec.push_back(read);
  }

  return 0;
}

// Find all reads before, inside and after a given STR locus.
int StrIndelSize::ScanStrLoc(std::vector<Read> &readVec, const Region strLoc)
{
  string refIdStr = strLoc.refIdStr;
  int refId;
  int refBeg = strLoc.refBeg - insertSizeMaxGrand - 1;
  int refEnd = strLoc.refEnd + insertSizeMaxGrand + 1;
  int result;

  list<Read> readList_tmp;
  list<Read>::iterator it1;
  string bamRegion = strLoc.refIdStr + ":" + lexical_cast<string>(refBeg) + "-"
      + lexical_cast<string>(refEnd);
  string idTag = "RG";
  BamAlignment bamAl;

  // Parse region.
  refId = bamReader->GetReferenceID(refIdStr);
  result = bamReader->SetRegion(refId, refBeg, refId, refEnd);

  if (result == false)
  {
    cout << "Invalid region in BAM file: " << bamRegion << "." << endl;
    return 1;
  }

  // Get reads.
  while (bamReader->GetNextAlignment(bamAl))
  {
    Read read;

    read.qname = bamAl.Name;
    bamAl.GetTag(idTag, read.id);
    read.length = bamAl.Length;
    read.insertSize = bamAl.InsertSize;
    read.pos = bamAl.Position;
    read.mpos = bamAl.MatePosition;
    read.mapq = bamAl.MapQuality;
    read.flag = bamAl.AlignmentFlag;
    readList_tmp.push_back(read);
  }

  // Save reads.
  //readVec.resize(readList_tmp.size());
  readVec.reserve(readList_tmp.size());
  it1 = readList_tmp.begin();

  for (int i0 = 0; i0 < int(readList_tmp.size()); i0++)
  {
    //readVec[i0] = *it1;
    readVec.push_back(*it1);
    it1++;
  }

  return 0;
}

int StrIndelSize::SetAnchorLen(int anchorLenNew)
{
  if (anchorLenNew < 1)
  {
    cout << "Anchor length must be at least 1." << endl;
    return 1;
  }

  this->anchorLen = anchorLenNew;
  strUtil.SetAnchorLen(anchorLenNew);

  return 0;
}

int StrIndelSize::SetBufferSize(int bufferSizeNew)
{
  if (bufferSize < 0 || bufferSize > 1000)
  {
    cout << "Buffer size is limited to the range: 0 < x < 1000." << endl;
    return 1;
  }

  bufferSize = bufferSizeNew;

  return 0;
}

int StrIndelSize::SetCheckNeighbor(bool checkNeighborNew)
{
  checkNeighbor = checkNeighborNew;
  return 0;
}

int StrIndelSize::SetClipping(bool clippingNew)
{
  clipping = clippingNew;
  return 0;
}

int StrIndelSize::SetInsertSize(int insertSizeMinNew, int insertSizeMaxNew)
{
  if (insertSizeMinNew < 1 || insertSizeMaxNew < 1)
  {
    cout << "Insert size minimum and maximum must be larger than 1 nucleotide."
        << endl;
    return 1;
  }

  if (insertSizeMinNew > insertSizeMaxNew)
  {
    cout << "insertSizeMin <= insertSizeMax." << endl;
    return 2;
  }

  this->insertSizeMin = insertSizeMinNew;
  this->insertSizeMax = insertSizeMaxNew;

  StrIndelSize::CalcMpersPdf();

  return 0;
}

int StrIndelSize::SetLibInsertSize(
    int libInsertSizeMaxNew,
    const string &libTag)
{
  map<string, int>::iterator it0 = libTagLibPosMap.find(libTag);
  int pos;

  if (libInsertSizeMaxNew < insertSizeMin)
  {
    cout
        << "The insert size of a library must be greater than or equal to the insert size minimum "
        << insertSizeMin << "." << endl;

    return 3;
  }

  // Update all libraries.
  if (libTag == "general")
  {
    for (int i0 = 0; i0 < int(libVec.size()); i0++)
      libVec[i0].libInsertSizeMax = libInsertSizeMaxNew;

    // Find new grand maximum.
    UpdateInsertSizeMaxGrand();

    return 0;
  }

  // Update a particular library.
  if (it0 == libTagLibPosMap.end())
  {
    cout << "SetLibInsertSize: The library \"" << libTag
        << "\" could not be found" << endl;
    return 2;
  }

  // Change the insert size maximum of a library.
  pos = it0->second;
  libVec[pos].libInsertSizeMax = libInsertSizeMaxNew;

  // Find new grand maximum.
  UpdateInsertSizeMaxGrand();

  // Recalculate the MPERS distributions.
  StrIndelSize::CalcMpersPdf();

  return 0;
}

int StrIndelSize::SetLibInsertSize(
    int libInsertSizeMaxNew,
    const char * libTag,
    int libTagLen)
{
  string libTag2(libTag, libTagLen);

  if (libTagLen < 1)
  {
    cout << "The length of the library tag must be at least 1." << endl;
    return 1;
  }

  int returnVal = StrIndelSize::SetLibInsertSize(libInsertSizeMaxNew, libTag2);

  return returnVal;
}

int StrIndelSize::SetLibInsertSizeFreq(
    const int * insertSizeArrNew,
    const int arrLen,
    const int libNo)
{
  int nRead = 0;

  if (libNo < 1 || libNo > int(libVec.size()) - 1)
  {
    cout << "\"libNo\" must be a number between 1 and the number of libraries ("
        << libVec.size() << ")." << endl;

    return 1;
  }

  libVec[libNo].insertSizeFreqMap.clear();

  for (int i0 = 0; i0 < arrLen; i0++)
  {
    int insertSize = insertSizeArrNew[2 * i0];
    int freq = insertSizeArrNew[2 * i0 + 1];

    if (insertSize < insertSizeMin || insertSize > insertSizeMax || freq < 1)
      continue;

    libVec[libNo].insertSizeFreqMap.insert(make_pair(insertSize, freq));
    nRead += freq;
  }

  libVec[libNo].nRead = nRead;

  this->CalcLibInsertSize();

  return 0;
}

int StrIndelSize::SetLibSizeMin(int libSizeMinNew)
{
  if (libSizeMinNew < 1)
  {
    cout
        << "The minimum size of a library for use in inferring the indel size must "
        << "be at least 1.";
    return 1;
  }

  this->libSizeMin = libSizeMinNew;

  return 0;
}

int StrIndelSize::SetLogOddsMin(T logOddsMinNew)
{
  if (logOddsMinNew < 0)
  {
    cout << "The minimum of the log-odds ratio is 0." << endl;
    return 1;
  }

  this->logOddsMin = logOddsMinNew;

  return 0;
}

int StrIndelSize::SetMapQMin(int mapQMinNew)
{
  if (mapQMinNew < 0 || mapQMinNew > 100)
  {
    cout << "0 <= mapQMin <= 100." << endl;
    return 1;
  }

  this->mapQMin = mapQMinNew;

  return 0;
}

int StrIndelSize::SetMapQMpersMin(int mapQMpersMinNew)
{
  if (mapQMpersMinNew < 0 || mapQMpersMinNew > 100)
  {
    cout << "0 <= mapQMpers <= 100." << endl;
    return 1;
  }

  this->mapQMpersMin = mapQMpersMinNew;

  return 0;
}

int StrIndelSize::SetMpersPdf(
    const char * libTag,
    int libTagLen,
    const T * mpersPdfNew,
    int nRow,
    int nCol)
{
  string libTagStr(libTag, libTagLen);
  bool libFound = false;
  int nRowMpersPdf = libVec[0].mpersPdfArr.shape()[0];
  int nColMpersPdf = libVec[0].mpersPdfArr.shape()[1];

  if (nRow != nRowMpersPdf || nCol != nColMpersPdf)
  {
    cout << "nRow and nCol must equal the dimensions of existing MPERS arrays ";
    cout << "(" << nRowMpersPdf << ", " << nColMpersPdf << ")." << endl;
    return 1;
  }

  if (libTagLen < 1)
  {
    cout << "The length of the library tag must be one character or longer."
        << endl;
    return 2;
  }

  for (int i0 = 0; i0 < int(libVec.size()); i0++)
  {
    if (libVec[i0].libTag == libTagStr)
    {

      for (int i1 = 0; i1 < nRowMpersPdf; i1++)
      {
        int elemNow = nColMpersPdf * i1;

        for (int i2 = 0; i2 < nColMpersPdf; i2++)
          libVec[i0].mpersPdfArr[i0][i1] = mpersPdfNew[elemNow + i2];
      }

      libFound = true;
      break;
    }
  }

  if (libFound == false)
  {
    cout << "SetMpers: The library \"" << libTag << "\" could not be found"
        << endl;
    return 3;
  }

  return 0;
}

int StrIndelSize::SetNSrpMin(int nSrpMinNew)
{
  if (nSrpMinNew < 1)
  {
    cout
        << "The number of spanning reads used to calculate the inferred indel size ";
    cout << "must be at least 1." << endl;
    return 1;
  }

  this->nSrpMin = nSrpMinNew;

  return 0;
}

int StrIndelSize::SetNVarIter(int nVarIterNew)
{
  if (nVarIterNew < 0)
  {
    cout << "nVarIter >= 0." << endl;
    return 1;
  }

  nVarIter = nVarIterNew;

  return 0;
}

int StrIndelSize::SetReadCountMax(int readCountMaxNew)
{
  if (readCountMaxNew < 1)
  {
    cout << "readCountMax > 0." << endl;
    return 1;
  }

  this->readCountMax = readCountMaxNew;

  return 0;
}

int StrIndelSize::SetRegion(Region region)
{
  bool regionValid;

  if (region.refIdStr.size() == 0)
  {
    regionValid = bamReader->Rewind();
  }
  else
  {
    region.refId = bamReader->GetReferenceID(region.refIdStr);

    if (region.refBeg > 0 && region.refEnd == 0)
    {
      regionValid = bamReader->Jump(region.refId, region.refBeg);
    }
    else if (region.refBeg > 0 && region.refEnd > 0
        && region.refBeg < region.refEnd)
    {
      regionValid = bamReader->SetRegion(region.refId, region.refBeg,
          region.refId, region.refEnd);
    }
    else
    {
      regionValid = bamReader->Jump(region.refId, 0);
    }
  }

  if (regionValid == false)
  {
    cout << "Could not find position of region " << region.refIdStr << ":"
        << region.refBeg << "-" << region.refEnd << "." << endl;
  }

  return regionValid;
}

int StrIndelSize::SetRefSize(int refSizeMinNew, int refSizeMaxNew)
{
  if (refSizeMinNew < 0 || refSizeMaxNew < 1)
  {
    cout << "Reference microsatellites must be one nucleotide or larger."
        << endl;
    return 1;
  }

  if (refSizeMaxNew < refSizeMinNew)
  {
    cout << "refSizeMax >= refSizeMin." << endl;
    return 2;
  }

  this->refSizeMin = refSizeMinNew;
  this->refSizeMax = refSizeMaxNew;

  StrIndelSize::CalcMpersPdf();

  return 0;
}

int StrIndelSize::UnpickLib(const string &pickedLibNew)
{
  bool isLibPicked = true;

  for (int i0 = 0; i0 < int(libVec.size()); i0++)
  {
    if (pickedLibNew == libVec[i0].libTag)
    {
      libVec[i0].isPicked = false;
      isLibPicked = false;
      break;
    }
  }

  if (isLibPicked == true)
  {
    cout << "UnpickLib: Library \"" << pickedLibNew
        << "\" could not be found for unpicking." << endl;
    return 1;
  }

  cout << "Library \"" << pickedLibNew << "\" is now unpicked." << endl;

  return 0;
}

int StrIndelSize::UnpickLib(const char * unpickedLibNew, int unpickedLibNewLen)
{
  string unpickedLibNew2(unpickedLibNew, unpickedLibNewLen);
  int returnVal = UnpickLib(unpickedLibNew2);

  return returnVal;
}

int StrIndelSize::UpdateInsertSizeMaxGrand()
{
  insertSizeMaxGrand = 100;

  for (int i0 = 1; i0 < int(libVec.size()); i0++)
  {
    if (libVec[i0].libInsertSizeMax > insertSizeMaxGrand)
    {
      insertSizeMaxGrand = libVec[i0].libInsertSizeMax;
      libVec[0].libInsertSizeMax = insertSizeMaxGrand;
    }
  }

  return 0;
}

int StrIndelSize::WriteInsertSize(const filesystem::path &f_path)
{
  int returnVal = StrOutput::OpenInsertSizeOutput(f_path);

  if (returnVal != 0)
  {
    return returnVal;
  }

  this->StrOutput::WriteInsertSize(libVec, mapQMpersMin, insertSizeMin,
      insertSizeMax);

  returnVal = StrOutput::CloseInsertSizeOutput();

  return 0;
}

int StrIndelSize::WriteMpersPdf(const filesystem::path &f_path)
{
  int returnVal = StrOutput::OpenMpersPdfOutput(f_path);

  if (returnVal != 0)
    return returnVal;

  for (int i0 = 0; i0 < int(libVec.size()); i0++)
  {
    this->StrOutput::WriteMpersPdf(libVec[i0].mpersPdfArr, libVec[i0].libTag);
  }

  returnVal = StrOutput::CloseMpersPdfOutput();

  return 0;
}

int StrIndelSize::WritePriorPdf(const filesystem::path &f_path)
{
  int returnVal = StrOutput::OpenPriorPdfOutput(f_path);

  if (returnVal != 0)
    return returnVal;

  this->StrOutput::WritePriorPdf(StrPd::priorPdf);

  returnVal = StrOutput::ClosePriorPdfOutput();

  return 0;
}

#endif /* STRINDELSIZE_HPP_ */
