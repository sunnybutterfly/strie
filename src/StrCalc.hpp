#ifndef STR_CALC_H
#define STR_CALC_H

// STL libraries
#include <algorithm>    // sort()
#include <cmath>        // log(), exp()
#include <cstdlib>      // abs(), rand
#include <fstream>      // open(), close()
#include <iomanip>      // Insertion and extraction operators.
#include <iostream>
#include <map>
#include <limits>
#include <set>
#include <string>
#include <utility>      // pair
#include <vector>

// External libraries
#include <boost/multi_array.hpp>

// Own files
#include "type.hpp"
#include "StrLibrary.hpp"

// Namespace
using namespace std;
using namespace boost;

// Class
class StrCalc
{
public:
  StrCalc();
  ~StrCalc();
  StrCalc(const StrCalc &rhs);

  int CalcCondProb(
      multi_array<T, 2> & indelLoglArr,
      multi_array<char, 2> & isVisitedArr,
      const Read & read,
      const multi_array<T, 2> & mpersPdfArr,
      const vector<int> & indelShift,
      const int center,
      const int strLen,
      const int nIndelStepLeft,
      const int nIndelStepRight,
      const int nIndelStepDown,
      const int nIndelStepUp,
      const int libInsertSizeMin,
      const int libInsertSizeMax,
      const int refSizeMin,
      const int refSizeMax,
      const int anchorLen,
      const int indelStepSize,
      const int nIndelElem,
      const string & libName,
      const Region & strLoc);

  int EstIndelSize(
      vector<multi_array<T, 2> > & indelLoglArr,
      multi_array<char, 2> & isVisitedArr,
      map<string*, int> & libCount,
      map<string*, multiset<int> > & libFreq,
      const vector<Read> & srpVec,
      const Region strLoc,
      vector<StrLibrary> & libVec,
      const map<string, int> & idLibPosMap,
      const vector<int> & indelShift,
      int nIndelElem,
      int indelStepSize,
      int insertSizeMin,
      int refSizeMin,
      int refSizeMax,
      int anchorLen);

  int EstStrSize(
      multi_array<T, 2> & indelLoglEstArr,
      const std::vector<multi_array<T, 2> > & indelLoglArr,
      //const multi_array <T, 2>                & priorPd,
      int nSrp,
      int nIndelElem);

  int EstStrSizeBootstrap1(
      vector<IndelPoint> & indelPointVec,
      const T ciProb,
      const IndelPoint indelPoint,
      const multi_array<T, 2> & indelLoglArr,
      const int nIndelElem);

  int EstStrSizeCredibleInterval(
      IndelInterval & indelInterval,
      const T ciProb,
      const IndelPoint indelPoint,
      const multi_array<T, 2> & indelLoglArr,
      const int nIndelElem);

  int EstStrSizePercentilePd(
      vector<IndelPoint> & indelPointVec,
      const T ciProb,
      const IndelPoint indelPoint,
      const multi_array<T, 2> & indelLoglArr,
      const std::vector<multi_array<T, 2> > & indelLoglArrVec,
      const multi_array<char, 2> & isVisistedArr,
      const multi_array<T, 2> & priorPd,
      const int nIndelElem);

  int EstStrSizePercentilePdError(
      IndelInterval & indelInterval,
      const T ciProb,
      const IndelPoint indelPoint,
      const multi_array<T, 2> & indelLoglArr,
      const std::vector<multi_array<T, 2> > & indelLoglArrVec,
      const int nIndelElem);

  int EstStrSizePercentilePdBiasCor(
      IndelInterval & indelInterval,
      const T ciProb,
      const IndelPoint indelPoint,
      const multi_array<T, 2> & indelLoglArr,
      const std::vector<multi_array<T, 2> > & indelLoglArrVec,
      const int nIndelElem);

  int EstStrSizePercentilePdBiasCorAccelerated(
      IndelInterval & indelInterval,
      const T ciProb,
      const IndelPoint indelPoint,
      const multi_array<T, 2> & indelLoglArr,
      const std::vector<multi_array<T, 2> > & indelLoglArrVec,
      const int nIndelElem);

  int EstStrSizeHypothesis(
      IndelInterval & indelInterval,
      const T ciProb,
      const IndelPoint indelPoint,
      const multi_array<T, 2> & indelLoglArr,
      const std::vector<multi_array<T, 2> > & indelLoglArrVec,
      const int nIndelElem);

  int FindProbMax(
      double & valMax,
      IndelPoint & indelPoint,
      const multi_array<T, 2> & indelLoglArr,
      const int nIndelElem);

  int FindValMax(
      double & valMax,
      IndelPoint & indelPoint,
      const multi_array<T, 2> & indelLoglArr,
      const multi_array<char, 2> & isVisitedArr,
      const int nIndelElem);

  int GenBootstrap(
      multi_array<T, 2> & indelLoglArrBootstrap,
      const std::vector<multi_array<T, 2> > & indelLoglArrVec,
      const int nIndelElem);

  int GetCiEstimatorRange(int & beg, int & end);
  std::vector<string> GetCiEstimatorVec();
  CiEstimator GetCiEstimator();
  //int GetCiEstimator();
  string GetCiEstimatorStr();
  T GetCiProb();
  int GetNBootstrapRep();
  int GetNVarSampMin();
  inline int GetNStepDownMax(
      const int nIndelStepRight,
      const int refSizeMax,
      const int strLen,
      const int indelStepSize);
  inline int GetNStepLeftMax(
      const int nIndelStepLeft,
      const int insertSize,
      const int libInsertSizeMin,
      const int strLen,
      const int indelStepSize,
      const int anchorLen,
      const int nStepUpMax);
  inline int GetNStepRightMax(
      const int nIndelStepRight,
      const int insertSize,
      const int libInsertSizeMax,
      const int indelStepSize);
  inline int GetNStepUpMax(
      const int nIndelStepLeft,
      const int refSizeMin,
      const int strLen,
      const int indelStepSize);

  int Normalize(
      multi_array<T, 2> & indelLoglArr,
      bool & isNan,
      const multi_array<T, 2> & priorPd,
      const multi_array<char, 2> & isVisitedArr,
      const double valMax,
      const int nIndelElem);

  int SetCiEstimator(CiEstimator ciEstimatorNew);
  int SetCiEstimatorInt(int ciEstimatorNew);
  int SetCiEstimatorStr(const string & ciEstimatorNew);
  int SetCiProb(T ciProbNew);
  int SetNBootstrapRep(int nBootstrapRepNew);
  int SetNVarSampMin(int nVarSampMinNew);

// Protected variables.
protected:
  CiEstimator ciEstimator;
  T ciProb;
  int nBootstrapRep;
  int nVarSampMin;
};

//StrCalc<T>::StrCalc()
StrCalc::StrCalc()
{
  ciEstimator = none;
  ciProb = 0.95;
  nBootstrapRep = 1000;
  nVarSampMin = 2;
}

StrCalc::~StrCalc()
{

}

StrCalc::StrCalc(const StrCalc & rhs)
{
  this->nBootstrapRep = rhs.nBootstrapRep;
  this->nVarSampMin = rhs.nVarSampMin;
  ciEstimator = rhs.ciEstimator;
  ciProb = rhs.ciProb;
}

int StrCalc::CalcCondProb(
    multi_array<T, 2> & indelLoglArr,
    multi_array<char, 2> & isVisitedArr,
    const Read & read,
    const multi_array<T, 2> & mpersPdfArr,
    const vector<int> & indelShift,
    const int center,
    const int strLen,
    const int nIndelStepLeft,
    const int nIndelStepRight,
    const int nIndelStepDown,
    const int nIndelStepUp,
    const int libInsertSizeMin,
    const int libInsertSizeMax,
    const int refSizeMin,
    const int refSizeMax,
    const int anchorLen,
    const int indelStepSize,
    const int nIndelElem,
    const string & libTag,
    const Region & strLoc)
{
  //string func = "SC::CalcCondProb ";

  // Calculate positions.
  int insertSize = read.insertSize;
  int mpersPdfRow = strLen - refSizeMin;
  int mpersPdfCol = insertSize - libInsertSizeMin;
  T probMin = numeric_limits < T > ::max();
  list < pair<int, int> > unvisited;

  for (int i0 = 0; i0 < center + nIndelStepUp + 1; i0++)
  {
    int indelNo1 = indelShift[i0];
    T val1;

    if (mpersPdfRow + indelNo1 < 0
        || mpersPdfRow + indelNo1 > refSizeMax - refSizeMin
        || mpersPdfCol + indelNo1 < 0
        || mpersPdfCol + indelNo1 > libInsertSizeMax - libInsertSizeMin)
    {
      // Store coordinates that will not be visited.
      for (int i1 = 0; i1 < i0 + 1; i1++)
        unvisited.push_back(make_pair(i0, i1));

      continue;
    }
    else
    {
      val1 = mpersPdfArr[mpersPdfRow + indelNo1][mpersPdfCol + indelNo1];
    }

    for (int i1 = 0; i1 < i0 + 1; i1++)
    {
      int indelNo2 = indelShift[i1];
      T val2;
      T valSum;

      if (mpersPdfRow + indelNo2 < 0
          || mpersPdfRow + indelNo2 > refSizeMax - refSizeMin
          || mpersPdfCol + indelNo2 < 0
          || mpersPdfCol + indelNo2 > libInsertSizeMax - libInsertSizeMin)
      {
        // Store coordinates that will not be visited.
        unvisited.push_back(make_pair(i0, i1));

        continue;
      }
      else
      {
        val2 = mpersPdfArr[mpersPdfRow + indelNo2][mpersPdfCol + indelNo2];
      }

      // Compute probability sum of indels.
      valSum = (val1 + val2) / 2;

      if (valSum == 0)
      {

        if (strLoc.refBeg == 40107926)
        {
          cout << read.qname << " " << libTag << " " << strLen << " "
              << insertSize << " " << indelNo1 << " " << indelNo2 << " " << val1
              << " " << val2 << " 3 ";
          cout << (valSum == 0) << endl;
        }

        // Store coordinates that will not be visited.
        unvisited.push_back(make_pair(i0, i1));

        continue;
      }

      // Find minimum probability.
      probMin = min(probMin, valSum);

      indelLoglArr[i0][i1] += log(valSum);
      isVisitedArr[i0][i1] = true;
    }

  }

  T probMinLog = log(probMin);
  for (list<pair<int, int> >::iterator it0 = unvisited.begin();
      it0 != unvisited.end(); it0++)
  {
    indelLoglArr[it0->first][it0->second] = probMinLog;
  }

  return 0;
}

int StrCalc::EstIndelSize(
    vector<multi_array<T, 2> > & indelLoglArrVec,
    multi_array<char, 2> & isVisitedArr,
    map<string*, int> & libCount,
    map<string*, multiset<int> > & libFreq,
    const std::vector<Read> & srpVec,
    const Region strLoc,
    vector<StrLibrary> & libVec,
    const map<string, int> & idLibPosMap,
    const std::vector<int> & indelShift,
    int nIndelElem,
    int indelStepSize,
    int insertSizeMin,
    int refSizeMin,
    int refSizeMax,
    int anchorLen)
{
  int strLen = strLoc.refEnd - strLoc.refBeg + 1;
  int nSrp = int(srpVec.size());
  int center = (nIndelElem - 1) / 2;
  int indelLeftMax = abs(indelShift[0]);
  int indelRightMax = abs(indelShift[nIndelElem - 1]);
  int nIndelStepLeft = indelLeftMax / indelStepSize;
  int nIndelStepRight = indelRightMax / indelStepSize;
  int nIndelStepDown = nIndelStepRight;
  int nIndelStepUp = nIndelStepLeft;
  int libInsertSizeMin = insertSizeMin;

  map<string *, int>::iterator it0;
  map<string *, multiset<int> >::iterator it1;

  for (int i0 = 0; i0 < nSrp; i0++)
  {

    // Find library.
    int libPosNo = srpVec[i0].libPosNo;
    int libInsertSizeMax = libVec[libPosNo].libInsertSizeMax;


    CalcCondProb(indelLoglArrVec[i0], isVisitedArr, srpVec[i0],
        libVec[libPosNo].mpersPdfArr, indelShift, center, strLen,
        nIndelStepLeft, nIndelStepRight, nIndelStepDown, nIndelStepUp,
        libInsertSizeMin, libInsertSizeMax, refSizeMin, refSizeMax, anchorLen,
        indelStepSize, nIndelElem, libVec[libPosNo].libTag, strLoc);

    // Save library frequencies.
    it0 = libCount.find(&libVec[libPosNo].libTag);
    it0->second += 1;
    it1 = libFreq.find(&libVec[libPosNo].libTag);
    it1->second.insert(srpVec[i0].insertSize);
  }

  return 0;
}

int StrCalc::EstStrSize(
    multi_array<T, 2> & indelLoglEstArr,
    const std::vector<multi_array<T, 2> > & indelLoglArr,
    int nSrp,
    int nIndelElem)
{
  // Compute posterior probabilities (read-pair probabilities).
  for (int i0 = 0; i0 < nSrp; i0++)
  {
    for (int i1 = 0; i1 < nIndelElem; i1++)
    {
      for (int i2 = 0; i2 < i1 + 1; i2++)
      {
        indelLoglEstArr[i1][i2] += indelLoglArr[i0][i1][i2];
      }
    }
  }

  return 0;
}

// Confidence interval using points of indels.
int StrCalc::EstStrSizeBootstrap1(
    vector<IndelPoint> & indelPointVec,
    const T ciProb,
    const IndelPoint indelPoint,
    const multi_array<T, 2> & indelLoglArr,
    const int nIndelElem)
{
  return 0;
}

// Confidence interval (credible interval).
int StrCalc::EstStrSizeCredibleInterval(
    IndelInterval & indelInterval,
    const T ciProb,
    const IndelPoint indelPoint,
    const multi_array<T, 2> & indelLoglArr,
    const int nIndelElem)
{
  vector<T> probArr(nIndelElem);
  vector<T> probCumArr(nIndelElem);

  // Save the points of the maximum.
  indelInterval.i1 = indelPoint.i1;
  indelInterval.i2 = indelPoint.i2;
  int indelSize[2] = { indelPoint.i2, indelPoint.i1 };
  int indelLow[2];
  int indelHigh[2];
  T probLow = (1 - ciProb) / 2;
  T probHigh = 1 - probLow;

  // Integrate to find interval of alleles.
  for (int i0 = 0; i0 < 2; i0++)
  {
    if (i0 == 1 && indelSize[0] == indelSize[1])
    {
      indelLow[1] = indelLow[0];
      indelHigh[1] = indelHigh[0];
      break;
    }

    int indel = indelSize[i0];
    T probCum = 0;

    // Get row.
    for (int i1 = 0; i1 < indel + 1; i1++)
    {
      probArr[i1] = indelLoglArr[indel][i1];
      probCum += indelLoglArr[indel][i1];
      probCumArr[i1] = probCum;
    }

    // Get column.
    for (int i1 = indel + 1; i1 < nIndelElem; i1++)
    {
      probArr[i1] = indelLoglArr[i1][indel];
      probCum += indelLoglArr[i1][indel];
      probCumArr[i1] = probCum;
    }

    // Normalize values.
    for (int i1 = 0; i1 < nIndelElem; i1++)
    {
      probArr[i1] /= probCum;
      probCumArr[i1] /= probCum;
    }

    // Find confidence interval (low).
    for (int i1 = 0; i1 < nIndelElem; i1++)
    {
      if (probCumArr[i1] >= probLow)
      {
        indelLow[i0] = i1;
        break;
      }
    }

    // Find confidence interval (high).
    for (int i1 = 0; i1 < nIndelElem; i1++)
    {
      if (probCumArr[i1] >= probHigh)
      {
        indelHigh[i0] = i1;
        break;
      }
    }
  }

  // Save values.
  indelInterval.i1Low = indelLow[0];
  indelInterval.i1High = indelHigh[0];
  indelInterval.i2Low = indelLow[1];
  indelInterval.i2High = indelHigh[1];

  return 0;
}

// Warning: isVisitedArr has not been implemented properly.
// Confidence interval (bootstrap, percentile of PDF).
int StrCalc::EstStrSizePercentilePd(
    vector<IndelPoint> & indelPointVec,
    const T ciProb,
    const IndelPoint indelPoint,
    const multi_array<T, 2> & indelLoglArr,
    const std::vector<multi_array<T, 2> > & indelLoglArrVec,
    const multi_array<char, 2> & isVisitedArr,
    const multi_array<T, 2> & priorPd,
    const int nIndelElem)
{
  set<pair<int, int> >::iterator it0;

  vector < multi_array<T, 2> > indelLoglArrBootstrap(nBootstrapRep);
  vector<T> pointMax(nBootstrapRep);
  vector<T> diff(nBootstrapRep);
  set < pair<int, int> > pointRecent;
  set < pair<int, int> > pointPrev;
  set < pair<int, int> > indelPointSet;

  // Initialize values.
  bool isNan;
  //IndelPoint indelPoint_cur = indelPoint;
  T propCritical = 1 - ciProb;

  // Initialize random seed.
  srand(0);

  // Calculate confidence interval.
  for (int i0 = 0; i0 < nBootstrapRep; i0++)
  {
    IndelPoint indelPoint_cur;
    T valMax = -1e10;

    indelLoglArrBootstrap[i0].resize(extents[nIndelElem][nIndelElem]);

    // Generate a bootstrap replicate.
    GenBootstrap(indelLoglArrBootstrap[i0], indelLoglArrVec, nIndelElem);

    // Find greatest nonzero posterior probability.
    FindValMax(valMax, indelPoint_cur, indelLoglArrBootstrap[i0], isVisitedArr,
        nIndelElem);

    // Exponentiate and normalize posterior probabilities, while avoiding numerical underflow
    // and overflow.
    Normalize(indelLoglArrBootstrap[i0], isNan, priorPd, isVisitedArr, valMax,
        nIndelElem);
  }

  // Compare indel sizes to that of the highest probability.

  // Insert the current highest probability.
  pointPrev.insert(make_pair(indelPoint.i1, indelPoint.i2));
  pointRecent.insert(make_pair(indelPoint.i1, indelPoint.i2));
  indelPointVec.push_back(indelPoint);

  // Include
  for (int i0 = 0; i0 < nBootstrapRep; i0++)
    pointMax[i0] = indelLoglArrBootstrap[i0][indelPoint.i1][indelPoint.i2];

  while (true)
  {

    std::vector<T>::iterator it1;

    set < pair<int, int> > pointNew;

    bool found = false;

    for (it0 = pointRecent.begin(); it0 != pointRecent.end(); it0++)
    {
      int pos1 = it0->first;
      int pos2 = it0->second;
      T rankProp;

      for (int i1 = pos1 - 1; i1 < pos1 + 2; i1++)
      {
        if (i1 < 0 || i1 >= nIndelElem)
          continue;

        for (int i2 = pos2 - 1; i2 < pos2 + 2; i2++)
        {
          if (pointPrev.find(make_pair(i1, i2)) != pointPrev.end() || i2 < 0
              || i2 >= nIndelElem)
          {
            continue;
          }

          // This indel has now been checked once, do not visit it again.
          pointPrev.insert(make_pair(i1, i2));
          pointNew.insert(make_pair(i1, i2));

          // Calculate difference.
          for (int i3 = 0; i3 < nBootstrapRep; i3++)
          {
            diff[i3] = pointMax[i3] - indelLoglArrBootstrap[i3][i1][i2];
          }

          // Sort values.
          sort(diff.begin(), diff.end());

          // Find proportion that is smaller than zero.
          it1 = lower_bound(diff.begin(), diff.end(), 0);
          rankProp = (it1 - diff.begin() + 0.5) / T(nBootstrapRep);

          // Check whether proportion is statistically insignificant.
          if (rankProp < propCritical)
            continue;

          // If statistically significant, save indel coordinate.
          indelPointSet.insert(make_pair(i1, i2));
          found = true;

        }
      }
    }

    pointRecent = pointNew;

    if (found == false)
      break;
  }

  // Copy indels that have been found.
  indelPointVec.clear();
  indelPointVec.reserve(indelPointSet.size() + 1);
  indelPointSet.insert(make_pair(indelPoint.i1, indelPoint.i2));

  for (it0 = indelPointSet.begin(); it0 != indelPointSet.end(); it0++)
  {
    //indelPointVec.push_back(IndelPoint{it0->first, it0->second});
    IndelPoint indelPoint;
    indelPoint.i1 = it0->first;
    indelPoint.i2 = it0->second;
    indelPointVec.push_back(indelPoint);
  }

  return 0;
}

// Confidence interval (bootstrap, percentile of error PDF).
int StrCalc::EstStrSizePercentilePdError(
    IndelInterval & indelInterval,
    const T ciProb,
    const IndelPoint indelPoint,
    const multi_array<T, 2> & indelLoglArr,
    const std::vector<multi_array<T, 2> > & indelLoglArrVec,
    const int nIndelElem)
{
  vector<T> probArr(nIndelElem);
  vector<T> probCumArr(nIndelElem);

  // Save the points of the maximum.
  indelInterval.i1 = indelPoint.i1;
  indelInterval.i2 = indelPoint.i2;

  for (int i0 = 0; i0 < nIndelElem; i0++)
  {

  }

  return 0;
}

// Confidence interval (bootstrap, percentile of PDF, bias-corrected).
int StrCalc::EstStrSizePercentilePdBiasCor(
    IndelInterval & indelInterval,
    const T ciProb,
    const IndelPoint indelPoint,
    const multi_array<T, 2> & indelLoglArr,
    const std::vector<multi_array<T, 2> > & indelLoglArrVec,
    const int nIndelElem)
{
  vector<T> probArr(nIndelElem);
  vector<T> probCumArr(nIndelElem);

  // Save the points of the maximum.
  indelInterval.i1 = indelPoint.i1;
  indelInterval.i2 = indelPoint.i2;

  for (int i0 = 0; i0 < nIndelElem; i0++)
  {

  }

  return 0;
}

// Confidence interval (bootstrap, percentile of PDF, accelerated bias-corrected).
int StrCalc::EstStrSizePercentilePdBiasCorAccelerated(
    IndelInterval & indelInterval,
    const T ciProb,
    const IndelPoint indelPoint,
    const multi_array<T, 2> & indelLoglArr,
    const std::vector<multi_array<T, 2> > & indelLoglArrVec,
    const int nIndelElem)
{
  vector<T> probArr(nIndelElem);
  vector<T> probCumArr(nIndelElem);

  // Save the points of the maximum.
  indelInterval.i1 = indelPoint.i1;
  indelInterval.i2 = indelPoint.i2;

  for (int i0 = 0; i0 < nIndelElem; i0++)
  {

  }

  return 0;
}

// Confidence interval (bootstrap, hypothesis).
int StrCalc::EstStrSizeHypothesis(
    IndelInterval & indelInterval,
    const T ciProb,
    const IndelPoint indelPoint,
    const multi_array<T, 2> & indelLoglArr,
    const std::vector<multi_array<T, 2> > & indelLoglArrVec,
    const int nIndelElem)
{
  vector<T> probArr(nIndelElem);
  vector<T> probCumArr(nIndelElem);

  // Save the points of the maximum.
  indelInterval.i1 = indelPoint.i1;
  indelInterval.i2 = indelPoint.i2;

  for (int i0 = 0; i0 < nIndelElem; i0++)
  {

  }

  return 0;
}

// Function to find maximum probability.
int StrCalc::FindProbMax(
    double & valMax,
    IndelPoint & indelPoint,
    const multi_array<T, 2> & indelLoglArr,
    const int nIndelElem)
{
  for (int i0 = 0; i0 < nIndelElem; i0++)
  {
    for (int i1 = 0; i1 < i0 + 1; i1++)
    {
      T val = indelLoglArr[i0][i1];

      if (valMax < val)
      {
        valMax = val;
        indelPoint.i1 = i0;
        indelPoint.i2 = i1;
      }
    }
  }

  return 0;
}

// Function to find the maximum value amongst the posterior probabilities.
int StrCalc::FindValMax(
    double & valMax,
    IndelPoint & indelPoint,
    const multi_array<T, 2> & indelLoglArr,
    const multi_array<char, 2> & isVisitedArr,
    const int nIndelElem)
{
  for (int i0 = 0; i0 < nIndelElem; i0++)
  {
    for (int i1 = 0; i1 < i0 + 1; i1++)
    {
      T val = indelLoglArr[i0][i1];

      if (isVisitedArr[i0][i1] == false)
        continue;

      if (valMax < val)
      {
        valMax = val;
        indelPoint.i1 = i0;
        indelPoint.i2 = i1;
      }
    }
  }

  return 0;
}

int StrCalc::GenBootstrap(
    multi_array<T, 2> & indelLoglArrBootstrap,
    const std::vector<multi_array<T, 2> > & indelLoglArrVec,
    const int nIndelElem)
{
  int nSamp = int(indelLoglArrVec.size());

  for (int i0 = 0; i0 < nSamp; i0++)
  {
    int sampRand = rand() % nSamp;

    for (int i1 = 0; i1 < nIndelElem; i1++)
      for (int i2 = 0; i2 < i1 + 1; i2++)
        indelLoglArrBootstrap[i1][i2] += indelLoglArrVec[sampRand][i1][i2];
  }

  return 0;
}

T StrCalc::GetCiProb()
{
  return ciProb;
}

CiEstimator StrCalc::GetCiEstimator()
{
  return ciEstimator;
}

string StrCalc::GetCiEstimatorStr()
{
  string ciEstimatorStr;

  if (ciEstimator == none)
    ciEstimatorStr = "none";
  else if (ciEstimator == Credible_Interval)
    ciEstimatorStr = "Credible Interval";
  else if (ciEstimator == Percentile_PDF)
    ciEstimatorStr = "Percentile PDF";

  return ciEstimatorStr;
}

int StrCalc::GetCiEstimatorRange(int & beg, int & end)
{
  beg = 0;
  end = 2;

  return 0;
}

std::vector<string> StrCalc::GetCiEstimatorVec()
{
  std::vector < string > ciEstimatorVec;
  ciEstimatorVec.push_back("none");
  ciEstimatorVec.push_back("Credible Interval");
  ciEstimatorVec.push_back("Percentile PDF");

  return ciEstimatorVec;
}

int StrCalc::GetNBootstrapRep()
{
  return nBootstrapRep;
}

int StrCalc::GetNVarSampMin()
{
  return nVarSampMin;
}

inline int StrCalc::GetNStepDownMax(
    const int nIndelStepRight,
    const int refSizeMax,
    const int strLen,
    const int indelStepSize)
{
  return min(nIndelStepRight, (refSizeMax - strLen) / indelStepSize);
}

inline int StrCalc::GetNStepLeftMax(
    const int nIndelStepLeft,
    const int insertSize,
    const int libInsertSizeMin,
    const int strLen,
    const int indelStepSize,
    const int anchorLen,
    const int nStepUpMax)
{
  int nStepLeftMax = min(nIndelStepLeft,
      (insertSize - strLen - 2 * anchorLen) / indelStepSize);

  nStepLeftMax = min(nStepLeftMax, nStepUpMax);

  return nStepLeftMax;
}

inline int StrCalc::GetNStepRightMax(
    const int nIndelStepRight,
    const int insertSize,
    const int libInsertSizeMax,
    const int indelStepSize)
{
  return min(nIndelStepRight, (libInsertSizeMax - insertSize) / indelStepSize);
}

inline int StrCalc::GetNStepUpMax(
    const int nIndelStepLeft,
    const int refSizeMin,
    const int strLen,
    const int indelStepSize)
{
  int mpersPdfRow = strLen - refSizeMin;
  return min(nIndelStepLeft, mpersPdfRow / indelStepSize);
}

// Exponentiate and normalize posterior probabilities, while avoiding numerical underflow and
// overflow.
int StrCalc::Normalize(
    multi_array<T, 2> & indelLoglArr,
    bool & isNan,
    const multi_array<T, 2> & priorPd,
    const multi_array<char, 2> & isVisitedArr,
    const double valMax,
    const int nIndelElem)
{
  T sum = 0;
  const T minLimit = 1e1 * numeric_limits < T > ::min();

  // Exponentiate values that are nonzero.
  for (int i0 = 0; i0 < nIndelElem; i0++)
  {
    for (int i1 = 0; i1 < i0 + 1; i1++)
    {
      T val = indelLoglArr[i0][i1];

      //if (abs(val) < minLimit)
      if (isVisitedArr[i0][i1] == false)
      {
        indelLoglArr[i0][i1] = 0;
        continue;
      }

      T valPrior = priorPd[i0][i1];

      val = exp(val - valMax + valPrior);
      sum += val;
      indelLoglArr[i0][i1] = val;
    }
  }

  // If there is no sum at all, then no spanning read-pairs have been used,
  // and all probabilities are zero.
  if (abs(sum) < minLimit)
  {
    isNan = true;
    return 1;
  }

  // Normalize values.
  for (int i0 = 0; i0 < nIndelElem; i0++)
  {
    for (int i1 = 0; i1 < i0 + 1; i1++)
    {
      indelLoglArr[i0][i1] /= sum;
    }
  }

  return 0;
}

int StrCalc::SetCiEstimator(CiEstimator ciEstimatorNew)
{
  ciEstimator = ciEstimatorNew;

  return 0;
}

int StrCalc::SetCiEstimatorInt(int ciEstimatorNew)
{
  if (ciEstimatorNew < 0 || ciEstimatorNew > 2)
  {
    cout << "Accepted confidence interval estimator IDs are:" << endl;
    cout << " 0: \"none\"," << endl;
    cout << " 1: \"Credible Interval\"," << endl;
    cout << " 2: \"Percentile PDF\"," << endl;

    return 1;
  }

  switch (ciEstimatorNew)
  {
  case 0:
    ciEstimator = none;
    break;
  case 1:
    ciEstimator = Credible_Interval;
    break;
  case 2:
    ciEstimator = Percentile_PDF;
    break;
  }

  return 0;
}

int StrCalc::SetCiEstimatorStr(const string & ciEstimatorNew)
{
  if (ciEstimatorNew == "none")
  {
    ciEstimator = none;
  }
  else if (ciEstimatorNew == "Credible Interval")
  {
    ciEstimator = Credible_Interval;
  }
  else if (ciEstimatorNew == "Percentile PDF")
  {
    ciEstimator = Percentile_PDF;
  }
  else
  {
    cout << "Accepted confidence interval estimators are:" << endl;
    cout << "  \"none\"," << endl;
    cout << "  \"Credible Interval\"," << endl;
    cout << "  \"Percentile PDF\"," << endl;
    return 1;
  }

  return 0;
}

int StrCalc::SetCiProb(T ciProbNew)
{
  if (ciProb <= 0 || ciProb >= 1)
  {
    cout
        << "Confidence interval estimator probability has the valid range 0 < x < 1."
        << endl;
    return 1;
  }

  ciProb = ciProbNew;

  return 0;
}

int StrCalc::SetNBootstrapRep(int nBootstrapRepNew)
{
  if (nBootstrapRepNew < 1)
  {
    cout << "The number of bootstrap replicates must be at least 1." << endl;
    return 1;
  }

  this->nBootstrapRep = nBootstrapRepNew;

  return 0;
}

int StrCalc::SetNVarSampMin(int nVarSampMinNew)
{
  if (nVarSampMinNew < 2)
  {
    cout << "The minimum sample size for variance estimate must be at least 2."
        << endl;
    return 1;
  }

  this->nVarSampMin = nVarSampMinNew;

  return 0;
}

#endif // STR_CALC_H
