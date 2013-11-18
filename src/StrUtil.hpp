#ifndef STRUTIL_HPP_
#define STRUTIL_HPP_

// STL libraries
#include <algorithm>    // swap(), max_element
#include <cmath>        // log(), exp()
#include <cstdlib>      // abs(), itoa()
#include <deque>
#include <fstream>      // open(), close()
#include <iostream>
#include <limits>
#include <string>
#include <vector>

// External libraries
#include "api/BamAlignment.h"
#include <boost/lexical_cast.hpp>
#include <boost/multi_array.hpp>
#include <boost/numeric/ublas/triangular.hpp>

// Namespace
using namespace std;
using namespace boost;

// Class
class StrUtil
{
public:
  StrUtil();
  ~StrUtil();
  StrUtil(const StrUtil &rhs);

  int AnchorOrFn(
      int &anchorLD,
      int &anchorLRWC,
      int &anchorLRWoC,
      int &anchorRD,
      int &anchorRRWC,
      int &anchorRRWoC,
      const Region &strLocDon,
      const Region &rLeftDon,
      const Region &rRightDon,
      const Region &strLocRef,
      const ReadPair &orp);
  int AnchorOrFp(
      int &anchorLD,
      int &anchorLRWC,
      int &anchorLRWoC,
      int &anchorRD,
      int &anchorRRWC,
      int &anchorRRWoC,
      const Region &strLocDon,
      const Region &rLeftDon,
      const Region &rRightDon,
      const Region &strLocRef,
      const ReadPair &orp);
  int AnchorOrTp(
      int &anchorLD,
      int &anchorLRWC,
      int &anchorLRWoC,
      int &anchorRD,
      int &anchorRRWC,
      int &anchorRRWoC,
      const Region &strLocDon,
      const Region &rLeftDon,
      const Region &rRightDon,
      const Region &strLocRef,
      const ReadPair &orp);
  int AnchorSrpFn(
      int &anchorLD,
      int &anchorLRWC,
      int &anchorLRWoC,
      int &anchorRD,
      int &anchorRRWC,
      int &anchorRRWoC,
      const Region &strLocDon,
      const Region &rpDon,
      const Region &rLeftDon,
      const Region &rRightDon,
      const Region &strLocRef,
      const ReadPair &orp);
  int AnchorSrpFp(
      int &anchorLD,
      int &anchorLRWC,
      int &anchorLRWoC,
      int &anchorRD,
      int &anchorRRWC,
      int &anchorRRWoC,
      const Region &strLocDon,
      const Region &rpDon,
      const Region &rLeftDon,
      const Region &rRightDon,
      const Region &strLocRef,
      const ReadPair &orp);
  int AnchorSrpTp(
      int &anchorLD,
      int &anchorLRWC,
      int &anchorLRWoC,
      int &anchorRD,
      int &anchorRRWC,
      int &anchorRRWoC,
      const Region &strLocDon,
      const Region &rpDon,
      const Region &rLeftDon,
      const Region &rRightDon,
      const Region &strLocRef,
      const ReadPair &orp);

  int ConvCoordTo0(Region &reg);
  int ConvCoordTo1(Region &reg);

  int FindNeighborNearest(
      int &strNearestOnLeft,
      int &strNearestOnRight,
      const Region &strLoc,
      const map<string, pair<std::vector<int>, std::vector<int> > > &neighborMap);

  string GenCigar(const std::vector<CigarOp> &cigarData);
  int GenNeighborMap(
      map<string, pair<std::vector<int>, std::vector<int> > > &neighborMap,
      const deque<Region> &strLocDeque,
      const set<string> &refIdSet);

  int GetAnchorLen();
  int GetLibPos(const Read &r, const map<string, int> &idLibPosMap);
  int GetRLen(const Read &r);

  bool IsCovered(const Region &read, const Region &strLoc);

  int QnameToRpPosMaq(Region &reg, const string &qname);

  int ReadLenWithClipping(const std::vector<CigarOp> &cigarData, int len);
  int ReadLenWithoutClipping(const std::vector<CigarOp> &cigarData, int len);
  int RectifyPosWithClip(Region &reg, const Read &read);
  int RectifyPosWoClip(Region &reg, const Read &read);

  int RpPosToRPosMaq(Region &regRead, const Region &regRp, int len, int mpers);

  int SetAnchorLen(int anchorLenNew);
protected:

private:
  int AnchorOrDonFp(
      int &anchorLD,
      int &anchorRD,
      const Region &strLocDon,
      const Region &rLeftDon,
      const Region &rRightDon,
      const ReadPair &orp);
  int AnchorOrDonTpFn(
      int &anchorLD,
      int &anchorRD,
      const Region &strLocDon,
      const Region &rLeftDon,
      const Region &rRightDon);
  int AnchorOrRefFn(
      int &anchorLR,
      int &anchorRR,
      const Region &rRef,
      const Region &strLocRef);
  int AnchorOrRefTpFp(
      int &anchorLR,
      int &anchorRR,
      const Region &rRef,
      const Region &strLocRef);
  int AnchorSrpDonFp(
      int &anchorLD,
      int &anchorRD,
      const Region &strLocDon,
      const Region &rLeftDon,
      const Region &rRightDon,
      const ReadPair &orp);
  int AnchorSrpDonTpFn(
      int &anchorLD,
      int &anchorRD,
      const Region &strLocDon,
      const Region &rpDon,
      const Region &rLeftDon,
      const Region &rRightDon);
  int AnchorSrpRefFn(
      int &anchorLR,
      int &anchorRR,
      const Region &rLRef,
      const Region &rRRef,
      const Region &strLocRef);
  int AnchorSrpRefTpFp(
      int &anchorLR,
      int &anchorRR,
      const Region &rLRef,
      const Region &rRRef,
      const Region &strLocRef);

  int anchorLen;
};

StrUtil::StrUtil()
{
  anchorLen = 1;
}

StrUtil::~StrUtil()
{
}

int StrUtil::AnchorOrDonFp(
    int &anchorLD,
    int &anchorRD,
    const Region &strLocDon,
    const Region &rLeftDon,
    const Region &rRightDon,
    const ReadPair &orp)
{
  int len = orp.r1.length;

  if (abs(rLeftDon.refBeg - strLocDon.refBeg)
      < abs(rRightDon.refBeg - strLocDon.refBeg))
  {
    anchorLD = max(strLocDon.refBeg - rLeftDon.refBeg, 0);
    anchorLD = min(anchorLD, len);
    anchorRD = max(rLeftDon.refEnd - strLocDon.refEnd, 0);
    anchorRD = min(anchorRD, len);
  }
  else
  {
    anchorLD = max(strLocDon.refBeg - rRightDon.refBeg, 0);
    anchorLD = min(anchorLD, len);
    anchorRD = max(rRightDon.refEnd - strLocDon.refEnd, 0);
    anchorRD = min(anchorRD, len);
  }

  return 0;
}

int StrUtil::AnchorOrDonTpFn(
    int &anchorLD,
    int &anchorRD,
    const Region &strLocDon,
    const Region &rLeftDon,
    const Region &rRightDon)
{
  // TP/FN.
  if (IsCovered(rLeftDon, strLocDon))
  {
    anchorLD = strLocDon.refBeg - rLeftDon.refBeg;
    anchorRD = rLeftDon.refEnd - strLocDon.refEnd;
  }
  // TP/FN.
  else if (IsCovered(rRightDon, strLocDon))
  {
    anchorLD = strLocDon.refBeg - rRightDon.refBeg;
    anchorRD = rRightDon.refEnd - strLocDon.refEnd;
  }

  return 0;
}

int StrUtil::AnchorOrRefFn(
    int &anchorLR,
    int &anchorRR,
    const Region &rRef,
    const Region &strLocRef)
{
  int len = rRef.refEnd - rRef.refBeg + 1;

  anchorLR = max(strLocRef.refBeg - rRef.refBeg, 0);
  anchorLR = min(anchorLR, len);
  anchorRR = max(rRef.refEnd - strLocRef.refEnd, 0);
  anchorRR = min(anchorRR, len);

  return 0;
}

int StrUtil::AnchorOrRefTpFp(
    int &anchorLR,
    int &anchorRR,
    const Region &rRef,
    const Region &strLocRef)
{
  anchorLR = max(strLocRef.refBeg - rRef.refBeg, 0);
  anchorRR = max(rRef.refEnd - strLocRef.refEnd, 0);

  return 0;
}

int StrUtil::AnchorOrFn(
    int &anchorLD,
    int &anchorLRWC,
    int &anchorLRWoC,
    int &anchorRD,
    int &anchorRRWC,
    int &anchorRRWoC,
    const Region &strLocDon,
    const Region &rLeftDon,
    const Region &rRightDon,
    const Region &strLocRef,
    const ReadPair &orp)
{
  Region rOverlap;

  AnchorOrDonTpFn(anchorLD, anchorRD, strLocDon, rLeftDon, rRightDon);
  RectifyPosWithClip(rOverlap, orp.r1);
  AnchorOrRefFn(anchorLRWC, anchorRRWC, rOverlap, strLocRef);
  RectifyPosWoClip(rOverlap, orp.r1);
  AnchorOrRefFn(anchorLRWoC, anchorRRWoC, rOverlap, strLocRef);

  return 0;
}

int StrUtil::AnchorOrFp(
    int &anchorLD,
    int &anchorLRWC,
    int &anchorLRWoC,
    int &anchorRD,
    int &anchorRRWC,
    int &anchorRRWoC,
    const Region &strLocDon,
    const Region &rLeftDon,
    const Region &rRightDon,
    const Region &strLocRef,
    const ReadPair &orp)
{
  Region rOverlap;

  AnchorOrDonFp(anchorLD, anchorRD, strLocDon, rLeftDon, rRightDon, orp);
  RectifyPosWithClip(rOverlap, orp.r1);
  AnchorOrRefTpFp(anchorLRWC, anchorRRWC, rOverlap, strLocRef);
  RectifyPosWoClip(rOverlap, orp.r1);
  AnchorOrRefTpFp(anchorLRWoC, anchorRRWoC, rOverlap, strLocRef);

  return 0;
}

int StrUtil::AnchorOrTp(
    int &anchorLD,
    int &anchorLRWC,
    int &anchorLRWoC,
    int &anchorRD,
    int &anchorRRWC,
    int &anchorRRWoC,
    const Region &strLocDon,
    const Region &rLeftDon,
    const Region &rRightDon,
    const Region &strLocRef,
    const ReadPair &orp)
{
  Region rOverlap;

  AnchorOrDonTpFn(anchorLD, anchorRD, strLocDon, rLeftDon, rRightDon);
  RectifyPosWithClip(rOverlap, orp.r1);
  AnchorOrRefTpFp(anchorLRWC, anchorRRWC, rOverlap, strLocRef);
  RectifyPosWoClip(rOverlap, orp.r1);
  AnchorOrRefTpFp(anchorLRWoC, anchorRRWoC, rOverlap, strLocRef);

  return 0;
}

int StrUtil::AnchorSrpFn(
    int &anchorLD,
    int &anchorLRWC,
    int &anchorLRWoC,
    int &anchorRD,
    int &anchorRRWC,
    int &anchorRRWoC,
    const Region &strLocDon,
    const Region &rpDon,
    const Region &rLeftDon,
    const Region &rRightDon,
    const Region &strLocRef,
    const ReadPair &srp)
{
  Region rLRef;
  Region rRRef;

  AnchorSrpDonTpFn(anchorLD, anchorRD, strLocDon, rpDon, rLeftDon, rRightDon);

  RectifyPosWithClip(rLRef, srp.r1);
  RectifyPosWithClip(rRRef, srp.r2);
  AnchorSrpRefFn(anchorLRWC, anchorRRWC, rLRef, rRRef, strLocRef);

  RectifyPosWoClip(rLRef, srp.r1);
  RectifyPosWoClip(rRRef, srp.r2);
  AnchorSrpRefFn(anchorLRWoC, anchorRRWoC, rLRef, rRRef, strLocRef);

  return 0;
}

int StrUtil::AnchorSrpFp(
    int &anchorLD,
    int &anchorLRWC,
    int &anchorLRWoC,
    int &anchorRD,
    int &anchorRRWC,
    int &anchorRRWoC,
    const Region &strLocDon,
    const Region &rpDon,
    const Region &rLeftDon,
    const Region &rRightDon,
    const Region &strLocRef,
    const ReadPair &srp)
{
  Region rLRef;
  Region rRRef;

  AnchorSrpDonFp(anchorLD, anchorRD, strLocDon, rLeftDon, rRightDon, srp);

  RectifyPosWithClip(rLRef, srp.r1);
  RectifyPosWithClip(rRRef, srp.r2);
  AnchorSrpRefTpFp(anchorLRWC, anchorRRWC, rLRef, rRRef, strLocRef);

  RectifyPosWoClip(rLRef, srp.r1);
  RectifyPosWoClip(rRRef, srp.r2);
  AnchorSrpRefTpFp(anchorLRWoC, anchorRRWoC, rLRef, rRRef, strLocRef);

  return 0;
}

int StrUtil::AnchorSrpTp(
    int &anchorLD,
    int &anchorLRWC,
    int &anchorLRWoC,
    int &anchorRD,
    int &anchorRRWC,
    int &anchorRRWoC,
    const Region &strLocDon,
    const Region &rpDon,
    const Region &rLeftDon,
    const Region &rRightDon,
    const Region &strLocRef,
    const ReadPair &srp)
{
  Region rLRef;
  Region rRRef;

  AnchorSrpDonTpFn(anchorLD, anchorRD, strLocDon, rpDon, rLeftDon, rRightDon);

  RectifyPosWithClip(rLRef, srp.r1);
  RectifyPosWithClip(rRRef, srp.r2);
  AnchorSrpRefTpFp(anchorLRWC, anchorRRWC, rLRef, rRRef, strLocRef);

  RectifyPosWoClip(rLRef, srp.r1);
  RectifyPosWoClip(rRRef, srp.r2);
  AnchorSrpRefTpFp(anchorLRWoC, anchorRRWoC, rLRef, rRRef, strLocRef);

  return 0;
}

int StrUtil::AnchorSrpDonFp(
    int &anchorLD,
    int &anchorRD,
    const Region &strLocDon,
    const Region &rLeftDon,
    const Region &rRightDon,
    const ReadPair &srp)
{
  int rLen = rLeftDon.refEnd - rLeftDon.refBeg + 1;

  anchorLD = max(strLocDon.refBeg - rLeftDon.refBeg, 0);
  if (anchorLD > rLen)
    anchorLD = rLen;

  rLen = rRightDon.refEnd - rRightDon.refBeg + 1;
  anchorRD = max(rRightDon.refEnd - strLocDon.refEnd, 0);
  if (anchorRD > rLen)
    anchorRD = rLen;

  return 0;
}

int StrUtil::AnchorSrpDonTpFn(
    int &anchorLD,
    int &anchorRD,
    const Region &strLocDon,
    const Region &rpDon,
    const Region &rLeftDon,
    const Region &rRightDon)
{
  // TP/FN.
  if (IsCovered(rpDon, strLocDon))
  {
    anchorLD = min(rLeftDon.refEnd + 1, strLocDon.refBeg) - rLeftDon.refBeg;
    anchorRD = rRightDon.refEnd - max(rRightDon.refBeg - 1, strLocDon.refEnd);
  }

  return 0;
}

int StrUtil::AnchorSrpRefFn(
    int &anchorLR,
    int &anchorRR,
    const Region &rLRef,
    const Region &rRRef,
    const Region &strLocRef)
{
  // Anchor in reference genome, left read.
  anchorLR = min(rLRef.refEnd + 1, strLocRef.refBeg) - rLRef.refBeg;
  if (anchorLR < 0)
    anchorLR = 0;

  // Anchor in reference genome, right read.
  anchorRR = rRRef.refEnd - max(rRRef.refBeg - 1, strLocRef.refEnd);
  if (anchorRR < 0)
    anchorRR = 0;

  return 0;
}

int StrUtil::AnchorSrpRefTpFp(
    int &anchorLR,
    int &anchorRR,
    const Region &rLRef,
    const Region &rRRef,
    const Region &strLocRef)
{
  // Anchor in reference genome, left read.
  anchorLR = min(rLRef.refEnd + 1, strLocRef.refBeg) - rLRef.refBeg;

  // Anchor in reference genome, right read.
  anchorRR = rRRef.refEnd - max(rRRef.refBeg - 1, strLocRef.refEnd);

  return 0;
}

inline int StrUtil::ConvCoordTo0(Region &reg)
{
  reg.refBeg--;
  reg.refEnd--;
  return 0;
}

int StrUtil::ConvCoordTo1(Region &reg)
{
  reg.refBeg++;
  reg.refEnd++;
  return 0;
}

int StrUtil::FindNeighborNearest(
    int &strNearestOnLeft,
    int &strNearestOnRight,
    const Region &strLoc,
    const map<string, pair<std::vector<int>, std::vector<int> > > &neighborMap)
{
  map<string, pair<std::vector<int>, std::vector<int> > >::const_iterator cit0 =
      neighborMap.find(strLoc.refIdStr);
  vector<int>::const_iterator cit1;

  // This should not need to be run, as "read" is on the same chromosome as "strLoc".
  //if (cit0 == neighborMap.end())
  //    return 1;

  // Find the nearest position to the left of the STR locus under investigation.
  cit1 = lower_bound(cit0->second.second.begin(), cit0->second.second.end(),
      strLoc.refBeg - 1);

  while (*cit1 >= strLoc.refBeg && cit1 != cit0->second.second.begin())
    cit1--;

  if (cit1 == cit0->second.second.begin())
    strNearestOnLeft = -1;
  else
    strNearestOnLeft = *cit1;

  // Find the nearest position to the right of the STR locus under investigation.
  cit1 = upper_bound(cit0->second.first.begin(), cit0->second.first.end(),
      strLoc.refEnd);

  if (cit1 == cit0->second.first.end())
    strNearestOnRight = numeric_limits<int>::max();
  else
    strNearestOnRight = *cit1;

  return 0;
}

string StrUtil::GenCigar(const std::vector<CigarOp> &cigarData)
{
  string cigar;
  int cigLen = int(cigarData.size());

  for (int i0 = 0; i0 < cigLen; i0++)
  {
    cigar += lexical_cast<string>(cigarData[i0].Length) + cigarData[i0].Type;
  }

  return cigar;
}

int StrUtil::GenNeighborMap(
    map<string, pair<std::vector<int>, std::vector<int> > > &neighborMap,
    const deque<Region> &strLocDeque,
    const set<string> &refIdSet)
{
  set<string>::const_iterator cit0;
  map<string, pair<std::vector<int>, std::vector<int> > >::iterator it1;
  std::vector<int> intVec;
  int nStrLoc = int(strLocDeque.size());

  intVec.reserve(nStrLoc);

  // Find chromosomes.
  for (cit0 = refIdSet.begin(); cit0 != refIdSet.end(); cit0++)
  {
    neighborMap.insert(make_pair(*cit0, make_pair(intVec, intVec)));
  }

  // Add STR locus coordinates to value (vector) of appropriate key (chromosome).
  for (int i0 = 0; i0 < nStrLoc; i0++)
  {
    it1 = neighborMap.find(strLocDeque[i0].refIdStr);

    if (it1 != neighborMap.end())
    {
      it1->second.first.push_back(strLocDeque[i0].refBeg);
      it1->second.second.push_back(strLocDeque[i0].refEnd);
    }
  }

  // Resize vectors and sort elements in them.
  for (it1 = neighborMap.begin(); it1 != neighborMap.end(); it1++)
  {
    it1->second.first.resize(it1->second.first.size());
    it1->second.second.resize(it1->second.second.size());

    sort(it1->second.first.begin(), it1->second.first.end());
    sort(it1->second.second.begin(), it1->second.second.end());
  }

  return 0;
}

int StrUtil::GetAnchorLen()
{
  return this->anchorLen;
}

int StrUtil::GetLibPos(const Read &r, const map<string, int> &idLibPosMap)
{
  map<string, int>::const_iterator it0;

  it0 = idLibPosMap.find(r.id);

  if (it0 == idLibPosMap.end())
    return -1;

  return it0->second;
}

int StrUtil::GetRLen(const Read &r)
{
  int cigLen = int(r.cigarData.size());
  int rLen = 0;

  for (int i0 = 0; i0 < cigLen; i0++)
  {
    char type = r.cigarData[i0].Type;
    if (type != 'D' && type != 'P')
      rLen += r.cigarData[i0].Length;
  }

  return rLen;
}

bool StrUtil::IsCovered(const Region &read, const Region &strLoc)
{
  if (read.refBeg <= strLoc.refBeg - anchorLen
      && read.refEnd >= strLoc.refEnd + anchorLen)
    return true;
  else
    return false;
}

int StrUtil::QnameToRpPosMaq(Region &reg, const string &qname)
{
  //const string func = "SU::QnameToRpPosMaq ";

  size_t posBeg = 0;
  size_t posEnd = qname.find('_');

  reg.refIdStr = qname.substr(posBeg, posEnd - posBeg);

  posBeg = posEnd + 1;
  posEnd = qname.find('_', posBeg);
  reg.refBeg = lexical_cast<int>(qname.substr(posBeg, posEnd - posBeg));

  posBeg = posEnd + 1;
  posEnd = qname.find('_', posBeg);
  reg.refEnd = lexical_cast<int>(qname.substr(posBeg, posEnd - posBeg));

  return 0;
}

int StrUtil::ReadLenWithClipping(
    const std::vector<CigarOp> &cigarData,
    int readLen)
{
  return readLen;
}

int StrUtil::ReadLenWithoutClipping(
    const vector<CigarOp> &cigarData,
    int readLen)
{
  for (int i0 = 0; i0 < int(cigarData.size()); i0++)
  {
    char base = cigarData[i0].Type;
    int opLen = cigarData[i0].Length;

    if (base == 'H' || base == 'S')
      readLen -= opLen;
  }

  return readLen;
}

int StrUtil::RectifyPosWithClip(Region &reg, const Read &read)
{
  int cigLen = read.cigarData.size();

  reg.refBeg = read.pos;
  reg.refEnd = read.pos + read.length - 1;

  if (cigLen > 0 && read.cigarData[0].Type == 'S')
  {
    reg.refBeg -= read.cigarData[0].Length;
    reg.refEnd -= read.cigarData[0].Length;
  }

  return 0;
}

int StrUtil::RectifyPosWoClip(Region &reg, const Read &read)
{
  int cigLen = read.cigarData.size();
  //int posBeg = read.pos;

  reg.refBeg = read.pos;
  reg.refEnd = read.pos + read.length - 1;

  if (cigLen > 1)
  {
    if (read.cigarData[0].Type == 'S')
      reg.refEnd -= read.cigarData[0].Length;
    if (read.cigarData[cigLen - 1].Type == 'S')
      reg.refEnd -= read.cigarData[cigLen - 1].Length;
  }

  return 0;
}

int StrUtil::RpPosToRPosMaq(
    Region &regRead,
    const Region &regRp,
    int len,
    int mpers)
{
  //const string func = "SU::RpPosToRPosMaq ";

  regRead.refId = regRp.refId;
  regRead.refIdStr = regRp.refIdStr;

  //cout << func << regRp.refBeg << " " << regRp.refEnd << " " << len << " " << mpers << " ";

  if (mpers < 0)
  {
    regRead.refBeg = regRp.refEnd - len + 1;
    regRead.refEnd = regRp.refEnd;
  }
  else
  {
    regRead.refBeg = regRp.refBeg;
    regRead.refEnd = regRead.refBeg + len - 1;
  }

  return 0;
}

int StrUtil::SetAnchorLen(int anchorLenNew)
{
  if (anchorLenNew < 1)
  {
    cout << "Anchor length must be at least 1." << endl;
    return 1;
  }

  this->anchorLen = anchorLenNew;

  return 0;
}

#endif /* STRUTIL_HPP_ */
