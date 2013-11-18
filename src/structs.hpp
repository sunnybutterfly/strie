#ifndef STRUCTS_HPP_
#define STRUCTS_HPP_

#include <list>

#include <boost/multi_array.hpp>
#include <api/BamAux.h>

// Namespaces
using namespace std;
using namespace BamTools;

// Structures
struct Read
{
    string  qname;
    string  id;
    int     refId;
    int     length;
    int     insertSize;
    int     pos;
    int     mpos;
    int     mapq;
    int     flag;
    int     libPosNo;
    std::vector<CigarOp> cigarData;
    int     overlapType;
    int     overlapIndelSizeDs;
    int     overlapIndelSize;
    int     overlapIndelSizeUs;
};



struct ReadPair
{
    Read r1;
    Read r2;
};



struct IndelPoint
{
    int i1;
    int i2;
};



struct IndelInterval
{
    int i1;
    int i2;
    int i1Low;
    int i1High;
    int i2Low;
    int i2High;
};



struct Region
{
    int refId;
    string refIdStr;
    int refBeg;
    int refEnd;
};



struct LibInsertSize
{
    string libTag;
    int insertSizeMin;
    int insertSizeMax;
};



enum CiEstimator
{
    none, Credible_Interval, Percentile_PDF
};



// To be used with "gettimeofday()" in library sys/time.h
struct TimeDur
{
    timeval beg;
    timeval end;
};




#endif /* STRUCTS_HPP_ */
